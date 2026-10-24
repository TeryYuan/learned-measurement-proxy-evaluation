"""Exact legacy context checks and 27 permitted inner-inner proxy fits."""
from b1_common import *
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestRegressor
PARAMS=dict(n_estimators=300,criterion='squared_error',min_samples_leaf=3,max_features=.3,bootstrap=True,max_depth=None,random_state=4109,n_jobs=4)
class ProxyModel:
    def __init__(self,allowed,context):self.allowed=set(allowed);self.context=context
    def fit(self,ids,x,z):
        guard(ids,self.allowed);assert len(ids)==len(set(ids));self.train_ids=list(ids);self.scaler=StandardScaler().fit(z);self.constant=np.std(z,axis=0)==0;self.constant_values=np.mean(z,axis=0)
        self.model=RandomForestRegressor(**PARAMS).fit(x,self.scaler.transform(z));return self
    def predict(self,ids,x):
        guard(ids,self.allowed)
        if set(ids)&set(self.train_ids):raise PermissionError('proxy training/query identity overlap')
        p=self.scaler.inverse_transform(self.model.predict(x));p[:,self.constant]=self.constant_values[self.constant];return p
class ProxyCache:
    def __init__(self,m,z,x):
        self.allowed=set(m.parent_id);self.z=z;self.x=x;self.contexts=[];self.manifest=[];self.cache={}
        self.boundaries=pd.read_csv(OLD/'pilot/training_query_boundaries.csv').set_index('context')
        self.targets=pd.read_csv(OLD/'pilot/proxy_training_target_support.csv')
        pr=pd.read_csv(OLD/'pilot/proxy_unique_parent_OOF.csv');guard(pr.parent_id,self.allowed)
        self.outer=pr.pivot(index='parent_id',columns='geneset',values='proxy').loc[:,z.columns]
        self.inner={f:pd.read_csv(OLD/f'pilot/outer{f}_training_proxy_crossfit.csv').set_index('parent_id').loc[:,z.columns] for f in range(3)}
    def reuse(self,context,S,Q,path,cached):
        guard(S+Q,self.allowed);assert not set(S)&set(Q)
        obj=load_legacy_model(path);b=self.boundaries.loc[context]
        assert obj.train_ids==S==b.training_ids.split(';') and Q==b.query_ids.split(';')
        assert obj.model.n_features_in_==self.x.shape[1]
        assert all(obj.model.get_params()[k]==v for k,v in PARAMS.items())
        target=self.z.loc[S].to_numpy();sc=StandardScaler().fit(target)
        assert np.allclose(obj.scaler.mean_,sc.mean_,atol=1e-13,rtol=0) and np.allclose(obj.scaler.var_,sc.var_,atol=1e-13,rtol=0)
        saved=self.targets.loc[self.targets.context.eq(context)].set_index('geneset').loc[self.z.columns]
        assert np.allclose(saved.training_mean,sc.mean_,atol=1e-13,rtol=0) and np.allclose(saved.training_variance,sc.var_,atol=1e-13,rtol=0)
        actual=obj.scaler.inverse_transform(obj.model.predict(self.x.loc[Q].to_numpy()));actual[:,obj.constant]=obj.constant_values[obj.constant]
        delta=float(abs(actual-cached.loc[Q].to_numpy()).max());assert delta<1e-12
        f=int(context[5]);ctx=json.loads((OLD/f'pilot/contexts/outer{f}_input_context.json').read_text(encoding='utf8'))
        assert ctx['proxy_input_columns']==self.x.columns.tolist() and ctx['geneset_order']==self.z.columns.tolist()
        row=dict(context=context,status='exact_legacy_reuse',training_ids=';'.join(S),query_ids=';'.join(Q),training_parents=len(S),query_parents=len(Q),outer_validation_overlap=0,object_path=str(path.relative_to(ROOT)),cached_prediction_max_error=delta,checks='training ID order; original input NPZ; feature order; targets; scaler mean/variance; RF config; direct object prediction cache equivalence')
        self.contexts.append(row);log('legacy_proxy_reuse',**row);return cached.loc[Q].copy()
    def outer_pair(self,f,S,Q):
        zq=self.reuse(f'outer{f}_proxy',S,Q,OLD/f'pilot/contexts/outer{f}_proxy.joblib',self.outer)
        assignment=folds(S)
        for j in range(3):
            U=[p for p in S if assignment[p]!=j];V=[p for p in S if assignment[p]==j]
            self.reuse(f'outer{f}_inner{j}',U,V,OLD/f'pilot/contexts/outer{f}_inner{j}_proxy.joblib',self.inner[f])
        return self.inner[f].loc[S].copy(),zq
    def inner_pair(self,f,j,U,V):
        # V predictions from saved g_U; U is never sliced from the old S-level crossfit.
        zout=pd.DataFrame(index=U,columns=self.z.columns,dtype=float);assignment=folds(U)
        for k in range(3):
            I=[p for p in U if assignment[p]!=k];J=[p for p in U if assignment[p]==k];key=f'outer{f}_inner{j}_crossfit{k}'
            assert not (set(I)|set(J))&set(V);assert set(I)|set(J)==set(U)
            p=E/f'models/{key}_proxy.joblib';predp=E/f'modeling/{key}_proxy_predictions.csv'
            if p.exists() and predp.exists():
                obj=joblib.load(p);assert obj.train_ids==I
                pred=pd.read_csv(predp).set_index('parent_id').loc[J,self.z.columns].to_numpy();status='valid_new_cache_reuse'
            else:
                start=time.time();obj=ProxyModel(self.allowed,key).fit(I,self.x.loc[I].to_numpy(),self.z.loc[I].to_numpy());pred=obj.predict(J,self.x.loc[J].to_numpy());joblib.dump(obj,p,compress=3)
                csv(pd.DataFrame(pred,index=J,columns=self.z.columns).reset_index(names='parent_id'),predp);status='new_fit_completed';log('proxy_RF_fit',context=key,training_ids=I,query_ids=J,params=PARAMS,elapsed_seconds=time.time()-start)
            zout.loc[J]=pred
            self.contexts.append(dict(context=key,status=status,training_ids=';'.join(I),query_ids=';'.join(J),training_parents=len(I),query_parents=len(J),outer_validation_overlap=0,inner_validation_overlap=0,object_path=str(p.relative_to(ROOT)),cached_prediction_max_error=0,checks='I/J within current U; excludes V and Q; current I target scaling'))
            self.manifest.append(dict(context=key,type='inner_inner_RF300',path=str(p.relative_to(E)),training_ids=';'.join(I),query_ids=';'.join(J),production_model=False))
        assert np.isfinite(zout.to_numpy()).all()
        csv(zout.reset_index(names='parent_id'),E/f'modeling/outer{f}_inner{j}_training_proxy_crossfit.csv')
        return zout,self.inner[f].loc[V].copy()
    def save(self):
        csv(self.contexts,E/'modeling/proxy_training_contexts.csv');csv(self.manifest,E/'modeling/new_proxy_model_manifest.csv')
