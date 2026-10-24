from revision_common import *
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
GRID=[.001,.1,10.,'DROP'];FORCED=list(itertools.product(GRID,GRID[:3]));OPTIONAL=list(itertools.product(GRID,GRID))
ROUTES={'C': [('DROP','DROP')],'X':[(a,'DROP') for a in GRID],'M_optional':OPTIONAL,'M_forced':FORCED,'P_forced':FORCED,'Perm_forced':FORCED}
def branch(route):return 'pred' if route=='P_forced' else 'perm' if route=='Perm_forced' else 'meas'
def run(repeat,outer):
    dest=E/f'alt/r{repeat}_o{outer}';dest.mkdir(exist_ok=True)
    if (dest/'completed.json').exists():return json.loads((dest/'completed.json').read_text())
    start=time.perf_counter();m,t,z,x=alt_data();allowed=set(m.parent_id);rs=2026101301+repeat
    af=read(E/'alt/outer_folds.csv');af=af.loc[af.repeat.eq(repeat)].set_index('parent_id').outer
    T=sorted(af.index[af.ne(outer)]);Q=sorted(af.index[af.eq(outer)]);ia=split(T,3,rs,outer,1)
    ts=t.loc[t.parent_id.isin(T)];tq=t.loc[t.parent_id.isin(Q)]
    contexts=[];perms=[];supports=[];selections=[];scores=[];preds=[];inner_preds=[];proxy_new=0
    def proxy(S,V,key):
        nonlocal proxy_new
        assert not set(S)&set(V) and set(S+V)<=allowed
        p=dest/f'proxy_{key}.joblib';qp=dest/f'proxy_{key}.csv'
        if p.exists() and qp.exists():
            obj=joblib.load(p);assert obj['train_ids']==S and obj['query_ids']==V and obj['columns']==x.columns.tolist() and obj['genesets']==z.columns.tolist()
            assert np.array_equal(obj['target'],z.loc[S].to_numpy()) and np.array_equal(obj['input'],x.loc[S].to_numpy())
            out=read(qp).set_index('parent_id').loc[V,z.columns]
        else:
            tick=time.perf_counter();sc=StandardScaler().fit(z.loc[S]);constant=np.std(z.loc[S].to_numpy(),axis=0)==0
            rf=RandomForestRegressor(n_estimators=300,min_samples_leaf=3,max_features=.3,random_state=4109,n_jobs=4).fit(x.loc[S].to_numpy(),sc.transform(z.loc[S]))
            pz=sc.inverse_transform(rf.predict(x.loc[V].to_numpy()));pz[:,constant]=np.mean(z.loc[S].to_numpy(),axis=0)[constant]
            obj=dict(model=rf,scaler=sc,constant=constant,train_ids=S,query_ids=V,columns=x.columns.tolist(),genesets=z.columns.tolist(),target=z.loc[S].to_numpy(),input=x.loc[S].to_numpy())
            joblib.dump(obj,p,compress=3);out=pd.DataFrame(pz,index=V,columns=z.columns);csv(out.reset_index(names='parent_id'),qp);proxy_new+=1
            event('ALT_proxy_fit',repeat=repeat,outer=outer,key=key,train_ids=S,query_ids=V,seconds=time.perf_counter()-tick)
        contexts.append(dict(key=key,train_ids=';'.join(S),query_ids=';'.join(V),n_train=len(S),n_query=len(V),overlap=0))
        return out
    def perm(S,j):
        src=np.array(S);ps=seed(rs,outer,j,2718);np.random.Generator(np.random.PCG64(ps)).shuffle(src)
        perms.extend(dict(context=j,seed=ps,parent_id=p,assigned_z_parent=v) for p,v in zip(S,src))
        return pd.DataFrame(z.loc[src].to_numpy(),index=S,columns=z.columns)
    def frame(tv,p,route,**kw):
        a=tv[['condition_id','parent_id','source_name','split','y']].copy();a['repeat']=repeat;a['outer']=outer;a['route_model']=route;a['prediction']=p
        for k,v in kw.items():a[k]=v
        return a
    ztcf=pd.DataFrame(index=T,columns=z.columns,dtype=float);zq=proxy(T,Q,'outer')
    pooled={(b,p):[] for b in ['meas','pred','perm'] for p in (OPTIONAL if b=='meas' else FORCED)}
    for j in range(3):
        U=[p for p in T if ia[p]!=j];V=[p for p in T if ia[p]==j];tu=ts.loc[ts.parent_id.isin(U)];tv=ts.loc[ts.parent_id.isin(V)]
        zv=proxy(U,V,f'inner{j}');ztcf.loc[V]=zv
        deep=split(U,3,rs,outer,j,2);zu=pd.DataFrame(index=U,columns=z.columns,dtype=float)
        for k in range(3):
            I=[p for p in U if deep[p]!=k];J=[p for p in U if deep[p]==k];zu.loc[J]=proxy(I,J,f'inner{j}_deep{k}')
        assert np.isfinite(zu).all().all()
        for b,zt,zvinput in [('meas',z,z),('pred',zu,zv),('perm',perm(U,j),z)]:
            q=qweights(tu.parent_id);rr=raw(tu,x,zt);tr=BlockTransform().fit(*rr,q);blocks=tr.transform(*rr);du=DualContext(blocks,tu.y.to_numpy(),q);qc=du.query_cache(tr.transform(*raw(tv,x,zvinput)))
            supports.extend(dict(context=f'inner{j}',branch=b,**v) for v in tr.support(blocks,q))
            for lx,lz in (OPTIONAL if b=='meas' else FORCED):
                pred,_=du.solve(lx,lz,qc);out=frame(tv,pred,b,inner=j,lambdaX=str(lx),lambdaZ=str(lz));pooled[b,(lx,lz)].append(out);inner_preds.append(out)
    chosen={}
    for route,candidates in ROUTES.items():
        sc={}
        for p in candidates:
            d=pd.concat(pooled[branch(route),p]);assert d.condition_id.is_unique and set(d.parent_id)==set(T)
            met=metrics(d.parent_id,d.y,d.prediction);sc[p]=met['MAE_cmp'];scores.append(dict(route_model=route,lambdaX=p[0],lambdaZ=p[1],**met))
        chosen[route]=select(sc)
    for b,r in [('meas','M'),('pred','P'),('perm','Perm')]:
        for lz in GRID[:3]:
            sc={p:metrics((d:=pd.concat(pooled[b,p])).parent_id,d.y,d.prediction)['MAE_cmp'] for p in FORCED if p[1]==lz}
            chosen[f'{r}_fixedZ{lz}']=select(sc)
    q=qweights(ts.parent_id);meanz=z.loc[T].mean(axis=0).to_numpy();modelctx={}
    for b,zt,zvinput in [('meas',z,z),('pred',ztcf,zq),('perm',perm(T,3),z)]:
        rr=raw(ts,x,zt);tr=BlockTransform().fit(*rr,q);blocks=tr.transform(*rr);du=DualContext(blocks,ts.y.to_numpy(),q);qc=du.query_cache(tr.transform(*raw(tq,x,zvinput)));modelctx[b]=(tr,du,qc,zvinput)
        supports.extend(dict(context='refit',branch=b,**v) for v in tr.support(blocks,q))
    reload_errors=[]
    for route,pair in chosen.items():
        b=('perm' if route.startswith('Perm') else 'pred' if route.startswith('P_') else 'meas');tr,du,qc,zvinput=modelctx[b]
        pred,coeff=du.solve(*pair,qc,True);model=BandedModel(tr,du,coeff,T,allowed,*pair,f'r{repeat}_o{outer}_{route}',meanz,x.columns[:-1].tolist(),z.columns.tolist())
        path=dest/f'{route}.joblib';joblib.dump(model,path,compress=3)
        pred2=model.predict(tq.parent_id,*raw(tq,x,zvinput));assert np.allclose(pred,pred2,atol=1e-10,rtol=0)
        if route=='M_forced':
            loaded=joblib.load(path);delta=float(abs(loaded.predict(tq.parent_id,*raw(tq,x,zvinput))-pred).max());assert delta<1e-10;reload_errors.append(delta)
        preds.append(frame(tq,pred,route));zrms=np.nan
        if route=='M_forced' or route.startswith('M_fixedZ'):
            zs=pd.DataFrame(np.tile(meanz,(len(Q),1)),index=Q,columns=z.columns)
            ps=model.predict(tq.parent_id,*raw(tq,x,zq));pm=model.predict(tq.parent_id,*raw(tq,x,zs))
            suffix='' if route=='M_forced' else '_'+route
            preds.extend([frame(tq,ps,'Swap'+suffix),frame(tq,pm,'MeanZ'+suffix)])
            zrms=float(np.sqrt(qweights(tq.parent_id)@(pred-pm)**2))
        selections.append(dict(route_model=route,lambdaC=.01,lambdaX=pair[0],lambdaZ=pair[1],betaZ_norm=float(np.linalg.norm(coeff.get('Z',[]))),query_RMS_measured_minus_meanZ=zrms,n_train=len(T),n_query=len(Q)))
    for r,v in [('zero_response',0.),('training_compound_mean',float(q@ts.y)),('training_weighted_median',median(ts.y.to_numpy(),q))]:preds.append(frame(tq,np.full(len(tq),v),r))
    csv(pd.concat(preds),dest/'predictions.csv');csv(pd.concat(inner_preds),dest/'inner_candidate_predictions.csv');csv(scores,dest/'inner_scores.csv');csv(selections,dest/'selected_parameters.csv');csv(contexts,dest/'proxy_contexts.csv');csv(perms,dest/'permutations.csv');csv(supports,dest/'block_support.csv')
    dump(dict(train_ids=T,query_ids=Q,inner_assignment=ia,proxy_columns=x.columns.tolist(),genesets=z.columns.tolist(),reload_errors=reload_errors),dest/'context.json')
    result=dict(completed_utc=now(),repeat=repeat,outer=outer,elapsed_seconds=time.perf_counter()-start,proxy_contexts=len(contexts),new_proxy_fits=proxy_new,selected_models=len(chosen),inner_unique_ridge_solves=3*(16+12+12))
    assert len(contexts)==13;dump(result,dest/'completed.json');event('ALT_context_complete',**result);print(result,flush=True);return result
if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser();a.add_argument('--repeat',type=int,required=True);a.add_argument('--outer',type=int,required=True);v=a.parse_args();run(v.repeat,v.outer)
