"""Authorized finite nested grouped development; exact checkpoints; <=4 model threads."""
from common7 import *
import argparse,joblib,warnings,math
from sklearn.metrics import average_precision_score
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import ConvergenceWarning
from prepare7A import threshold_function

def inputs():
 a=np.load(A/'inputs.npz');m=read(A/'members.csv');c=json.loads((A/'fold_contexts.json').read_text('utf8'));sp=json.loads((A/'model_protocol_before_timing.json').read_text('utf8'));return a,m,c,sp
def prob(model,x):
 i=np.flatnonzero(np.asarray(model.classes_)==1);assert len(i)==1
 return model.predict_proba(x)[:,i[0]]
def arrays(a,family,branch):
 x=np.column_stack([a['counts'],a['descriptors']]).astype(float) if family=='B3' else a['descriptors'].astype(float)
 if branch!='base':x=np.column_stack([x,a['real' if branch=='real' else 'legacy_perm']])
 return x
def fit_one(family,params,sp,x,y,smiles,train,query):
 start=time.perf_counter()
 if family=='B3':
  from lightgbm import LGBMClassifier
  model=LGBMClassifier(**sp['B3_fixed'],**params)
  with warnings.catch_warnings(record=True) as ws:
   warnings.simplefilter('always');model.fit(x[train],y[train]);p=prob(model,x[query])
  used=model.booster_.params
  assert used['deterministic'] and used['force_col_wise'] and used['num_threads']==4
  context=dict(actual_booster_params=used,iterations=model.n_iter_,warnings=[str(w.message) for w in ws])
 else:
  from chemprop7 import fit_graph
  model,p,context=fit_graph(smiles,x,y,train,query,params['dropout'])
 assert len(p)==len(query) and np.isfinite(p).all() and ((p>=0)&(p<=1)).all()
 import psutil
 return model,p,dict(context,seconds=time.perf_counter()-start,rss_after_fit_bytes=psutil.Process().memory_info().rss,memory_scope='process snapshot after fit; not a measured peak')

def technical(family):
 a,m,cs,sp=inputs();c=cs[0];tr=np.array(c['tr']);inner=np.array(c['inner']);U=tr[inner!=0];V=tr[inner==0];branch='legacy_perm';grid=sp[family+'_grid'];ci=len(grid)-1
 folder=A/'cache'/f"{family}_r00_f0_{branch}";folder.mkdir(parents=True,exist_ok=True);cache=folder/f'candidate{ci:02d}_inner0.npz';out=A/f'{family}_technical_timing.json'
 if out.exists() and cache.exists():print('technical cache already exists',family);return
 assert remaining()>3600
 event('A_technical','start',family=family,train_ids=m.iloc[U].member_id.tolist(),query_ids=m.iloc[V].member_id.tolist(),outer_queries_excluded=True,candidate=ci)
 model,p,info=fit_one(family,grid[ci],sp,arrays(a,family,branch),a['y'],a['smiles'],U,V)
 np.savez_compressed(cache,train_idx=U,query_idx=V,score=p);joblib.dump(model,folder/'technical_model.joblib',compress=3)
 nominal={'B3':5550,'B4':1050}[family]
 rec=dict(family=family,status='technical_pass',candidate=ci,branch=branch,repeat=0,outer_fold=0,inner_fold=0,train_ids=m.iloc[U].member_id.tolist(),query_ids=m.iloc[V].member_id.tolist(),nominal_fits=nominal,conservative_estimate_seconds=info['seconds']*nominal*1.5,score_not_used_for_admission=True,outer_predictions_generated=0,**info)
 dump(rec,out);event('A_technical','fit_complete',family=family,new_fit=1,query_rows=len(V),cache_path=str(cache),**info);print(json.dumps(rec,ensure_ascii=False,default=str)[:1200])

def run_family(family):
 a,m,cs,sp=inputs();adm=json.loads((A/'frozen_admission.json').read_text('utf8'));assert family in adm['admitted_families'];grid=sp[family+'_grid'];threshold=threshold_function();y=a['y'];smiles=a['smiles'];count=0
 for c in cs:
  tr=np.array(c['tr']);va=np.array(c['va']);inner=np.array(c['inner']);r=c['repeat'];f=c['outer_fold']
  for branch in ['base','real','legacy_perm']:
   key=f'{family}_r{r:02d}_f{f}_{branch}';folder=A/'cache'/key;folder.mkdir(parents=True,exist_ok=True)
   if (folder/'completed.json').exists():continue
   x=arrays(a,family,branch);context=dict(family=family,branch=branch,repeat=r,outer_fold=f,train_ids=m.iloc[tr].member_id.tolist(),query_ids=m.iloc[va].member_id.tolist(),training_y=y[tr].tolist(),inner_folds=inner.tolist(),input_dim=x.shape[1],grid=grid,seed=4109,post_hoc=True)
   dump(context,folder/'context.json');pis=[];scores=[]
   for ci,params in enumerate(grid):
    pi=np.full(len(tr),np.nan);aps=[]
    for j in range(3):
     U=tr[inner!=j];V=tr[inner==j];cache=folder/f'candidate{ci:02d}_inner{j}.npz'
     if cache.exists():
      z=np.load(cache);assert np.array_equal(z['train_idx'],U) and np.array_equal(z['query_idx'],V);p=z['score'];event('A_train','cache_reused',key=key,candidate=ci,inner=j)
     else:
      if remaining()<3600:
       dump(dict(status='COMPUTATION_INCOMPLETE_BUDGET',next_key=key,candidate=ci,inner=j),A/'checkpoint.json');return False
      event('A_train','fit_started',key=key,candidate=ci,inner=j)
      model,p,info=fit_one(family,params,sp,x,y,smiles,U,V);np.savez_compressed(cache,train_idx=U,query_idx=V,score=p)
      event('A_train','fit_complete',key=key,candidate=ci,inner=j,new_fit=1,query_rows=len(V),**info);del model
     pi[inner==j]=p;aps.append(float(average_precision_score(y[V],p)))
    assert np.isfinite(pi).all();pis.append(pi);scores.append(dict(candidate=ci,params=params,mean_inner_AP=float(np.mean(aps)),inner_AP=aps))
   best=max(s['mean_inner_AP'] for s in scores);chosen=next(s for s in scores if best-s['mean_inner_AP']<=1e-12);ci=chosen['candidate'];pi=pis[ci]
   dump(scores,folder/'candidate_scores.json');np.savez_compressed(folder/'all_candidate_inner_oof.npz',score=np.column_stack(pis),member_ids=m.iloc[tr].member_id.to_numpy(str),inner_fold=inner,y=y[tr])
   save(dict(member_id=m.iloc[tr].member_id,y=y[tr],raw_score=pi,inner_fold=inner),folder/'selected_inner_oof.csv')
   if remaining()<3600:dump(dict(status='COMPUTATION_INCOMPLETE_BUDGET',next_key=key,step='outer_refit'),A/'checkpoint.json');return False
   model,p,info=fit_one(family,grid[ci],sp,x,y,smiles,tr,va);event('A_train','fit_complete',key=key,candidate=ci,inner='outer_refit',new_fit=1,query_rows=len(va),**info)
   th,tinfo=threshold(y[tr],pi)
   with warnings.catch_warnings(record=True) as ws:
    warnings.simplefilter('always');cal=LogisticRegression(C=np.inf,l1_ratio=0.,solver='lbfgs',class_weight=None,max_iter=2000,tol=1e-9).fit(pi.reshape(-1,1),y[tr])
   ok=not any(issubclass(w.category,ConvergenceWarning) for w in ws);calp=prob(cal,p.reshape(-1,1)) if ok else np.full(len(va),np.nan)
   calibration=dict(a=float(cal.coef_[0,0]),b=float(cal.intercept_[0]),converged=ok,n_iter=cal.n_iter_.tolist(),warnings=[str(w.message) for w in ws]);event('A_train','calibration_fit',key=key,new_fit=1,**calibration)
   models=A/'models';models.mkdir(exist_ok=True);joblib.dump(dict(model=model,calibrator=cal if ok else None,threshold=th,context=context,selected_candidate=ci,params=grid[ci],fit_info=info),models/f'{key}.joblib',compress=3)
   oof=m.iloc[va][['member_id','fold_group','y']].copy();oof['raw_score']=p;oof['calibrated_score']=calp;oof['binary_pred']=(p>=th).astype(int);oof['threshold']=th;oof['repeat']=r;oof['outer_fold']=f;oof['family']=family;oof['branch']=branch;oof['selected_candidate']=ci
   save(oof,folder/'outer_predictions.csv');dump(dict(status='complete',chosen=chosen,threshold=th,threshold_info=tinfo,calibration=calibration,n_query=len(va),model_path=str(models/f'{key}.joblib')),folder/'completed.json');count+=1
   print(f'{now()} {key} complete ({count} new outer-branch contexts this run)',flush=True)
 done=list((A/'cache').glob(f'{family}_*/completed.json'));assert len(done)==150
 save(pd.concat([read(p.parent/'outer_predictions.csv') for p in done]),A/f'{family}_OOF.csv');dump(dict(status='complete',outer_branch_contexts=150,nominal_model_fits=150*(3*len(grid)+1),new_auxiliary_fits=0),A/f'{family}_completion.json');return True
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('action',choices=['technical','run']);ap.add_argument('family',choices=['B3','B4']);args=ap.parse_args()
 try:technical(args.family) if args.action=='technical' else run_family(args.family)
 except Exception as e:
  event('A','failure',action=args.action,family=args.family,error=repr(e),trace=traceback.format_exc());raise
