"""Phase 1 primitives. Model access is restricted to verified A and V1-labelled O."""
import functools,json,pathlib,time,datetime,warnings,math
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors,rdFingerprintGenerator
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.ensemble import RandomForestClassifier,RandomForestRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import average_precision_score,roc_auc_score,mean_absolute_error,mean_squared_error,r2_score
from dili_phase0.core import DESCRIPTORS
from dili_phase1_log import append_jsonl

ROOT=pathlib.Path(__file__).resolve().parents[1]
DATA=ROOT/'data/phase1';LOG=ROOT/'logs/phase1'
TASKS=['BSEPi','MTX_MP','ARE'];MODELS=['M0','M1','Mperm','Mchem']
FP=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)

def log_event(kind,**kw):
    LOG.mkdir(exist_ok=True,parents=True)
    append_jsonl(LOG/'events.jsonl',{'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'event':kind,**kw})

def structural_features(smiles):
    m=Chem.MolFromSmiles(smiles)
    if m is None:raise ValueError('invalid two-dimensional input')
    if any(a.GetChiralTag()!=Chem.ChiralType.CHI_UNSPECIFIED for a in m.GetAtoms()):raise ValueError('main_2d contains assigned stereo')
    b=np.asarray(FP.GetFingerprintAsNumPy(m),dtype=np.uint8)
    d=np.array([getattr(Descriptors,k)(m) for k in DESCRIPTORS],dtype=np.float64)
    if not np.isfinite(d).all():raise ValueError('nonfinite descriptor')
    return b,d

def require_role(role):
    if role not in ['O']+['A_'+t for t in TASKS]:raise PermissionError('Phase1 allows A/BSEPi,MTX_MP,ARE and O only; N and dedicated C remain closed')

@functools.lru_cache(maxsize=4)
def load_dataset(role):
    require_role(role)
    ix=pd.read_csv(DATA/f'inputs/{role}_index.csv')
    with np.load(DATA/f'inputs/{role}_base.npz',allow_pickle=False) as f:b=f['b'].copy();d=f['d'].copy();ids=f['ids'].copy()
    assert np.array_equal(ids,ix.member_id.astype(str)) and not ix.member_id.duplicated().any()
    return ix,b,d

def validate_access(role,ids):
    require_role(role);ids=np.asarray(ids,dtype=str)
    allowed=load_dataset(role)[0]
    if len(set(ids))!=len(ids) or not set(ids)<=set(allowed.member_id):raise PermissionError('Training/prediction IDs outside authorized role')
    return allowed.set_index('member_id').loc[ids]

def fit_forest(X,y,ids,role,key,params,branch='real'):
    ix=validate_access(role,ids);X=np.asarray(X);y=np.asarray(y)
    assert len(X)==len(y)==len(ix) and np.isfinite(X).all()
    if role=='O':
        assert X.shape[1] in {2056,2059} and np.array_equal(y,ix.y_v1)
    else:
        assert X.shape[1]==2048
        if branch=='real':assert np.array_equal(y,ix.y)
        elif branch=='perm':assert np.array_equal(np.sort(y),np.sort(ix.y))
        elif branch=='chem':
            _,_,des=load_dataset(role);target={'BSEPi':'MolLogP','MTX_MP':'MolWt','ARE':'TPSA'}[role[2:]]
            rows=ix.feature_row.to_numpy(int);assert np.allclose(y,des[rows,DESCRIPTORS.index(target)])
        else:raise ValueError('unknown auxiliary branch')
    log_event('fit_start',key=key,role=role,branch=branch,ids=list(map(str,ids)),n=len(y),dimension=X.shape[1],params=params)
    model=RandomForestRegressor(**params) if branch=='chem' else RandomForestClassifier(**params)
    started=time.perf_counter()
    try:model.fit(X,y)
    except Exception as e:
        log_event('fit_failed',key=key,role=role,error=repr(e));raise
    log_event('fit_complete',key=key,role=role,branch=branch,n=len(y),seconds=time.perf_counter()-started,trees=len(model.estimators_))
    return model

def predict_allowed(model,X,ids,role,key,regression=False):
    validate_access(role,ids)
    assert len(X)==len(ids)
    pred=model.predict(X) if regression else model.predict_proba(X)[:,1]
    log_event('predict_complete',key=key,role=role,ids=list(map(str,ids)),n=len(ids),regression=regression)
    return pred

def permutation_seed(task_index,fold):return 2718+1000*(task_index+1)+(fold if fold>=0 else 100)
def permute_training(y,indices,seed):
    local=np.asarray(y)[np.asarray(indices,dtype=int)].copy()
    if not np.isfinite(local).all():raise ValueError('Permutation requires selected observed training cells')
    return np.random.default_rng(seed).permutation(local)

def assign_folds(y,groups,k,minimum_groups=1):
    y=np.asarray(y);groups=np.asarray(groups);result=np.full(len(y),-1,int)
    if len(np.unique(y))!=2:raise ValueError('requires both classes')
    for f,(tr,va) in enumerate(StratifiedGroupKFold(n_splits=k,shuffle=True,random_state=4109).split(np.zeros(len(y)),y,groups)):
        if set(groups[tr]) & set(groups[va]):raise AssertionError('group overlap')
        if len(np.unique(y[tr]))!=2 or len(np.unique(y[va]))!=2:raise ValueError('fold lacks a class')
        if min(len(set(groups[tr])),len(set(groups[va])))<minimum_groups:raise ValueError('fold has insufficient distinct groups')
        result[va]=f
    assert (result>=0).all();return result

def split_support(d,prefix):
    sizes=d.groupby('fold_group').size();y=d.y if 'y' in d else d.y_v1
    return {prefix+'n':len(d),prefix+'positive':int(y.sum()),prefix+'negative':int(len(d)-y.sum()),prefix+'groups':len(sizes),
      prefix+'positive_bearing_groups':d.loc[y.eq(1),'fold_group'].nunique(),prefix+'negative_bearing_groups':d.loc[y.eq(0),'fold_group'].nunique(),
      prefix+'largest_group_n':int(sizes.max()),prefix+'largest_group_fraction':float(sizes.max()/len(d))}

def choose_grid(rows):
    best=max(r['mean_AP'] for r in rows)
    tied=[r for r in rows if abs(r['mean_AP']-best)<=1e-12]
    return sorted(tied,key=lambda r:(-r['min_samples_leaf'],0 if r['max_features']=='sqrt' else 1))[0]

def confusion(y,pred):
    y=np.asarray(y);pred=np.asarray(pred,dtype=int)
    tp=int(((y==1)&(pred==1)).sum());tn=int(((y==0)&(pred==0)).sum())
    fp=int(((y==0)&(pred==1)).sum());fn=int(((y==1)&(pred==0)).sum());reasons={}
    out={'TP':tp,'TN':tn,'FP':fp,'FN':fn}
    for name,num,den in [('MCC',tp*tn-fp*fn,math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))),('sensitivity',tp,tp+fn),('specificity',tn,tn+fp),('PPV',tp,tp+fp),('NPV',tn,tn+fn)]:
        out[name]=num/den if den else None
        if not den:reasons[name]='zero denominator'
    out['undefined_reasons']=json.dumps(reasons);return out

def select_threshold(y,r):
    y=np.asarray(y);r=np.asarray(r);rows=[]
    for t in np.r_[np.unique(r),np.nextafter(r.max(),np.inf)]:
        m=confusion(y,r>=t)['MCC']
        if m is not None:rows.append((float(t),m))
    if not rows:return .5,{'fallback':True,'reason':'all MCC candidates undefined','candidate_count':0}
    best=max(m for t,m in rows);threshold=max(t for t,m in rows if abs(m-best)<=1e-12)
    return threshold,{'fallback':False,'training_internal_MCC':best,'candidate_count':len(rows)}

def calibrate(r,y,ids,key):
    validate_access('O',ids)
    start=time.perf_counter();info={'method':'sigmoid a*raw+b','solver':'lbfgs','regularization':'none C=inf','class_weight':None,'max_iter':2000,'tol':1e-9,'constant_raw':bool(np.ptp(r)==0)}
    try:
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter('always')
            model=LogisticRegression(C=np.inf,l1_ratio=0.0,solver='lbfgs',max_iter=2000,tol=1e-9,class_weight=None).fit(np.asarray(r).reshape(-1,1),y)
        for w in ws:log_event('calibration_warning',key=key,message=str(w.message),category=w.category.__name__)
        converged=not any(issubclass(w.category,ConvergenceWarning) for w in ws)
        info.update(converged=converged,a=float(model.coef_[0,0]),b=float(model.intercept_[0]),n_iter=int(model.n_iter_[0]),warnings=[str(w.message) for w in ws])
        if not converged:model=None
    except Exception as e:model=None;info.update(converged=False,error=repr(e))
    log_event('calibration_fit',key=key,role='O',ids=list(map(str,ids)),seconds=time.perf_counter()-start,**info)
    return model,info

def classification_metrics(y,r,ids,calibrated=None,pred=None):
    y=np.asarray(y,dtype=int);r=np.asarray(r);n=len(y)
    out={'n':n,'positive':int(y.sum()),'negative':int(n-y.sum()),'positive_fraction':float(y.mean()),
      'AP':float(average_precision_score(y,r)) if len(np.unique(y))==2 else None,'ROC_AUC':float(roc_auc_score(y,r)) if len(np.unique(y))==2 else None,
      'Brier_raw':float(np.mean((y-r)**2))}
    if calibrated is not None:out['Brier_calibrated']=float(np.mean((y-calibrated)**2)) if np.isfinite(calibrated).all() else None
    if pred is not None:out.update(confusion(y,pred))
    k=max(1,math.ceil(.1*n));ix=np.lexsort((np.asarray(ids,dtype=str),-r))[:k];hits=int(y[ix].sum())
    out.update(top_k=k,top_hits=hits,precision_at_q=hits/k,recall_at_q=hits/int(y.sum()) if y.sum() else None,EF_at_q=(hits/k)/float(y.mean()) if y.sum() else None)
    return out

def regression_metrics(y,r):
    constant=bool(np.ptp(y)==0)
    return {'n':len(y),'MAE':float(mean_absolute_error(y,r)),'RMSE':float(math.sqrt(mean_squared_error(y,r))),
      'R2':None if constant or len(y)<2 else float(r2_score(y,r)),'constant_target':constant,'target_mean':float(np.mean(y)),
      'undefined_reason':'constant target or n<2' if constant or len(y)<2 else ''}
