"""Finite Phase 2BC primitives. No N model query and no method selection."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='1'
import pathlib,json,datetime,time,contextlib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier,RandomForestRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import average_precision_score,roc_auc_score
from dili_phase1_log import append_jsonl
from dili_phase1 import structural_features,confusion,TASKS
from dili_phase2a import positive_probability,bootstrap_indices,percentile_interval

ROOT=pathlib.Path(__file__).resolve().parents[1]
LOG=ROOT/'logs/phase2bc'; RC=ROOT/'reports/phase2b_C'; RE=ROOT/'reports/phase2c_BSEP'
RN=ROOT/'reports/phase2a_N_posthoc'; DC=ROOT/'data/phase2b_C'; DE=ROOT/'data/phase2c_BSEP'
MODELS=['C_M0','C_M1','C_Mperm','C_Mchem']
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def read(path):return pd.read_csv(ROOT/path,float_precision='round_trip')
def js(path):return json.loads((ROOT/path).read_text(encoding='utf-8-sig'))
def dump(obj,path):
    p=ROOT/path;p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,ensure_ascii=False,indent=2,default=lambda x:x.item() if hasattr(x,'item') else str(x)),encoding='utf8')
def save(df,path):
    p=ROOT/path;p.parent.mkdir(parents=True,exist_ok=True);df.to_csv(p,index=False,encoding='utf-8-sig')
def event(kind,**kw):
    LOG.mkdir(parents=True,exist_ok=True);append_jsonl(LOG/'events.jsonl',{'utc':utc(),'event':kind,**kw})
def metric(y,p):
    y=np.asarray(y);p=np.asarray(p);n=len(y);both=len(np.unique(y))==2
    return dict(n=n,positive=int(y.sum()),negative=int(n-y.sum()),positive_fraction=float(y.mean()) if n else None,
        AP=float(average_precision_score(y,p)) if both else None,ROC_AUC=float(roc_auc_score(y,p)) if both else None,
        Brier_raw=float(np.mean((y-p)**2)) if n else None,rank_metric_NA_reason='' if both else 'single_class_or_empty')
def brier_contributions(y1,y2,pold,pnew):
    y1=np.asarray(y1);y2=np.asarray(y2);pold=np.asarray(pold);pnew=np.asarray(pnew);n=len(y1)
    return (y2-y1)*(1-2*pold)/n,((pnew-y2)**2-(pold-y2)**2)/n
def c_arrays():
    m=read('data/phase2b_C/C_membership.csv')
    with np.load(DC/'C_fixed_inputs.npz',allow_pickle=False) as f:
        assert f['ids'].tolist()==m.member_id.tolist()
        x={name:f[name].copy() for name in MODELS}
    return m,x
def validate_c_fit(m,train_ids,query_ids,version,model,params):
    if model not in MODELS or version not in ['v1','v2']:raise PermissionError('Unapproved C strategy/version')
    if params!=js('config/phase2bc_execution.json')['C']['rf']:raise PermissionError('C parameters changed')
    tr=list(train_ids);qu=list(query_ids)
    if len(set(tr))!=len(tr) or len(set(qu))!=len(qu):raise PermissionError('Repeated identity')
    if not set(tr)<=set(m.member_id) or not set(qu)<=set(m.member_id):raise PermissionError('Non-C identity')
    iq=m.set_index('member_id').loc[qu];folds=iq.fold.unique()
    if len(folds)!=1:raise PermissionError('Exactly one inherited C query fold required')
    expected_train=m.loc[m.fold.ne(folds[0]),'member_id'].tolist()
    expected_query=m.loc[m.fold.eq(folds[0]),'member_id'].tolist()
    if tr!=expected_train or qu!=expected_query:raise PermissionError('Must use exact ordered C complement and fold')
    it=m.set_index('member_id').loc[tr]
    if set(it.fold_group)&set(iq.fold_group):raise PermissionError('Structure group leakage')
    return int(folds[0])
@contextlib.contextmanager
def no_fit(scope,record=True):
    """Reject learning calls for saved-N arithmetic and fixed SI004 inference."""
    patches=[];attempts=[]
    def blocked(*args,**kw):
        attempts.append('fit');raise PermissionError('No learning authorized in '+scope)
    for cls in [RandomForestClassifier,RandomForestRegressor,LogisticRegression,StandardScaler,Pipeline]:
        for name in ['fit','fit_transform','partial_fit']:
            if hasattr(cls,name):
                inherited=name not in cls.__dict__;prior=getattr(cls,name)
                patches.append((cls,name,prior,inherited));setattr(cls,name,blocked)
    try:yield attempts
    finally:
        for cls,name,prior,inherited in reversed(patches):
            if inherited:delattr(cls,name)
            else:setattr(cls,name,prior)
        if record:event('no_fit_guard_complete',scope=scope,blocked_attempts=len(attempts))
def deduplicate_si(rows):
    """Only observed, in-domain records; any discordant observed binary labels isolate identity."""
    groups=[]
    for key,g in rows.groupby('identity_group',sort=True,dropna=False):
        observed=g.loc[g.BSEPi.isin([0,1])]
        labels=observed.BSEPi.unique();valid=g.eligible_2d.all() and isinstance(key,str) and bool(key)
        reason='out_of_domain' if not valid else 'label_conflict' if len(labels)>1 else 'missing_or_invalid_label' if len(labels)!=1 else ''
        groups.append(dict(identity_group=key,source_excel_rows='|'.join(g.excel_row.astype(str)),source_names='|'.join(g.Name.astype(str)),
            raw_n=len(g),observed_n=len(observed),y_BSEPi=int(labels[0]) if len(labels)==1 else None,eligible_base=not reason,base_exclusion_reason=reason))
    return pd.DataFrame(groups)
