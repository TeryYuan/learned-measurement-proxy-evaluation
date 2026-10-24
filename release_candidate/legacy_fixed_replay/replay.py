"""Recompute existing fixed-score estimates using saved indices; NO RNG/models."""
from pathlib import Path
import json,time
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score,roc_auc_score
from scipy.stats.contingency import odds_ratio
P=Path(__file__).resolve().parent; D=P/'data';O=P/'output';O.mkdir(exist_ok=True)
checks=[];start=time.time()
def read(n):return pd.read_csv(D/n,float_precision='round_trip')
def check(label,value,expected):
    a=np.asarray(value,float);b=np.asarray(expected,float)
    passed=bool(np.allclose(a,b,atol=2e-12,rtol=0,equal_nan=True))
    checks.append(dict(check=label,passed=passed,max_absolute_difference=float(np.nanmax(abs(a-b))) if a.size else 0))
    if not passed:print('MISMATCH',label,value,expected,flush=True)
def met(y,p):return np.array([average_precision_score(y,p),roc_auc_score(y,p),np.mean((y-p)**2)])
def pct(a,coverage=.95):
    tail=(1-coverage)/2;return np.quantile(a,[tail,1-tail],axis=0,method='linear')
print('Reading fixed DILI N scores and original indices',flush=True)
n=read('N_predictions.csv');saved=read('N_metrics_expected.csv').set_index('model')
ix=np.load(D/'N_bootstrap_indices.npz',allow_pickle=False);ids=ix['ids'];idx=ix['indices'];y=ix['labels'];nb={};points={}
for model,g in n.groupby('model',sort=False):
    g=g.set_index('member_id').loc[ids];check('N labels '+model,g.y_v2,y);p=g.raw_score.to_numpy();m=met(y,p);points[model]=m
    check('N point metrics '+model,m,saved.loc[model,['AP','ROC_AUC','Brier_raw']].to_numpy(float))
    check('N calibrated Brier '+model,np.mean((y-g.calibrated_score.to_numpy())**2),saved.loc[model,'Brier_calibrated'])
    # Only the previously designated main AP intervals are recomputed here.
    if model in ['M0','M1','B_star','R_gate']:
        nb[model]=np.array([average_precision_score(y[s],p[s]) for s in idx])
for r in read('N_primary_expected.csv').itertuples():
    a,b=r.comparison.split(' - ');check('N '+r.comparison+' point',points[a][0]-points[b][0],r.delta_AP)
    check('N '+r.comparison+' 95',pct(nb[a]-nb[b]),[r.CI_95_low,r.CI_95_high])
    check('N '+r.comparison+' 97.5',pct(nb[a]-nb[b],.975),[r.CI_97_5_low,r.CI_97_5_high])
print('Replaying E201/E172 historical paired indices',flush=True)
bsep=read('BSEP_predictions.csv');expected=read('BSEP_intervals_expected.csv')
for co in ['E_BSEP','E_no_anyA']:
    a=np.load(D/(co+'_bootstrap_indices.npz'),allow_pickle=False);idx=a['indices'];ids=a['ids'];y=a['labels'];bs={};pts={}
    for model,g in bsep.groupby('model'):
        q=g.set_index('member_id').loc[ids];check(co+' labels '+model,q.y_BSEPi,y);p=q.raw_score.to_numpy()
        bs[model]=np.array([met(y[s],p[s]) for s in idx]);pts[model]=met(y,p)
    bs['real-perm']=bs['g_BSEPi_real']-bs['g_BSEPi_perm'];pts['real-perm']=pts['g_BSEPi_real']-pts['g_BSEPi_perm']
    for r in expected[expected.cohort.eq(co)].itertuples():
        k=['AP','ROC_AUC','Brier_raw'].index(r.metric);check(f'{co} {r.target} {r.metric} point',pts[r.target][k],r.point)
        check(f'{co} {r.target} {r.metric} CI',pct(bs[r.target][:,k],r.coverage),[r.low,r.high])
print('Checking C matched four cells from saved OOF scores',flush=True)
c=read('C_OOF.csv');ce=read('C_cells_expected.csv');ce=ce[ce.scope.eq('pooled_OOF')]
for r in ce.itertuples():
    q=c[c.model.eq(r.model)&c.training_version.eq(r.training_version)]
    assert len(q)==624 and q.member_id.nunique()==624
    check(f'C {r.model} {r.training_version}/{r.scoring_version}',met(q['y_'+r.scoring_version].to_numpy(),q.raw_score.to_numpy()),[r.AP,r.ROC_AUC,r.Brier_raw])
print('Replaying original NAM whole-parent draws',flush=True)
h=read('NAM_predictions.csv');hi=read('NAM_indices.csv').iloc[:,1:].to_numpy(int);ctx=json.loads((D/'NAM_context.json').read_text());parents=ctx['parent_ID_order'];models={};hb={}
target=h[h.route_model.eq('D_cond')].groupby('parent_id').y.agg(['mean',lambda x:np.mean(x*x)]).loc[parents].to_numpy()
def nammetrics(q,indices):
    e=q.assign(a=abs(q.prediction-q.y),s=(q.prediction-q.y)**2).groupby('parent_id')[['a','s']].mean().loc[parents].to_numpy()
    loss=e[indices].mean(axis=1);t=target[indices].mean(axis=1);var=t[:,1]-t[:,0]**2
    r2=1-loss[:,1]/var;r2[var<=16*np.finfo(float).eps*np.maximum(1,t[:,1])]=np.nan
    return np.column_stack([loss[:,0],np.sqrt(loss[:,1]),r2])
for model,q in h.groupby('route_model'):
    models[model]=nammetrics(q,np.arange(24)[None,:])[0];hb[model]=nammetrics(q,hi)
for r in read('NAM_metrics_expected.csv').itertuples():check('NAM '+r.route_model+' points',models[r.route_model],[r.MAE_cmp,r.RMSE_cmp,r.R2_cmp])
for r in read('NAM_single_intervals_expected.csv').itertuples():
    k=['MAE','RMSE','R2'].index(r.metric);check('NAM single '+r.route_model+' '+r.metric,pct(hb[r.route_model][:,k],r.confidence),[r.ci_low,r.ci_high])
for r in read('NAM_intervals_expected.csv').itertuples():
    k={'delta_MAE':0,'delta_RMSE':1}[r.metric];check('NAM '+r.comparison+' '+r.metric+' point',models[r.model_a][k]-models[r.model_b][k],r.estimate)
    check('NAM '+r.comparison+' '+r.metric+' '+str(r.confidence),pct(hb[r.model_a][:,k]-hb[r.model_b][:,k],r.confidence),[r.ci_low,r.ci_high])
pivot=h.pivot(index='condition_id',columns='route_model',values='prediction');identical=['D_cond','D_struct','D_meas_only','D_meas','D_swap','D_meanZ']
for i,a in enumerate(identical):
    for b in identical[i+1:]:
        assert np.array_equal(pivot[a],pivot[b]);check('Exact NAM function identity '+a+'/'+b,pivot[a],pivot[b])
print('Replaying post-hoc bridge original indices and association',flush=True)
bridge=read('bridge_membership.csv');be=read('bridge_intervals_expected.csv');bb={};bp={}
for outcome,filename in [('SI004_BSEP','bridge_indices_BSEP.npz'),('DILI_V2','bridge_indices_DILI.npz')]:
    f=np.load(D/filename,allow_pickle=False);q=bridge.set_index('identity_id').loc[f['identity_id']];y=f['stratification_label'];idx=f['indices'];z={};point={}
    for model,col in [('g_BSEPi_real','BSEP_score_real'),('g_BSEPi_perm','BSEP_score_perm')]:
        p=q[col].to_numpy();z[model]=np.array([met(y[s],p[s])[:2] for s in idx]);point[model]=met(y,p)[:2]
    z['difference']=z['g_BSEPi_real']-z['g_BSEPi_perm'];point['difference']=point['g_BSEPi_real']-point['g_BSEPi_perm']
    for r in be[be.outcome.eq(outcome)&be.metric.isin(['AP','ROC_AUC'])].itertuples():
        name='difference' if r.statistic_type=='paired_difference' else r.model_A;k=['AP','ROC_AUC'].index(r.metric)
        check('Bridge '+outcome+' '+name+' '+r.metric+' point',point[name][k],r.point);check('Bridge '+outcome+' '+name+' '+r.metric+' CI',pct(z[name][:,k]),[r.low,r.high])
    if outcome=='SI004_BSEP':
        dy=q.DILI_y_v2.to_numpy();by=q.y_BSEPi.to_numpy();rd=np.array([dy[s][by[s]==1].mean()-dy[s][by[s]==0].mean() for s in idx])
        r=be[be.metric.eq('RD_BSEP1_minus_BSEP0')].iloc[0];check('Bridge RD interval',pct(rd),[r.low,r.high])
q=bridge[bridge.bridge_eligible];table=np.array([[sum((q.y_BSEPi==b)&(q.DILI_y_v2==y)) for y in [1,0]] for b in [1,0]])
a=read('bridge_association_expected.csv').iloc[0];check('Bridge 2x2',table.ravel(),[a.a,a.b,a.c,a.d]);orr=odds_ratio(table,kind='conditional');ci=orr.confidence_interval()
check('Bridge conditional OR',[orr.statistic,ci.low,ci.high],[a.conditional_OR,a.OR_low,a.OR_high])
pd.DataFrame(checks).to_csv(O/'replay_checks.csv',index=False)
summary={'checks':len(checks),'passed':sum(c['passed'] for c in checks),'failed':sum(not c['passed'] for c in checks),'elapsed_seconds':time.time()-start,'new_random_draws':0,'new_model_queries':0,'comparison_scope':'previously defined contrasts only','absolute_tolerance':2e-12}
(O/'replay_summary.json').write_text(json.dumps(summary,indent=2));print(summary,flush=True)
assert summary['failed']==0
