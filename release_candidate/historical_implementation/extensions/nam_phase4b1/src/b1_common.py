"""Phase4B1 only: legacy inputs read-only; learned operations restricted to dev."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_k]='1'
from pathlib import Path
import sys,json,time,itertools,traceback
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import joblib
E=Path(__file__).resolve().parents[1]
ROOT=E.parents[1]
OLD=ROOT/'extensions/nam_phase4b0'
PACK=ROOT/'NAM_Phase4B1_信号诊断与有限修正_提示词包/NAM_Phase4B1_信号诊断与有限修正'
def now():return datetime.now(timezone.utc).isoformat()
def dump(x,p):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else str(v)),encoding='utf-8')
def csv(d,p):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(d).to_csv(p,index=False,encoding='utf-8-sig')
def log(event,**kw):
    (E/'logs').mkdir(parents=True,exist_ok=True)
    with (E/'logs/events.jsonl').open('a',encoding='utf8') as f:f.write(json.dumps(dict(time=now(),event=event,**kw),ensure_ascii=False,default=str)+'\n')
def guard(ids,allowed):
    if not set(ids)<=set(allowed):raise PermissionError('Holdout/unknown/failed ID prohibited')
def qweights(ids):
    ids=np.asarray(ids);u,n=np.unique(ids,return_counts=True);counts=dict(zip(u,n));w=np.array([1/counts[x] for x in ids]);return w/w.sum()
def metrics(ids,y,p):
    y=np.asarray(y);p=np.asarray(p);q=qweights(ids);den=np.sum(q*(y-q@y)**2)
    return dict(MAE_cmp=float(q@abs(p-y)),RMSE_cmp=float(np.sqrt(q@(p-y)**2)),R2_cmp=float(1-(q@(p-y)**2)/den) if den>0 else np.nan,n_compounds=len(set(ids)),n_conditions=len(y))
def median(y,q):
    order=np.argsort(y,kind='stable');s=np.cumsum(q[order]/q.sum());return float(np.asarray(y)[order[min(np.searchsorted(s,.5-1e-14),len(y)-1)]])
def folds(ids):
    a=np.array(sorted(ids),dtype=str);np.random.Generator(np.random.PCG64(4109)).shuffle(a);return {x:i%3 for i,x in enumerate(a)}
def load_legacy_model(p):
    # Unpickle class definitions only. Never call legacy fit/predict (those write old logs).
    if str(OLD/'src') not in sys.path:sys.path.append(str(OLD/'src'))
    return joblib.load(p)
def data():
    m=pd.read_csv(E/'metadata/dev_members.csv');ids=m.parent_id.tolist()
    t=pd.read_csv(E/'derived/dev_queries.csv');guard(t.parent_id,ids)
    z=pd.read_csv(E/'derived/dev_PHH50.csv').set_index('parent_id')
    a=np.load(E/'derived/dev_structure.npz',allow_pickle=False)
    x=pd.DataFrame(a['X'],index=a['parent_id'],columns=a['columns']);guard(x.index,ids);guard(z.index,ids)
    return m,t,z,x
def old_predict(obj,numeric,route):
    design=np.column_stack([obj.scaler.transform(numeric),obj.encoder.transform(np.asarray(route).reshape(-1,1))])
    return obj.model.predict(design)
def initialize():
    if (E/'config/execution.json').exists():raise RuntimeError('Initialization already saved; do not overwrite pre-fit rules')
    for n in ['config','metadata','derived','diagnostics','modeling','models','reports','figures/source_data','logs','tests']: (E/n).mkdir(parents=True,exist_ok=True)
    spec=json.loads((PACK/'03_phase4b1_spec.json').read_text(encoding='utf-8-sig'))
    spec.update(status='authorized_rules_saved_before_diagnostics_and_fits',activated_utc=now(),permutation_context_seed='2718 + 100*outer_fold + inner_fold; outer refit inner_code=3',inner_inner_split='sorted parent IDs shuffled PCG64(4109), position modulo 3',quantiles='lower observed weighted inverse empirical CDF',constant_columns='exact zero weighted population variance; no query-based rule')
    dump(spec,E/'config/execution.json')
    (E/'config/method_amendment.md').write_text('# Phase4B1 方法修订（拟合前保存）\n\n本轮在观察Phase4B0负向结果后，经作者授权进行一次事后开发。不是原始预登记或原创算法证明。旧结果与对象只读。沿用实际合格NAM_dev、原折、ALT和50集合；留出只读ID排除。\n\n方法设定：母体等权、训练专属中心化/标准化、结构及生物块一阶剂量交互与能量缩放、lambdaC=0.01及两块各0.001/0.1/10/DROP，六路线、两项同对象消融。内部proxy在当前U再交叉拟合；旧g仅经上下文核对后复用。选择按汇总内层母体等权MAE；并列容差1e-10，依次少块、强Z、强X。置换种子2718+100*outer+inner，外重拟合inner=3，唯一母体整向量置换。\n\n目标分布分位数用加权经验CDF达到指定累计概率的较小观测；中位数同规则。全部规则不因结果修改，不增加seed/集合/候选，不拟合全dev，不推断NAM_holdout。\n',encoding='utf8')
    inventory=[]
    paths=[ROOT/'AGENTS.md',ROOT/'STUDY_SPEC.md',*PACK.iterdir(),OLD/'README.md',*list((OLD/'config').glob('*')),*list((OLD/'reports').glob('*.md'))]
    needed=['metadata/new_NAM_membership_and_split.csv','metadata/shared_development_folds.csv','metadata/fold_support_and_proxy_inner_members.csv','metadata/ALT_eligibility_and_controls.csv','metadata/AllAttribute_exact_condition_dose_evidence.csv','metadata/PHH_panel_selection.csv','metadata/Hallmark_fixed_mapped_members.csv','metadata/Hallmark_platform_support.csv','derived/ALT_individual_long.csv','derived/ALT_condition_responses.csv','derived/common_ALT_PHH_query_table.csv','derived/PHH_fixed_geneset_responses.csv','derived/structure_and_PHH_condition.npz']
    paths += [OLD/p for p in needed]+list((OLD/'pilot').glob('*.csv'))+list((OLD/'pilot/contexts').glob('*'))+list((OLD/'src').glob('*.py'))+list((OLD/'logs').glob('*'))
    for p in dict.fromkeys(paths):
        if p.is_file():
            if p.suffix in ('.md','.txt','.json','.py') and '/logs/' not in p.as_posix():p.read_text(encoding='utf-8-sig')
            s=p.stat();inventory.append(dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns))
    csv(inventory,E/'metadata/readonly_input_snapshot.csv')
    members=pd.read_csv(OLD/needed[0]);valid=members.technical_final_eligibility.eq(True)
    dev=members.loc[valid & members.split.eq('NAM_dev')].copy()
    excluded=members.loc[~members.parent_id.isin(dev.parent_id),['parent_id','split','technical_final_eligibility']].copy()
    csv(dev,E/'metadata/dev_members.csv');csv(excluded,E/'metadata/excluded_ids_only.csv')
    ids=set(dev.parent_id)
    t=pd.read_csv(OLD/'derived/common_ALT_PHH_query_table.csv');t=t.loc[t.parent_id.isin(ids)].copy();guard(t.parent_id,ids)
    z=pd.read_csv(OLD/'derived/PHH_fixed_geneset_responses.csv');z=z.loc[z.parent_id.isin(ids)].copy()
    a=np.load(OLD/'derived/structure_and_PHH_condition.npz',allow_pickle=False);ix=np.isin(a['parent_id'],list(ids))
    np.savez_compressed(E/'derived/dev_structure.npz',parent_id=a['parent_id'][ix],X=a['X'][ix],columns=a['columns'])
    csv(t,E/'derived/dev_queries.csv');csv(z,E/'derived/dev_PHH50.csv')
    for p in ['metadata/shared_development_folds.csv','derived/ALT_condition_responses.csv','metadata/ALT_eligibility_and_controls.csv']:
        d=pd.read_csv(OLD/p);d=d.loc[d.parent_id.isin(ids)] if 'parent_id' in d else d.loc[d.condition_id.isin(t.condition_id)];csv(d,E/'metadata'/Path(p).name)
    animal_ids=set(';'.join(t.valid_treatment_sample_ids.tolist()+t.valid_control_sample_ids.tolist()).split(';'))
    alt=pd.read_csv(OLD/'derived/ALT_individual_long.csv');alt=alt.loc[alt.sample_id.isin(animal_ids)].copy();csv(alt,E/'derived/dev_ALT_individuals.csv')
    assert set(t.parent_id)==ids==set(z.parent_id)
    inner=pd.read_csv(OLD/'metadata/fold_support_and_proxy_inner_members.csv');csv(inner,E/'metadata/fold_support_and_proxy_inner_members.csv')
    summary=dict(dev_parents=len(dev),dev_conditions=len(t),fixed_sets=z.shape[1]-1,holdout_ids_only=int((excluded.split.eq('NAM_holdout')&excluded.technical_final_eligibility).sum()),failed_ids=excluded.loc[~excluded.technical_final_eligibility,'parent_id'].tolist(),stage='inputs_initialized_not_models_completed')
    dump(summary,E/'metadata/actual_support.json');log('initialize',**summary);print(summary)
