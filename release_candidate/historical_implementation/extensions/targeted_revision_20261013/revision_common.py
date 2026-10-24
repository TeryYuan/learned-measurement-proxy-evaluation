"""Scoped revision; every new output is rooted here, legacy sources read only."""
import os
for _k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_k]='1'
os.environ['PYTHONDONTWRITEBYTECODE']='1'
import sys,json,time,itertools,traceback,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
import pandas as pd
import joblib
E=Path(__file__).resolve().parent
ROOT=E.parents[1]
B0=ROOT/'extensions/nam_phase4b0'
B1=ROOT/'extensions/nam_phase4b1'
sys.path.insert(0,str(B1/'src'))
from banded import BlockTransform,DualContext,BandedModel,select
from b1_common import qweights,metrics,median
def now():return datetime.now(timezone.utc).isoformat()
def dump(o,p):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+'.tmp')
    tmp.write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x)),encoding='utf8');tmp.replace(p)
def csv(d,p):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);pd.DataFrame(d).to_csv(p,index=False,encoding='utf-8-sig')
def read(p,**kw):return pd.read_csv(p,float_precision='round_trip',**kw)
def event(kind,**kw):
    p=E/'logs/events.jsonl';p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('a',encoding='utf8') as f:f.write(json.dumps(dict(utc=now(),event=kind,**kw),default=str,ensure_ascii=False)+'\n')
def seed(*parts):return int(np.random.SeedSequence(list(parts)).generate_state(1)[0])
def split(ids,k,*parts):
    a=np.array(sorted(ids),dtype=str);np.random.Generator(np.random.PCG64(np.random.SeedSequence(list(parts)))).shuffle(a)
    return {p:f for f,group in enumerate(np.array_split(a,k)) for p in group}
def alt_data():
    m=read(B0/'metadata/new_NAM_membership_and_split.csv');m=m.loc[m.technical_final_eligibility.eq(True)].copy();ids=set(m.parent_id)
    t=read(B0/'derived/common_ALT_PHH_query_table.csv');t=t.loc[t.parent_id.isin(ids)].copy()
    z=read(B0/'derived/PHH_fixed_geneset_responses.csv').set_index('parent_id')
    a=np.load(B0/'derived/structure_and_PHH_condition.npz');x=pd.DataFrame(a['X'],index=a['parent_id'],columns=a['columns']).loc[sorted(ids)]
    assert len(m)==121 and len(t)==370 and t.condition_id.is_unique
    assert ids==set(t.parent_id)==set(z.index)==set(x.index) and z.shape[1]==50 and x.shape[1]==2057
    assert np.isfinite(z).all().all() and np.isfinite(x).all().all() and np.isfinite(t.y).all()
    return m,t,z,x
def raw(t,x,z):return t[['log10_PHH_uM','log10_rat_mg_kg']].to_numpy(),x.loc[t.parent_id].to_numpy()[:,:-1],z.loc[t.parent_id].to_numpy(),t.route.to_numpy()
def load_dili(role):
    if role=='N':
        ix=read(ROOT/'data/phase2a_N/N_membership.csv');a=np.load(ROOT/'data/phase2a_N/N_basic_features.npz')
    else:
        ix=read(ROOT/f'data/phase1/inputs/{role}_index.csv');a=np.load(ROOT/f'data/phase1/inputs/{role}_base.npz')
    assert np.array_equal(a['ids'],ix.member_id.astype(str))
    return ix,a['b'].copy(),a['d'].copy()
