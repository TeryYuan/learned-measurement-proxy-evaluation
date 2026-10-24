"""Finite artificial-data checks only. Never reads participant tables."""
import tempfile
import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from runtime import ROOT, OUT, read_json, write_json, now
from modeling import FeatureMap, ProxyModel, OutcomeModel, group_folds, select_proxy, replace_z, pick_candidate
from check_statistics import main as check_statistics

def main():
    cfg=read_json(ROOT/'protocol/analysis_config.json')
    rng=np.random.default_rng(913)
    n=180
    frame=pd.DataFrame({x:rng.normal(size=n) for x in cfg['features']['numeric_X']})
    for x in cfg['features']['categorical_X']:
        frame[x]=rng.choice(['a','b','c'],n)
    frame['psu_id']=[f'ARTIFICIAL_PSU_{i%12:02d}' for i in range(n)]
    frame['person_id']=[f'ARTIFICIAL_{i}' for i in range(n)]
    frame['WTMEC2YR']=rng.uniform(1,3,n)
    frame['y_5y']=(rng.uniform(size=n)<.3).astype(int)
    z=np.column_stack([13+.3*frame.iloc[:,0]+rng.normal(size=n),14+.5*frame.iloc[:,1]+rng.normal(size=n)])
    frame=replace_z(frame,z,cfg)
    frame.loc[0,cfg['features']['numeric_X'][0]]=np.nan
    frame.loc[1,cfg['features']['categorical_X'][0]]=None
    checks=[]
    with threadpool_limits(limits=8):
        fmap=FeatureMap(cfg).fit(frame)
        a=fmap.transform(frame)
        altered=frame.copy(); altered['y_5y']=1-altered.y_5y; altered['WTMEC2YR']=1e9
        assert np.array_equal(a,fmap.transform(altered))
        assert a[:,7:14].shape==(n,7) and a[0,7]==1
        novel=frame.iloc[:5].copy(); novel[cfg['features']['categorical_X'][0]]='UNSEEN'
        assert np.isfinite(fmap.transform(novel)).all()
        checks+=['whitelist ignores outcome and weights','seven fixed missing indicators','unknown category handled']
        for cand in cfg['proxy']['candidates']:
            model=ProxyModel(cfg,cand).fit(frame,z,frame.WTMEC2YR.to_numpy())
            pred=model.predict(frame)
            assert pred.shape==z.shape and np.isfinite(pred).all()
            assert np.allclose(model.mean_,np.average(z,axis=0,weights=frame.WTMEC2YR))
            assert np.allclose(model.scale_**2,np.average((z-model.mean_)**2,axis=0,weights=frame.WTMEC2YR))
            checks.append('proxy original units and weighted target scaling: '+cand['id'])
        fold,mapping=group_folds(frame,cfg['random_state'])
        assert np.array_equal(fold,group_folds(frame,cfg['random_state'])[0])
        outer_train=frame.loc[fold!=0]
        chosen,record=select_proxy(outer_train,cfg,1,'ARTIFICIAL_NESTED_CHECK')
        held=set(frame.loc[fold==0,'psu_id'])
        for row in record['fold_records']:
            assert not held & (set(row['fit_psus'])|set(row['heldout_psus']))
            assert not set(row['fit_psus'])&set(row['heldout_psus'])
        for row in record['candidates']:
            parts=[r for r in record['fold_records'] if r['candidate']==row['candidate']]
            assert np.isclose(row['score'],np.average([r['score'] for r in parts],weights=[r['validation_weight_total'] for r in parts]))
        checks+=['deterministic PSU folds','inner selection excludes outer heldout PSU','validation folds aggregate by heldout weight']
        for idx in [0,4]:
            cand=cfg['downstream']['candidates'][idx]
            m=OutcomeModel(cfg,cand,True).fit(frame,frame.y_5y.to_numpy(),frame.WTMEC2YR.to_numpy())
            p=m.predict_proba(frame)
            assert np.array_equal(p,m.predict_proba(replace_z(frame,z,cfg)))
            for queried in [replace_z(frame,np.broadcast_to(z.mean(axis=0),z.shape),cfg),novel.assign(**{q:frame[q].iloc[:5].values for q in cfg['features']['measured_Z']})]:
                q=m.predict_proba(queried)
                assert ((q>=0)&(q<=1)).all()
            with tempfile.TemporaryDirectory(prefix='q3_artificial_') as d:
                path=d+'/model.joblib'; joblib.dump(m,path)
                assert np.array_equal(p,joblib.load(path).predict_proba(frame))
            checks.append('downstream serialization and identical measured-query roundtrip: '+cand['family'])
        assert pick_candidate([1.,1.-5e-13])==0 and pick_candidate([1.,.9])==1
        checks.append('candidate order breaks numerical ties')
    check_statistics()
    checks.append('provided hand-calculated survey and Fieller boundary checks')
    record=dict(utc=now(),status='passed',data='artificial fixtures only; no study rows read',checks=checks,
                study_fits=0,scientific_configuration_unchanged=True)
    write_json(OUT/'checks/artificial_modeling.json',record)
    print(record)

if __name__=='__main__':
    main()
