"""Bounded input and environment checks; never fit a study model."""
import importlib.metadata as metadata
import platform
import sys
import numpy as np
import pandas as pd
from runtime import ROOT, OUT, write_json, read_json, now, event
from fetch_inputs import looks_complete
from prepare_cohort import endpoint_5y

def main():
    destination = OUT/'checks/preflight.json'
    if destination.exists():
        print('REUSE passed preflight:', destination)
        return read_json(destination)
    cfg = read_json(ROOT/'protocol/analysis_config.json')
    manifest = read_json(ROOT/'sources/source_manifest.json')
    files = manifest['data_files']
    bad = [x['file'] for x in files if not looks_complete(ROOT/x['local_path'], x)]
    assert len(files) == 40 and not bad, bad
    dtype = {x:str for x in ['person_id','psu_id','stratum_id']}
    c = pd.read_csv(ROOT/'data/analysis_cohort.csv', dtype=dtype)
    f = pd.read_csv(ROOT/'data/design_frame.csv', dtype=dtype)
    x = cfg['features']['numeric_X']+cfg['features']['categorical_X']
    z = cfg['features']['measured_Z']
    assert len(x)==16 and not set(x)&set(cfg['features']['forbidden_X'])
    assert len(c)==13236 and len(f)==39608
    assert not c.person_id.duplicated().any() and not f.person_id.duplicated().any()
    pd.testing.assert_frame_equal(c.reset_index(drop=True), f.loc[f.is_eligible].reset_index(drop=True), check_dtype=False)
    assert np.isfinite(c[z].to_numpy()).all() and (c[z]>0).all().all()
    assert c.RIDAGEYR.ge(40).all() and c.age_cap80.eq(c.RIDAGEYR.clip(upper=80)).all()
    assert c.ELIGSTAT.eq(1).all()
    assert np.array_equal(c.y_5y.to_numpy(), endpoint_5y(c.MORTSTAT,c.PERMTH_EXM).to_numpy())
    assert np.allclose(c.weight, c.WTMEC2YR/2) and np.allclose(f.weight, f.WTMEC2YR/2)
    expected = {'train':(2955,370),'validation':(2804,288),'test':(7477,648)}
    counts = {}
    for role, (n,e) in expected.items():
        sub=c.loc[c.role.eq(role)]
        assert len(sub)==n and sub.y_5y.sum()==e
        wanted = cfg['cycles'][role]
        assert set(sub.cycle)==set(wanted)
        counts[role]=dict(n=n, events=e, cycles=wanted)
    t=f.loc[f.role.eq('test')]
    psu=t[['stratum_id','psu_id']].drop_duplicates()
    assert len(psu)==63 and psu.stratum_id.nunique()==31
    assert (psu.groupby('stratum_id').size()>=2).all()
    env=dict(python=sys.version, executable=sys.executable, platform=platform.platform(),
             packages={k:metadata.version(k) for k in ['numpy','pandas','scipy','scikit-learn','matplotlib','joblib','threadpoolctl','psutil']},
             CPU_threads_max=8, outer_parallel_jobs=1, GPU=False)
    write_json(OUT/'logs/environment.json',env)
    locked='\n'.join(f'{k}=={v}' for k,v in env['packages'].items())+'\n'
    (ROOT/'requirements.lock.txt').write_text(locked,encoding='utf-8')
    result=dict(utc=now(), status='passed', raw_files=40, raw_checks='recorded size and format marker',
        raw_downloads=0, cohort_rebuilds=0, study_fits=0, analysis_n=len(c), full_frame_n=len(f),
        counts_by_role=counts, test_frame_n=len(t), test_strata=31, test_psus=63, test_df=32,
        X_whitelist=x, measured_Z=z, provided_descriptives_reused=True,
        missing_X_by_role={r:{col:int(sub[col].isna().sum()) for col in x} for r,sub in c.groupby('role')},
        data_files={name:dict(size=(ROOT/name).stat().st_size,mtime_ns=(ROOT/name).stat().st_mtime_ns)
                    for name in ['data/analysis_cohort.csv','data/design_frame.csv']})
    write_json(destination,result)
    event('preflight_passed', analysis_n=len(c), study_fits=0)
    print(dict(status='passed',raw_files=40,analysis_n=len(c),full_frame_n=len(f),counts=counts))
    return result

if __name__=='__main__':
    main()
