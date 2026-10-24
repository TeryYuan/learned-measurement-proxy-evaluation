"""Run the frozen five-route benchmark, resuming only compatible stage caches."""
# Reporting must branch before imports of fitting libraries and before any cache contract.
# The strict contract below remains applicable to actual fitting/resumption only.
if __name__ == '__main__':
    import sys
    if any(flag in sys.argv[1:] for flag in ('--outputs-only', '--report-only')):
        from saved_reports import main as report_main
        report_main()
        raise SystemExit(0)

from runtime import limit_threads
limit_threads()
import copy
import time
import joblib
import numpy as np
import pandas as pd
import psutil
import importlib.metadata
from threadpoolctl import threadpool_limits, threadpool_info
from runtime import ROOT, OUT, now, read_json, write_json, event, require_remote_record, mark_first_fit
from modeling import ProxyModel, OutcomeModel, group_folds, select_proxy, timed_fit, pick_candidate, replace_z

IMPLEMENTATION_VERSION='v2.1'

def fit_callback(state,label,frame,model,seconds):
    if state=='start':
        mark_first_fit()
        event('fit_start', label=label, n=len(frame), model=type(model).__name__)
    else:
        mem=psutil.Process().memory_info()
        event('fit_finished', label=label, seconds=seconds, rss_bytes=mem.rss,
              peak_working_set_bytes=getattr(mem,'peak_wset',None),
              convergence_extension=getattr(model,'convergence_messages_',[]))
        if not (OUT/'logs/first_fit_timing.json').exists():
            write_json(OUT/'logs/first_fit_timing.json',dict(utc=now(),label=label,seconds=seconds,
                note='First necessary study fit only; no separate timing pilot. Different estimator families need different times.'))
            print(f'FIRST_REQUIRED_FIT {label}: {seconds:.3f} seconds',flush=True)

def contract(cfg):
    science=copy.deepcopy(cfg)
    science['scope'].pop('timebox',None)
    science.pop('status_at_delivery',None)
    return dict(implementation=IMPLEMENTATION_VERSION,science=science,
        packages={k:importlib.metadata.version(k) for k in ['numpy','pandas','scipy','scikit-learn','joblib','threadpoolctl']},
        inputs=read_json(OUT/'checks/preflight.json')['data_files'])

def save_model(path,model):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.tmp')
    joblib.dump(model,tmp,compress=3)
    tmp.replace(path)

def cached_proxy_unit(frame,predict_frame,cfg,pool_k,label,path):
    record_path=path.with_suffix('.json')
    if path.exists() and record_path.exists():
        saved=joblib.load(path)
        assert saved['fit_ids']==frame.person_id.tolist()
        return saved['model'].predict(predict_frame),read_json(record_path),saved['model']
    candidate,record=select_proxy(frame,cfg,pool_k,label,fit_callback)
    model=timed_fit(ProxyModel(cfg,candidate),frame,frame[cfg['features']['measured_Z']].to_numpy(float),
                    frame.WTMEC2YR.to_numpy(float)/pool_k,label+'/selected_refit',fit_callback)
    record.update(fit_psus=sorted(frame.psu_id.unique().tolist()),
                  predicted_psus=sorted(predict_frame.psu_id.unique().tolist()),
                  fit_n=len(frame),predict_n=len(predict_frame))
    assert not set(record['fit_psus'])&set(record['predicted_psus'])
    save_model(path,dict(model=model,fit_ids=frame.person_id.tolist()))
    write_json(record_path,record)
    return model.predict(predict_frame),record,model

def proxy_stage(frame,holdout,cfg,pool_k,stage):
    folder=OUT/'cache'/stage
    folder.mkdir(parents=True,exist_ok=True)
    fold,mapping=group_folds(frame,cfg['random_state'],cfg['proxy']['oof_folds'])
    write_json(folder/'outer_group_to_fold.json',mapping)
    oof=np.full((len(frame),2),np.nan)
    records=[]
    for k in range(3):
        train=fold!=k; val=~train
        pred,record,_=cached_proxy_unit(frame.loc[train],frame.loc[val],cfg,pool_k,
                                       stage+f'/outer{k}',folder/f'outer{k}.joblib')
        # Each inner-selection PSU is strictly within this outer training set.
        allowed=set(frame.loc[train,'psu_id'])
        held=set(frame.loc[val,'psu_id'])
        for row in record['fold_records']:
            assert set(row['fit_psus'])<=allowed and set(row['heldout_psus'])<=allowed
            assert not set(row['fit_psus'])&held and not set(row['heldout_psus'])&held
            assert not set(row['fit_psus'])&set(row['heldout_psus'])
        oof[val]=pred
        record['outer_fold']=k
        records.append(record)
        print(f'PROXY {stage} outer {k+1}/3 complete ({record["selected"]})',flush=True)
    pred,record,model=cached_proxy_unit(frame,holdout,cfg,pool_k,stage+'/full',folder/'full.joblib')
    records.append(record)
    assert np.isfinite(oof).all() and np.isfinite(pred).all()
    result=frame[['person_id','psu_id']].copy()
    result['outer_fold']=fold
    for j,z in enumerate(cfg['features']['measured_Z']):
        result['proxy_'+z]=oof[:,j]
    result.to_csv(folder/'strict_oof_proxies.csv',index=False)
    write_json(folder/'proxy_selection_records.json',records)
    return oof,pred,model

def main():
    require_remote_record()  # Precedes all study-data fitting and model-dependent work.
    cfg=read_json(ROOT/'protocol/analysis_config.json')
    expected=contract(cfg)
    cp=OUT/'cache/contract.json'
    if cp.exists():
        if read_json(cp)!=expected:
            raise RuntimeError('Cache contract differs: diagnose affected stages; do not silently reuse or refit.')
    else:
        write_json(cp,expected)
    for name,stat in expected['inputs'].items():
        p=ROOT/name
        assert p.stat().st_size==stat['size'] and p.stat().st_mtime_ns==stat['mtime_ns'], 'Input changed since preflight.'
    if (OUT/'logs/benchmark_complete.json').exists():
        print('REUSE completed benchmark; no new fits.',flush=True)
        return
    started=time.perf_counter()
    event('benchmark_started',implementation=IMPLEMENTATION_VERSION)
    c=pd.read_csv(ROOT/'data/analysis_cohort.csv',dtype={x:str for x in ['person_id','psu_id','stratum_id']})
    train=c.loc[c.role.eq('train')].reset_index(drop=True)
    val=c.loc[c.role.eq('validation')].reset_index(drop=True)
    dev=c.loc[c.period.eq('development')].reset_index(drop=True)
    test=c.loc[c.role.eq('test')].reset_index(drop=True)
    for p in ['models','tables','predictions','checks']:
        (OUT/p).mkdir(parents=True,exist_ok=True)
    with threadpool_limits(limits=4):
        write_json(OUT/'logs/threadpools.json',threadpool_info())
        oof,proxyval,_=proxy_stage(train,val,cfg,1,'train')
        rows=[]; selected={}
        for route in ['B','M','P']:
            tr=replace_z(train,oof,cfg) if route=='P' else train
            va=replace_z(val,proxyval,cfg) if route=='P' else val
            route_rows=[]
            for candidate in cfg['downstream']['candidates']:
                path=OUT/'cache/validation'/f'{route}_{candidate["id"]}.json'
                if path.exists():
                    row=read_json(path)
                else:
                    t=time.perf_counter()
                    m=timed_fit(OutcomeModel(cfg,candidate,route!='B'),tr,tr.y_5y.to_numpy(int),
                                tr.WTMEC2YR.to_numpy(float),f'validation/{route}/{candidate["id"]}',fit_callback)
                    p=m.predict_proba(va)
                    loss=float(np.average((p-va.y_5y.to_numpy())**2,weights=va.WTMEC2YR))
                    row=dict(route=route,candidate=candidate['id'],weighted_Brier=loss,
                             seconds=time.perf_counter()-t,training_n=len(tr),validation_n=len(va),
                             convergence_extension=m.convergence_messages_)
                    write_json(path,row)
                route_rows.append(row)
                print(f'VALIDATION {route} {candidate["id"]}: {row["weighted_Brier"]:.8f}',flush=True)
            best=pick_candidate([x['weighted_Brier'] for x in route_rows],cfg['downstream']['tie_absolute_tolerance'])
            selected[route]=cfg['downstream']['candidates'][best]
            for i,row in enumerate(route_rows):
                row['selected']=i==best
            rows.extend(route_rows)
        pd.DataFrame(rows).to_csv(OUT/'tables/validation_candidates.csv',index=False)
        write_json(OUT/'tables/validation_candidates.json',rows)
        write_json(OUT/'tables/selected_configurations.json',selected)
        event('downstream_selection_complete',selected=selected,final_test_performance_seen=False)
        oofdev,proxytest,g=proxy_stage(dev,test,cfg,2,'development')
        save_model(OUT/'models/g_final.joblib',g)
        models={}
        for route in ['B','M','P']:
            path=OUT/'models'/f'{route}_final.joblib'
            if path.exists():
                models[route]=joblib.load(path)
            else:
                tr=replace_z(dev,oofdev,cfg) if route=='P' else dev
                m=timed_fit(OutcomeModel(cfg,selected[route],route!='B'),tr,tr.y_5y.to_numpy(int),
                            tr.WTMEC2YR.to_numpy(float)/2,f'development/{route}/selected_refit',fit_callback)
                save_model(path,m)
                models[route]=joblib.load(path)
        mean=np.average(dev[cfg['features']['measured_Z']],axis=0,weights=dev.WTMEC2YR/2)
        inputs=replace_z(test,proxytest,cfg)
        mean_inputs=replace_z(test,np.broadcast_to(mean,(len(test),2)),cfg)
        p={
            'B':models['B'].predict_proba(test),
            'M':models['M'].predict_proba(test),
            'S':models['M'].predict_proba(inputs),
            'MEAN':models['M'].predict_proba(mean_inputs),
            'P':models['P'].predict_proba(inputs),
        }
        original_again=models['M'].predict_proba(replace_z(test,test[cfg['features']['measured_Z']],cfg))
        assert np.array_equal(p['M'],original_again)
        assert all(np.isfinite(v).all() and ((v>=0)&(v<=1)).all() for v in p.values())
        predictions=test[['person_id']].copy()
        for route,v in p.items():
            predictions['p_'+route]=v
        assert len(predictions)==7477 and predictions.person_id.is_unique
        predictions.to_csv(OUT/'predictions/test_predictions.csv',index=False)
        truth=test[['person_id']+cfg['features']['measured_Z']].copy()
        for j,z in enumerate(cfg['features']['measured_Z']):
            truth['proxy_'+z]=proxytest[:,j]
            truth['development_mean_'+z]=mean[j]
        truth.to_csv(OUT/'predictions/test_measurement_proxies.csv',index=False)
        usage={}
        for a,b in [('M','S'),('M','MEAN'),('S','MEAN')]:
            diff=p[a]-p[b]
            usage[a+'_vs_'+b]=dict(max_abs=float(np.max(abs(diff))),RMS=float(np.sqrt(np.mean(diff**2))),
                                  fraction_abs_over_1e_minus_12=float(np.mean(abs(diff)>1e-12)))
        write_json(OUT/'checks/input_use.json',dict(status='passed',observed_Z_roundtrip_max_abs=0.,
            queries=usage,scope='observed test queries; not a proof of global function independence',
            shared_M_model=True,S_fits=0,MEAN_fits=0,development_Z_mean=mean.tolist(),
            proxy_ranges={z:dict(min=float(proxytest[:,j].min()),max=float(proxytest[:,j].max()))
                          for j,z in enumerate(cfg['features']['measured_Z'])}))
    result=dict(utc=now(),seconds_this_invocation=time.perf_counter()-started,n_test=len(test),
                status='five_routes_complete',peak_working_set_bytes=getattr(psutil.Process().memory_info(),'peak_wset',None))
    write_json(OUT/'logs/benchmark_complete.json',result)
    event('benchmark_completed',**result)
    print(result,flush=True)

if __name__=='__main__':
    main()
