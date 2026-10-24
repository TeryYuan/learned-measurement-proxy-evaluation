from revision_common import *
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score,roc_auc_score
sys.path.insert(0,str(ROOT/'src'))
import dili_phase1 as old
# Reuse exact selection/calibration/threshold rules with logs scoped to this extension.
old.log_event=lambda kind,**kw:event('seed_'+kind,**kw)
def run(root):
    dest=E/f'seeds/seed{root}';dest.mkdir(exist_ok=True)
    if (dest/'completed.json').exists():return json.loads((dest/'completed.json').read_text())
    tick=time.perf_counter();cfg=json.loads((E/'protocol.json').read_text())['DILI']['original_config'];o,ob,od=load_dili('O');n,nb,nd=load_dili('N')
    ids=o.member_id.to_numpy(str);y=o.y_v1.to_numpy(int);nids=n.member_id.to_numpy(str)
    f=read(E/'seeds/fixed_O_inner_folds.csv');assert np.array_equal(f.member_id,ids);fold=f.fold.to_numpy(int)
    po=o[['member_id','identity_group']].copy();pn=n[['member_id','identity_group']].copy();fit_count=0;new_fits=0
    for j,task in enumerate(['BSEPi','MTX_MP','ARE']):
        a,ab,_=load_dili('A_'+task);par={**cfg['aux_classifier'],'random_state':seed(root,1,j),'n_jobs':4};p=dest/f'aux_{task}.joblib'
        if p.exists():
            obj=joblib.load(p);assert obj['train_ids']==a.member_id.tolist() and obj['params']==par;model=obj['model']
        else:
            start=time.perf_counter();model=RandomForestClassifier(**par).fit(ab,a.y.to_numpy(int));new_fits+=1
            joblib.dump(dict(model=model,train_ids=a.member_id.tolist(),params=par,task=task),p,compress=3)
            event('seed_RF_fit',root=root,context='aux_'+task,n=len(a),params=par,seconds=time.perf_counter()-start)
        fit_count+=1;po['z_'+task]=model.predict_proba(ob)[:,1];pn['z_'+task]=model.predict_proba(nb)[:,1];del model
    csv(po,dest/'O_proxies.csv');csv(pn,dest/'N_proxies.csv')
    xo=np.column_stack([ob,od]);xn=np.column_stack([nb,nd]);tasks=['z_BSEPi','z_MTX_MP','z_ARE'];allrows=[];settings=[];predictions=[]
    for name in ['M0','M1']:
        X=xo if name=='M0' else np.column_stack([xo,po[tasks]]);NX=xn if name=='M0' else np.column_stack([xn,pn[tasks]])
        rows=[];gridpred={};foldrows=[]
        for gi,g in enumerate(cfg['grid']):
            pp=np.full(len(o),np.nan);aps=[]
            for j in sorted(set(fold)):
                tr=np.flatnonzero(fold!=j);va=np.flatnonzero(fold==j);par={**cfg['downstream'],**g,'random_state':seed(root,2,gi,int(j)),'n_jobs':4}
                ck=dest/f'{name}_grid{gi}_inner{j}.npz';meta=ck.with_suffix('.json');expected=dict(params=par,train_ids=ids[tr].tolist(),query_ids=ids[va].tolist())
                if ck.exists() and meta.exists():
                    assert json.loads(meta.read_text())==expected;pv=np.load(ck)['prediction']
                else:
                    start=time.perf_counter();model=RandomForestClassifier(**par).fit(X[tr],y[tr]);pv=model.predict_proba(X[va])[:,1];new_fits+=1
                    np.savez_compressed(ck,prediction=pv);dump(expected,meta);del model
                    event('seed_RF_fit',root=root,context=f'{name}_grid{gi}_inner{j}',params=par,seconds=time.perf_counter()-start)
                fit_count+=1;pp[va]=pv;ap=float(average_precision_score(y[va],pv));aps.append(ap);foldrows.append(dict(model=name,candidate=gi,inner=int(j),AP=ap,random_state=par['random_state']))
            assert np.isfinite(pp).all();gridpred[gi]=pp;rows.append(dict(model=name,candidate=gi,mean_AP=float(np.mean(aps)),**g))
            print(f'seed {root} {name} grid {gi+1}/6 complete',flush=True)
        winner=old.choose_grid(rows);gi=winner['candidate'];internal=gridpred[gi]
        cal,info=old.calibrate(internal,y,ids,f'revision_seed{root}_{name}')
        threshold,tinfo=old.select_threshold(y,internal)
        par={**cfg['downstream'],**cfg['grid'][gi],'random_state':seed(root,2,99),'n_jobs':4};path=dest/f'{name}_final.joblib'
        if path.exists():
            obj=joblib.load(path);assert obj['train_ids']==ids.tolist() and obj['params']==par;rf=obj['rf']
        else:
            start=time.perf_counter();rf=RandomForestClassifier(**par).fit(X,y);new_fits+=1
            joblib.dump(dict(rf=rf,calibrator=cal,threshold=threshold,train_ids=ids.tolist(),params=par,calibration=info,feature_dimension=X.shape[1]),path,compress=3)
            event('seed_RF_fit',root=root,context=f'{name}_final',params=par,seconds=time.perf_counter()-start)
        fit_count+=1
        raw=rf.predict_proba(NX)[:,1];calp=cal.predict_proba(raw[:,None])[:,1] if cal is not None else np.full(len(n),np.nan)
        out=pd.DataFrame(dict(root_seed=root,model=name,member_id=nids,raw_score=raw,calibrated_score=calp,threshold=threshold,predicted_class=(raw>=threshold).astype(int)))
        predictions.append(out);csv(out,dest/f'{name}_N_predictions_before_scoring.csv')
        # N outcomes are read only for scoring after all model choices for this branch.
        yn=n.y_v2.to_numpy(int);met=dict(root_seed=root,model=name,AP=float(average_precision_score(yn,raw)),AUC=float(roc_auc_score(yn,raw)),Brier_raw=float(np.mean((yn-raw)**2)),Brier_calibrated=float(np.mean((yn-calp)**2)) if np.isfinite(calp).all() else None,candidate=gi,**cfg['grid'][gi])
        allrows.append(met);settings.append(dict(model=name,selected=winner,params=par,calibration=info,threshold=threshold,threshold_info=tinfo))
        csv(pd.DataFrame(dict(member_id=ids,inner=fold,y_v1=y,raw_score=internal)),dest/f'{name}_selected_O_inner_predictions.csv');csv(rows,dest/f'{name}_grid_scores.csv');csv(foldrows,dest/f'{name}_fold_scores.csv')
        if root==2026101401:
            loaded=joblib.load(path);err=float(abs(loaded['rf'].predict_proba(NX)[:,1]-raw).max());assert err<1e-12;dump(dict(model=name,reload_max_abs_error=err),dest/f'{name}_reload.json')
    assert fit_count==41
    csv(allrows,dest/'metrics.csv');csv(pd.concat(predictions),dest/'N_predictions.csv');dump(settings,dest/'selected_parameters.json')
    result=dict(root_seed=root,completed_utc=now(),elapsed_seconds=time.perf_counter()-tick,nominal_fits=fit_count,new_fits=new_fits,DeltaAP=allrows[1]['AP']-allrows[0]['AP'])
    dump(result,dest/'completed.json');event('seed_complete',**result);print(result,flush=True);return result
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--seed',type=int,required=True);a=p.parse_args();run(a.seed)
