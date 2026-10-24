"""new_random=yes grouped splits/permutations; model_fit=yes RF/LR/calibration; model_query=yes used O+N; used_evaluation_data=yes post-hoc DILI."""
from common import *
import warnings,itertools,concurrent.futures,joblib,psutil
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import average_precision_score,roc_auc_score
from rdkit import Chem,RDLogger
RDLogger.DisableLog('rdApp.warning')
E=P/'pooled_cv';TASKS=['BSEPi','MTX_MP','ARE'];MODELS=['M0','M1','Mperm','Mchem','L2LR']
ROLE='post-hoc internal grouped reevaluation of previously used O+N; not external validation'

def directory_bytes(root):
    """Exact logical file bytes with cached directory-entry attributes; no model/data transformation."""
    total=0;stack=[str(root)]
    while stack:
        with os.scandir(stack.pop()) as entries:
            for entry in entries:
                if entry.is_dir(follow_symlinks=False):stack.append(entry.path)
                elif entry.is_file():total+=entry.stat().st_size
    return total

def metric(y,p,c=None):
    return dict(n=len(y),positive=int(np.sum(y)),positive_fraction=float(np.mean(y)),AP=float(average_precision_score(y,p)),ROC_AUC=float(roc_auc_score(y,p)),Brier_raw=float(np.mean((y-p)**2)),Brier_calibrated=float(np.mean((y-c)**2)) if c is not None and np.isfinite(c).all() else np.nan)

def prepare():
    E.mkdir(exist_ok=True);(E/'cache').mkdir(exist_ok=True);(E/'models').mkdir(exist_ok=True)
    members=readcsv(ROOT/'data/phase1/tables/main_membership.csv').set_index('member_id')
    o=readcsv(ROOT/'data/phase1/inputs/O_index.csv');n=readcsv(ROOT/'data/phase2a_N/N_membership.csv')
    ids=o.member_id.tolist()+n.member_id.tolist();assert len(ids)==len(set(ids));m=members.loc[ids].copy().reset_index();assert m.analysis.eq('main_2d').all()
    m['original_role']=['O']*len(o)+['N']*len(n);m['mixed_y']=np.r_[o.y_v1.to_numpy(int),n.y_v2.to_numpy(int)]
    m['full_v2_eligible']=m.C.astype(bool)|m.original_role.eq('N');m['label_source']=np.where(m.original_role.eq('O'),'V1','V2');m['evidence_role']=ROLE
    common=m.full_v2_eligible.to_numpy();assert m.loc[common,'y_v2'].notna().all()
    b0=np.load(ROOT/'data/phase1/inputs/O_base.npz');bn=np.load(ROOT/'data/phase2a_N/N_basic_features.npz');assert np.array_equal(b0['ids'],o.member_id) and np.array_equal(bn['ids'],n.member_id)
    b=np.vstack([b0['b'],bn['b']]);d=np.vstack([b0['d'],bn['d']]);X=np.column_stack([b,d])
    p0=readcsv(ROOT/'data/phase1/features/O_proxy_features.csv').set_index('member_id');pn=readcsv(ROOT/'data/phase2a_N/N_proxy_features.csv').set_index('member_id');proxy=pd.concat([p0,pn]).loc[ids]
    Xs={'M0':X,'L2LR':X}
    for name,prefix in [('M1','z'),('Mperm','z_perm'),('Mchem','z_chem')]:Xs[name]=np.column_stack([X,proxy[[f'{prefix}_{k}' for k in TASKS]].to_numpy()])
    overlaps=[];A={};pkeys=set(m.identity_group);pinchi={Chem.MolToInchi(Chem.MolFromSmiles(s)) for s in m.model_input_smiles}
    for k in TASKS:
        ix=readcsv(ROOT/f'data/phase1/inputs/A_{k}_index.csv');z=np.load(ROOT/f'data/phase1/inputs/A_{k}_base.npz');assert np.array_equal(ix.member_id,z['ids'])
        ii={Chem.MolToInchi(Chem.MolFromSmiles(s)) for s in ix.model_input_smiles};ko=set(ix.identity_group)&pkeys;io=ii&pinchi
        overlaps.append(dict(task=k,A_identities=len(ix),positive=int(ix.y.sum()),key_intersection=len(ko),standard_inchi_intersection=len(io),all_A_labels_observed=bool(ix.y.isin([0,1]).all()),status='verified' if not ko and not io else 'identity_conflict'))
        assert not ko and not io,'A/pooled identity conflict: do not silently change membership'
        ax=z['b'].copy();ap=E/'cache'/f'A_{k}_source.npz'
        if ap.exists():
            previous=np.load(ap);assert np.array_equal(previous['ids'],ix.member_id) and np.array_equal(previous['X'],ax) and np.array_equal(previous['y'],ix.y)
        else:np.savez_compressed(ap,ids=ix.member_id.to_numpy(str),X=ax,y=ix.y.to_numpy(int))
        A[k]=(ix,ax)
    savecsv(overlaps,E/'A_pooled_identity_check.csv',ROLE);savecsv(m,E/'label_source_by_identity.csv',ROLE)
    cohorts={name:dict(rows=np.flatnonzero(mask),y=m.loc[mask,'y_v2' if name=='P_v2_common' else 'mixed_y'].to_numpy(int)) for name,mask in [('P_mixed_full',np.ones(len(m),bool)),('P_mixed_common',common),('P_v2_common',common)]}
    savecsv([dict(cohort=name,n=len(q['rows']),positive=int(q['y'].sum()),negative=int(len(q['y'])-q['y'].sum()),evidence_role=ROLE) for name,q in cohorts.items()],E/'cohort_support.csv',ROLE)
    for name,q in cohorts.items():savecsv(m.iloc[q['rows']].assign(y=q['y']),E/f'{name}_members.csv',ROLE)
    current=dict(ids=np.array(ids),b=b,structural_X=X,**{name:v for name,v in Xs.items() if name not in ['M0','L2LR']});cache=E/'cache/pooled_inputs.npz'
    if cache.exists():
        previous=np.load(cache);assert all(np.array_equal(previous[k],v) for k,v in current.items()),'Input matrices changed; no silent cache reuse'
        event('P3','input_matrix_equivalence_checked',arrays=list(current),all_exact=True)
    else:np.savez_compressed(cache,**current)
    return m,Xs,b,A,cohorts

def split(y,groups,k,seed):
    out=np.full(len(y),-1,int)
    for f,(tr,va) in enumerate(StratifiedGroupKFold(n_splits=k,shuffle=True,random_state=seed).split(np.zeros(len(y)),y,groups)):out[va]=f
    def valid(a):return all(len(np.unique(y[a==f]))==2 and len(np.unique(y[a!=f]))==2 for f in range(k))
    if valid(out):return out,[]
    # One deterministic allocation-repair pass, based solely on group/category support.
    original=out.copy();unique=sorted(set(groups));gix={g:np.flatnonzero(groups==g) for g in unique}
    def objective(a):
        missing=sum(2-len(np.unique(y[a==f])) for f in range(k));sizes=np.array([(a==f).sum() for f in range(k)]);pos=np.array([y[a==f].sum() for f in range(k)]);return missing,float(np.sum((sizes-len(y)/k)**2)/len(y)**2+np.sum((pos-y.mean()*sizes)**2)/len(y)**2)
    best=objective(out);choice=None
    for g,h in itertools.combinations(unique,2):
        fg=out[gix[g][0]];fh=out[gix[h][0]]
        if fg==fh:continue
        temp=out.copy();temp[gix[g]]=fh;temp[gix[h]]=fg;obj=objective(temp)
        if obj<best:best=obj;choice=(g,h,temp)
    changes=[]
    if choice is not None:
        g,h,out=choice;changes=[dict(group_a=g,group_b=h,original_fold_a=int(original[gix[g][0]]),original_fold_b=int(original[gix[h][0]]),reason='one deterministic support-only swap')]
    if not valid(out):raise RuntimeError('Grouped folds still lack a class after one allowed support-only repair')
    return out,changes

def folds(m,cohorts):
    saved=E/'fold_contexts.json'
    if saved.exists():return json.loads(saved.read_text('utf8'))
    groups=m.fold_group.to_numpy(str);contexts=[];outerrows=[];innerrows=[];support=[];repairs=[];cfg=spec()['P3']
    for r in range(10):
        seed=cfg['split_seed_base']+r;of,fix=split(m.mixed_y.to_numpy(int),groups,5,seed);repairs.extend(dict(repeat=r,stage='outer',**x) for x in fix)
        outerrows.extend(dict(repeat=r,seed=seed,member_id=m.member_id[i],fold_group=groups[i],outer_fold=int(f)) for i,f in enumerate(of))
        for f in range(5):
            maps={}
            for name,c in cohorts.items():
                rows=c['rows'];tr=rows[of[rows]!=f];va=rows[of[rows]==f];seedinner=cfg['split_seed_base']+1000+100*r+f
                if name=='P_v2_common':inn=maps['P_mixed_common'].copy()
                else:
                    inn,fix=split(m.mixed_y.to_numpy(int)[tr],groups[tr],3,seedinner);repairs.extend(dict(repeat=r,outer=f,cohort=name,stage='inner',**x) for x in fix);maps[name]=inn
                y=m.y_v2.to_numpy()[tr].astype(int) if name=='P_v2_common' else m.mixed_y.to_numpy(int)[tr]
                yq=m.y_v2.to_numpy()[va].astype(int) if name=='P_v2_common' else m.mixed_y.to_numpy(int)[va]
                assert len(np.unique(y))==2 and len(np.unique(yq))==2 and not set(groups[tr])&set(groups[va])
                for j in range(3):
                    a=tr[inn!=j];v=tr[inn==j];yy=y[inn==j]
                    assert len(np.unique(yy))==2 and len(np.unique(y[inn!=j]))==2 and not set(groups[a])&set(groups[v])
                    support.append(dict(cohort=name,repeat=r,outer_fold=f,inner_fold=j,train_n=len(a),validation_n=len(v),positive=int(yy.sum()),negative=len(v)-int(yy.sum()),training_groups=len(set(groups[a])),validation_groups=len(set(groups[v])),identity_overlap=0,group_overlap=0))
                innerrows.extend(dict(cohort=name,repeat=r,outer_fold=f,member_id=m.member_id[i],inner_fold=int(j),seed=seedinner) for i,j in zip(tr,inn))
                contexts.append(dict(cohort=name,repeat=r,outer_fold=f,tr=tr.tolist(),va=va.tolist(),inner=inn.tolist(),label='y_v2' if name=='P_v2_common' else 'mixed_y'))
    savecsv(outerrows,E/'outer_folds.csv','support-fixed post-hoc grouped folds before predictions');savecsv(innerrows,E/'inner_folds.csv',ROLE);savecsv(support,E/'fold_support.csv',ROLE);savecsv(pd.DataFrame(repairs,columns=['repeat','stage','outer','cohort','group_a','group_b','original_fold_a','original_fold_b','reason']),E/'fold_repairs.csv',ROLE)
    d=pd.DataFrame(outerrows);unique=len({tuple(g.sort_values('member_id').outer_fold) for _,g in d.groupby('repeat')});dump(dict(requested_repeats=10,distinct_labelled_assignments=unique,seed_rule='61031000+repeat; inner 61031000+1000+100*repeat+outer',class_support_verified_before_predictions=True),E/'fold_summary.json');dump(contexts,saved);return contexts

def fit_inner(job):
    X,y,tr,va,params,kind=job;start=time.time()
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        model=RandomForestClassifier(**params) if kind=='RF' else make_pipeline(StandardScaler(),LogisticRegression(C=params['C'],solver='lbfgs',class_weight=None,max_iter=5000,tol=1e-6))
        model.fit(X[tr],y[tr]);pred=model.predict_proba(X[va])[:,1]
    return pred,time.time()-start,[str(w.message) for w in ws],not any(issubclass(w.category,ConvergenceWarning) for w in ws)

def context_fit(m,X,ctx,name,key,executor,save_model=False):
    folder=E/'cache'/key;folder.mkdir(parents=True,exist_ok=True)
    done=folder/'completed.json'
    if done.exists():return readcsv(folder/'outer_predictions.csv'),json.loads(done.read_text('utf8'))
    tr=np.array(ctx['tr']);va=np.array(ctx['va']);inn=np.array(ctx['inner']);y=m[ctx['label']].fillna(-99).to_numpy(int);start=time.time();kind='LR' if name=='L2LR' else 'RF'
    oldcfg=json.loads((ROOT/'config/phase1_execution.json').read_text('utf8'));grids=[dict(C=c) for c in spec()['P3']['LR_C']] if kind=='LR' else [{**oldcfg['downstream'],**p,'n_jobs':1} for p in oldcfg['grid']]
    expected=dict(train_ids=m.member_id.iloc[tr].tolist(),query_ids=m.member_id.iloc[va].tolist(),training_labels=y[tr].tolist(),inner=inn.tolist(),input_dimension=X.shape[1],kind=kind,grid=grids,random_state=4109)
    meta=folder/'context.json'
    if meta.exists():assert json.loads(meta.read_text('utf8'))==expected
    else:dump(expected,meta)
    pending=[];innerpred={i:np.full(len(tr),np.nan) for i in range(len(grids))};scores=[];fitrecords=[];allvalid={i:True for i in range(len(grids))}
    for gi,params in enumerate(grids):
        for j in range(3):
            it=tr[inn!=j];iv=tr[inn==j];cp=folder/f'g{gi}_inner{j}.npz'
            if cp.exists():
                a=np.load(cp);innerpred[gi][inn==j]=a['prediction'];allvalid[gi]&=bool(a['valid']);event('P3','inner_cache_reused',key=key,candidate=gi,inner=j)
            else:pending.append((gi,j,cp,executor.submit(fit_inner,(X,y,it,iv,params,kind))))
    for gi,j,cp,future in pending:
        pred,seconds,warnings_,valid=future.result();innerpred[gi][inn==j]=pred;allvalid[gi]&=valid;np.savez_compressed(cp,prediction=pred,valid=valid)
        rec=dict(key=key,candidate=gi,inner=j,kind=kind,training_n=int((inn!=j).sum()),query_n=int((inn==j).sum()),elapsed_seconds=seconds,valid=valid,warnings=warnings_);event('P3','fit_and_query_completed',**rec);fitrecords.append(rec)
    for gi,params in enumerate(grids):
        p=innerpred[gi];assert np.isfinite(p).all();ap=[float(average_precision_score(y[tr[inn==j]],p[inn==j])) for j in range(3)]
        scores.append(dict(candidate=gi,mean_AP=float(np.mean(ap)),fold_AP=ap,valid=allvalid[gi],params=params))
    valid=[s for s in scores if s['valid']]
    if not valid:raise RuntimeError('All candidates invalid: '+key)
    best=max(s['mean_AP'] for s in valid);tied=[s for s in valid if best-s['mean_AP']<=1e-12]
    win=min(tied,key=lambda s:s['params']['C'] if kind=='LR' else (-s['params']['min_samples_leaf'],0 if s['params']['max_features']=='sqrt' else 1));gi=win['candidate'];pi=innerpred[gi]
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always');cal=LogisticRegression(C=np.inf,l1_ratio=0.,solver='lbfgs',class_weight=None,max_iter=2000,tol=1e-9).fit(pi.reshape(-1,1),y[tr]);calok=not any(issubclass(w.category,ConvergenceWarning) for w in ws)
    event('P3','calibrator_fit_completed',key=key,training_n=len(tr),valid=calok,warnings=[str(w.message) for w in ws])
    params=win['params'];model=RandomForestClassifier(**params) if kind=='RF' else make_pipeline(StandardScaler(),LogisticRegression(C=params['C'],solver='lbfgs',class_weight=None,max_iter=5000,tol=1e-6))
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always');model.fit(X[tr],y[tr]);raw=model.predict_proba(X[va])[:,1]
    if any(issubclass(w.category,ConvergenceWarning) for w in ws):raise RuntimeError('Selected refit did not converge: '+key)
    calibrated=cal.predict_proba(raw.reshape(-1,1))[:,1] if calok else np.full(len(va),np.nan)
    event('P3','selected_fit_and_query_completed',key=key,kind=kind,train_n=len(tr),query_n=len(va),params=params,calibration_query=calok)
    out=pd.DataFrame(dict(member_id=m.member_id.iloc[va].to_numpy(),fold_group=m.fold_group.iloc[va].to_numpy(),y=y[va],raw_score=raw,calibrated_score=calibrated,repeat=ctx['repeat'],outer_fold=ctx['outer_fold'],model=name,cohort=ctx['cohort'],evidence_role=ROLE))
    savecsv(out,folder/'outer_predictions.csv',ROLE);savecsv(pd.DataFrame(dict(member_id=m.member_id.iloc[tr].to_numpy(),y=y[tr],raw_score=pi,inner_fold=inn)),folder/'selected_inner_oof.csv','current-training OOF used for selection and calibration; not independent assessment')
    dump(scores,folder/'candidate_scores.json')
    if save_model:joblib.dump(dict(model=model,calibrator=cal if calok else None,context=expected,selected_candidate=win,role=ROLE),E/'models'/f'{key}.joblib',compress=3)
    summary=dict(status='complete',key=key,kind=kind,selected=win,calibration_valid=calok,elapsed_seconds=time.time()-start,new_inner_fits=len(pending),new_selected_fits=1,new_calibration_fits=1,metrics=metric(y[va],raw,calibrated),rss_gib=psutil.Process().memory_info().rss/2**30)
    dump(summary,done);return out,summary

def perm_proxy(m,b,A,number,only_first=False):
    file=E/'cache'/f'permutation_{number:03d}_proxy.npz'
    if file.exists():
        a=np.load(file);assert np.array_equal(a['ids'],m.member_id)
        for k,task in enumerate(TASKS):
            saved=np.load(E/'cache'/f'permutation_{number:03d}_{task}.npz');ix,_=A[task]
            assert np.array_equal(saved['training_ids'],ix.member_id) and int(saved['seed'])==spec()['P3']['permutation_seed_base']+100*number+k and np.array_equal(saved['query_ids'],m.member_id)
            assert np.array_equal(np.sort(saved['permuted_y']),np.sort(ix.y)) and np.array_equal(a['proxy'][:,k],saved['prediction'])
        return a['proxy']
    values=[];maps=[];oldcfg=json.loads((ROOT/'config/phase1_execution.json').read_text('utf8'))
    for k,task in enumerate(TASKS):
        p=E/'cache'/f'permutation_{number:03d}_{task}.npz';ix,X=A[task];seed=spec()['P3']['permutation_seed_base']+100*number+k
        if p.exists():a=np.load(p);pred=a['prediction'];yp=a['permuted_y'];event('P3','permuted_auxiliary_reused',b=number,task=task)
        else:
            start=time.time();yp=np.random.Generator(np.random.PCG64(seed)).permutation(ix.y.to_numpy(int));params=oldcfg['aux_classifier'];model=RandomForestClassifier(**params).fit(X,yp);pred=model.predict_proba(b)[:,1]
            np.savez_compressed(p,prediction=pred,permuted_y=yp,training_ids=ix.member_id.to_numpy(str),query_ids=m.member_id.to_numpy(str),seed=seed)
            joblib.dump(model,E/'models'/f'permutation_{number:03d}_{task}.joblib',compress=3)
            event('P3','new_permuted_auxiliary_fit_query',b=number,task=task,seed=seed,training_n=len(ix),query_n=len(m),elapsed_seconds=time.time()-start,trees=500)
            if number==1 and k==0:dump(dict(task=task,b=number,elapsed_seconds=time.time()-start,trees=500,rss_gib=psutil.Process().memory_info().rss/2**30),E/'technical_auxiliary_timing.json')
        assert np.array_equal(np.sort(yp),np.sort(ix.y));values.append(pred);maps.extend(dict(b=number,task=task,seed=seed,member_id=id_,original_label=int(y),permuted_label=int(p_)) for id_,y,p_ in zip(ix.member_id,ix.y,yp))
        if only_first:return None
    savecsv(maps,E/f'permutation_{number:03d}_labels.csv','task-wise global observed A-label permutation; joint three-task reference')
    values=np.column_stack(values);np.savez_compressed(file,proxy=values,ids=m.member_id.to_numpy(str));return values

def run(resume=False,estimate=False):
    start=time.time();m,Xs,b,A,cohorts=prepare();contexts=folds(m,cohorts);cfg=spec()['P3'];budget=spec()['resources']['P3_wall_budget_hours']*3600
    protocol=dict(locked_utc=now(),role=ROLE,config=cfg,execution_workers=4,per_forest_workers=1,auxiliary_workers_when_pool_idle=4,sequence='main three cohorts (all five routes), then b=1..100 with complete 50-outer jobs',calibration='selected current-training inner OOF only; AP/AUC raw',threshold='not selected: P3 ranking and probability metrics only',permutation_seed='61032000+100*b+task_index; b1..100, task order BSEPi MTX_MP ARE',inner_seed='61031000+1000+100*repeat+outer',unchanged_historical_outputs=True)
    if not (E/'protocol.json').exists():dump(protocol,E/'protocol.json')
    if not (E/'technical_auxiliary_timing.json').exists():perm_proxy(m,b,A,1,only_first=True)
    statepath=E/'resource_clock.json'
    if statepath.exists():resource=json.loads(statepath.read_text('utf8'));previous_elapsed=resource.get('consumed_active_seconds',0)
    else:resource=dict(started_utc=now(),consumed_active_seconds=0);previous_elapsed=0
    def limit():
        scanstart=time.perf_counter();disk=directory_bytes(P)/2**30;scanseconds=time.perf_counter()-scanstart;rss=psutil.Process().memory_info().rss/2**30
        elapsed=previous_elapsed+time.time()-start;resource.update(consumed_active_seconds=elapsed,updated_utc=now(),disk_gib=disk,rss_gib=rss,last_disk_scan_seconds=scanseconds,scanner='os.scandir exact logical-byte traversal');dump(resource,statepath)
        reason='wall_budget' if elapsed>=budget else 'disk_soft_limit' if disk>=20 else 'memory_soft_limit' if rss>=24 else ''
        if reason:dump(dict(reason=reason,elapsed_seconds=elapsed,disk_gib=disk,rss_gib=rss),E/'resource_checkpoint.json')
        return bool(reason)
    allmain=[];complete_main=0;completed_b=[];interrupted=False
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for ctx in contexts:
            for name in MODELS:
                key=f"{ctx['cohort']}_r{ctx['repeat']:02d}_f{ctx['outer_fold']}_{name}"
                out,summary=context_fit(m,Xs[name],ctx,name,key,pool,save_model=True);allmain.append(out);complete_main+=1
                if not (E/'technical_downstream_timing.json').exists():
                    dump(dict(**summary,complete_grid='6 x 3 inner fits + selected 500-tree refit, train-only sigmoid',permutation_contexts_expected=5000,naive_context_time_projection_hours=summary['elapsed_seconds']*5000/3600,not_result_driven=True),E/'technical_downstream_timing.json')
                    if estimate:return
                if complete_main%5==0:
                    dump(dict(completed_main_contexts=complete_main,expected_main_contexts=750,completed_permutations=completed_b,updated_utc=now()),E/'progress.json');print(f'P3 main contexts {complete_main}/750',flush=True)
                if limit():interrupted=True;break
            if interrupted:break
        if complete_main==750:
            main=pd.concat(allmain,ignore_index=True);savecsv(main,E/'pooled_oof_predictions.csv',ROLE)
            for number in range(1,101):
                bp=E/f'permutation_{number:03d}_completed.json'
                if bp.exists():completed_b.append(number);continue
                if limit():interrupted=True;break
                proxy=perm_proxy(m,b,A,number);X=np.column_stack([Xs['M0'],proxy]);outs=[]
                for ctx in [c for c in contexts if c['cohort']=='P_mixed_full']:
                    key=f"perm{number:03d}_r{ctx['repeat']:02d}_f{ctx['outer_fold']}";out,summary=context_fit(m,X,ctx,f'Mperm_{number:03d}',key,pool);outs.append(out)
                    if limit():interrupted=True;break
                if interrupted:break
                dd=pd.concat(outs,ignore_index=True);savecsv(dd,E/f'permutation_{number:03d}_oof.csv','100-set post-hoc random auxiliary supervision reference; not confirmatory increment P value')
                dump(dict(status='complete',b=number,outer_contexts=50,auxiliary_models=3,rows=len(dd),updated_utc=now()),bp);completed_b.append(number)
                dump(dict(completed_main_contexts=750,expected_main_contexts=750,completed_permutations=completed_b,expected_permutations=100,updated_utc=now()),E/'progress.json');print(f'P3 permutations {len(completed_b)}/100',flush=True)
    if allmain:savecsv(pd.concat(allmain,ignore_index=True),E/'available_main_oof.csv',ROLE)
    status='models_complete_statistics_pending' if len(completed_b)==100 else 'resource_checkpoint_incomplete'
    finish('P3',status,dict(completed_main_contexts=complete_main,expected_main_contexts=750,completed_permutations=completed_b,expected_permutations=100,elapsed_seconds=previous_elapsed+time.time()-start,current_invocation_elapsed_seconds=time.time()-start,used_N_as_training=True,original_N_scores_unchanged=True),f'# P3 实际执行\n\n主分析上下文{complete_main}/750；完整联合置换{len(completed_b)}/100。所有新增结果属于已用O/N的后验内部再评价。原N固定主评价未改变。\n\n模型/校准及抽样调用见events.jsonl；未完成计算不写成完整100套。\n')
    if complete_main==750:
        from p3_statistics import run_statistics
        run_statistics(m,cohorts,completed_b)
