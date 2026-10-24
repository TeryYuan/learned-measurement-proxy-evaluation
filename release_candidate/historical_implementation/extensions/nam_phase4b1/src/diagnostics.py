from b1_common import *

def run():
    m,t,z,x=data();ids=set(m.parent_id);q=qweights(t.parent_id);y=t.y.to_numpy();diag=E/'diagnostics'
    # Response distribution with explicitly defined empirical weighted quantiles.
    out=[]
    for name,w in [('condition_equal',np.ones(len(t))/len(t)),('compound_equal',q)]:
        order=np.argsort(y,kind='stable');cum=np.cumsum(w[order]);quant={f'p{p}':float(y[order[min(np.searchsorted(cum,p/100-1e-14),len(y)-1)]]) for p in [0,5,25,50,75,95,100]}
        out.append(dict(weighting=name,mean=w@y,SD=np.sqrt(w@(y-w@y)**2),mean_abs=w@abs(y),root_mean_square=np.sqrt(w@y**2),**quant))
    csv(out,diag/'target_distribution.csv')
    ym=t.groupby('parent_id').y.mean();within=t.y-t.parent_id.map(ym);total=q@(y-q@y)**2;between=np.mean((ym-ym.mean())**2);vwithin=q@within.to_numpy()**2
    tr=set(';'.join(t.valid_treatment_sample_ids).split(';'));co=set(';'.join(t.valid_control_sample_ids).split(';'))
    summary=dict(n_parents=len(ids),n_conditions=len(t),n_unique_treated=len(tr),n_unique_controls=len(co),V_total=total,V_between=between,V_within=vwithin,identity_residual=total-between-vwithin)
    csv([summary],diag/'target_support_and_variance.csv')
    # Source animal uncertainty and the exact shared-control graph.
    alt=pd.read_csv(E/'derived/dev_ALT_individuals.csv').set_index('sample_id');assert alt.index.is_unique
    rows=[];edges=[];mapping=[];sets={}
    for r in t.itertuples():
        ti=r.valid_treatment_sample_ids.split(';');ci=r.valid_control_sample_ids.split(';');sets[r.condition_id]=set(ci)
        a=alt.loc[ti,'ALT_log2'].to_numpy();b=alt.loc[ci,'ALT_log2'].to_numpy();nt=len(a);nc=len(b)
        assert np.isfinite(a).all() and np.isfinite(b).all() and min(nt,nc)>=2
        vt=np.var(a,ddof=1);vc=np.var(b,ddof=1);vhat=vt/nt+vc/nc
        rows.append(dict(condition_id=r.condition_id,parent_id=r.parent_id,source_name=r.source_name,y=r.y,n_treated=nt,n_control=nc,s2_treated=vt,s2_control=vc,control_covariance_component=vc/nc,vhat=vhat,estimated_SE=np.sqrt(vhat),mean_log_reconstruction_error=a.mean()-b.mean()-r.y,treated_ids=';'.join(ti),control_ids=';'.join(ci)))
        for role,samples in [('treated',ti),('control',ci)]:
            for sid in samples:mapping.append(dict(parent_id=r.parent_id,condition_id=r.condition_id,role=role,sample_id=sid,ALT_source_file=alt.loc[sid,'ALT_source_file'],ALT_source_row=alt.loc[sid,'ALT_source_row']))
    df=pd.DataFrame(rows);csv(df,diag/'ALT_sampling_variance_and_shared_controls.csv');csv(mapping,diag/'ALT_condition_animal_source_mapping.csv')
    dfi=df.set_index('condition_id')
    for a,b in itertools.combinations(sets,2):
        shared=sets[a]&sets[b]
        if shared:edges.append(dict(condition_a=a,condition_b=b,parent_a=dfi.loc[a,'parent_id'],parent_b=dfi.loc[b,'parent_id'],shared_control_ids=';'.join(sorted(shared)),n_shared=len(shared),complete_same_control=sets[a]==sets[b],covariance_component=dfi.loc[a,'control_covariance_component'] if sets[a]==sets[b] else np.nan,partial_overlap_reason='' if sets[a]==sets[b] else 'actual IDs retained; complete-control covariance not assigned'))
    csv(edges,diag/'shared_control_edges.csv')
    reuse=pd.DataFrame(mapping);reuse=reuse[reuse.role.eq('control')].groupby('sample_id',as_index=False).agg(n_conditions=('condition_id','nunique'),condition_ids=('condition_id',lambda s:';'.join(s)),n_parents=('parent_id','nunique'));csv(reuse,diag/'control_reuse_patterns.csv')
    contribution=t[['parent_id','condition_id','source_name','y']].copy();contribution['q']=q;contribution['contribution_to_Vtotal']=q*(y-q@y)**2;contribution['contribution_to_zero_MSE']=q*y**2
    csv(contribution.sort_values('contribution_to_Vtotal',ascending=False),diag/'all_target_contributions_descriptive_not_exclusions.csv')
    # All fixed source repetitions. No CEL, normalization, mapping or set selection repeated.
    members=pd.read_csv(OLD/'metadata/Hallmark_fixed_mapped_members.csv',dtype=str);gene_lists=members.groupby('geneset').Entrez_gene.apply(list).to_dict()
    halves=[];cancel=[];source=[];readfiles=[]
    for pid in sorted(ids):
        p=OLD/f'derived/expression_by_parent/{pid}_gene_log2.csv.gz';g=pd.read_csv(p,dtype={'Entrez_gene':str}).set_index('Entrez_gene')
        sp=OLD/f'metadata/archive_directories/{pid}_acquired.csv';s=pd.read_csv(sp,dtype=str)
        gp=OLD/f'derived/expression_by_parent/{pid}_gene_response.csv';saved=pd.read_csv(gp,dtype={'Entrez_gene':str}).set_index('Entrez_gene')
        halves_role={}
        for role in ['treated','control']:
            rs=s.loc[s.role.eq(role)].sort_values('sample_id',kind='stable');assert len(rs)>=2
            split=np.array_split(np.arange(len(rs)),2);halves_role[role]=[rs.iloc[a].BARCODE.tolist() for a in split]
            for j,ix in enumerate(split):
                for r in rs.iloc[ix].itertuples():source.append(dict(parent_id=pid,role=role,half=j+1,sample_id=r.sample_id,BARCODE=r.BARCODE,source_file=r.source_file,source_row=r.source_row,gene_matrix=str(p.relative_to(ROOT))))
        r1=g[halves_role['treated'][0]].mean(axis=1)-g[halves_role['control'][0]].mean(axis=1)
        r2=g[halves_role['treated'][1]].mean(axis=1)-g[halves_role['control'][1]].mean(axis=1)
        full=g[sum(halves_role['treated'],[])].mean(axis=1)-g[sum(halves_role['control'],[])].mean(axis=1)
        assert np.allclose(full,saved.loc[g.index,'response'],atol=1e-12,rtol=0)
        for name in z.columns:
            genes=gene_lists[name];a=float(r1.loc[genes].mean());b=float(r2.loc[genes].mean());r=full.loc[genes].to_numpy();ab=float(abs(r.mean()));ma=float(abs(r).mean())
            halves.append(dict(parent_id=pid,geneset=name,half1=a,half2=b,difference=a-b,n_mapped_genes=len(genes)))
            cancel.append(dict(parent_id=pid,geneset=name,abs_mean_response=ab,mean_absolute_gene_response=ma,cancellation_ratio=ab/ma if ma>0 else np.nan,full_response=r.mean(),reconstruction_error=r.mean()-z.loc[pid,name]))
        for fp in [p,sp,gp]:
            st=fp.stat();readfiles.append(dict(path=str(fp),bytes=st.st_size,mtime_ns=st.st_mtime_ns,scope='NAM_dev_only_actual_read'))
    csv(readfiles,E/'metadata/PHH_actual_read_inventory.csv');hf=pd.DataFrame(halves);ca=pd.DataFrame(cancel)
    csv(hf,diag/'PHH_repeat_half_responses.csv');csv(source,diag/'PHH_repeat_half_source_mapping.csv');csv(ca,diag/'Hallmark_cancellation_per_parent.csv')
    rep=[]
    for name,g in hf.groupby('geneset'):
        rep.append(dict(geneset=name,n_parents=len(g),half_response_correlation=np.corrcoef(g.half1,g.half2)[0,1],mean_half_difference=g.difference.mean(),mean_absolute_half_difference=abs(g.difference).mean(),RMS_half_difference=np.sqrt(np.mean(g.difference**2)),half1_SD=g.half1.std(ddof=0),half2_SD=g.half2.std(ddof=0),limitation='source repeats; one treated and one control per half; donors not established'))
    csv(rep,diag/'PHH_repeat_response_stability.csv');csv(ca.groupby('geneset',as_index=False).agg(n_parents=('parent_id','size'),mean_abs_mean_response=('abs_mean_response','mean'),mean_absolute_gene_response=('mean_absolute_gene_response','mean'),mean_cancellation_ratio=('cancellation_ratio','mean'),median_cancellation_ratio=('cancellation_ratio','median')),diag/'Hallmark_cancellation_summary.csv')
    # Read all saved OOF and actual old models, checking input context and baseline reconciliation.
    old=pd.read_csv(OLD/'pilot/NAM_dev_OOF_predictions.csv');guard(old.parent_id,ids)
    pr=pd.read_csv(OLD/'pilot/proxy_unique_parent_OOF.csv');guard(pr.parent_id,ids)
    proxy_metrics=pd.read_csv(OLD/'pilot/proxy_per_geneset_metrics.csv');old_errors=pd.read_csv(OLD/'pilot/per_compound_errors.csv');guard(old_errors.parent_id,ids)
    manifest=pd.read_csv(OLD/'pilot/model_manifest.csv');csv(manifest,E/'metadata/legacy_model_manifest_read.csv')
    baselines=[];checks=[];newpred=[];modelrows=[];meanrefs=[]
    for f in range(3):
        ctx=json.loads((OLD/f'pilot/contexts/outer{f}_input_context.json').read_text(encoding='utf8'))
        S=ctx['training_parents'];Q=ctx['validation_parents'];guard(S+Q,ids);assert not set(S)&set(Q)
        assert ctx['geneset_order']==z.columns.tolist() and ctx['proxy_input_columns']==x.columns.tolist()
        ts=t.loc[t.parent_id.isin(S)];tq=t.loc[t.parent_id.isin(Q)];w=qweights(ts.parent_id)
        mu=w@ts.y.to_numpy();med=median(ts.y.to_numpy(),w)
        base=np.column_stack([x.loc[tq.parent_id].to_numpy()[:,:-1],tq[['log10_PHH_uM','log10_rat_mg_kg']]])
        predz=pr.loc[pr.outer_fold.eq(f)].pivot(index='parent_id',columns='geneset',values='proxy').loc[Q,z.columns]
        for name,val in [('training_compound_mean',mu),('training_weighted_median',med)]:
            r=tq.copy();r['route_model']=name;r['prediction']=val
            baselines.append(dict(outer_fold=f,reference=name,value=val,training_parents=len(S),query_parents=len(Q),training_ids=';'.join(S),**metrics(tq.parent_id,tq.y,np.full(len(tq),val))))
            if name=='training_weighted_median':newpred.append(r)
            else:assert np.allclose(old.loc[(old.outer_fold.eq(f))&old.route_model.eq(name),'prediction'],val,atol=1e-14,rtol=0)
        assert np.allclose(old.loc[(old.outer_fold.eq(f))&old.route_model.eq('zero_response'),'prediction'],0)
        for route in ['F_cond','F_struct','F_pred','F_meas']:
            obj=load_legacy_model(OLD/f'pilot/contexts/outer{f}_{route}.joblib');assert obj.train_ids==ts.parent_id.tolist()
            numeric=tq[['log10_PHH_uM','log10_rat_mg_kg']].to_numpy() if route=='F_cond' else base
            if route in ['F_pred','F_meas']:numeric=np.column_stack([numeric,(predz if route=='F_pred' else z).loc[tq.parent_id]])
            p=old_predict(obj,numeric,tq.route);savedold=old.loc[old.outer_fold.eq(f)&old.route_model.eq(route)].set_index('condition_id').loc[tq.condition_id]
            assert np.allclose(p,savedold.prediction,atol=1e-11,rtol=0)
            dose_index=ctx['numeric_feature_names'][route].index('log10_rat_mg_kg');slope=obj.model.coef_[dose_index]/obj.scaler.scale_[dose_index]
            for (pid,rt),r in savedold.groupby(['parent_id','route']):
                d=np.log10(r.dose_mg_kg.to_numpy());p=r.prediction.to_numpy();res=(p-p[0])-slope*(d-d[0]);checks.append(dict(outer_fold=f,route_model=route,parent_id=pid,administration_route=rt,n_conditions=len(r),dose_slope_from_object=slope,max_dose_slope_residual=abs(res).max()))
            modelrows.append(dict(outer_fold=f,route_model=route,actual_class=type(obj.model).__name__,n_numeric=obj.scaler.n_features_in_,alpha=obj.model.alpha,dose_index=dose_index,dose_slope=slope,max_saved_OOF_reconstruction_error=abs(p-savedold.prediction.to_numpy()).max() if len(p)==len(savedold) else float(abs(old_predict(obj,numeric,tq.route)-savedold.prediction.to_numpy()).max())))
            if route=='F_meas':
                meanz=w@z.loc[ts.parent_id].to_numpy();assert np.allclose(meanz,obj.scaler.mean_[-50:],atol=1e-12,rtol=0)
                ans=old_predict(obj,np.column_stack([base,np.tile(meanz,(len(tq),1))]),tq.route)
                r=tq.copy();r['route_model']='F_meanZ';r['prediction']=ans;newpred.append(r)
        for name in z.columns:
            for pid in Q:meanrefs.append(dict(parent_id=pid,outer_fold=f,geneset=name,training_mean=z.loc[S,name].mean()))
    csv(baselines,diag/'training_median_baseline.csv');csv(modelrows,diag/'actual_old_model_properties.csv')
    chk=pd.DataFrame(checks)
    wide=old.pivot(index='condition_id',columns='route_model',values='prediction');d=t[['condition_id','parent_id','source_name','outer_fold','y']].copy();d['e']=wide.loc[d.condition_id,'F_meas'].to_numpy()-d.y;d['a']=(wide.loc[d.condition_id,'F_swap']-wide.loc[d.condition_id,'F_meas']).to_numpy()
    d['two_e_a']=2*d.e*d.a;d['a_squared']=d.a**2;d['actual_squared_loss_difference']=(d.e+d.a)**2-d.e**2;d['identity_residual']=d.actual_squared_loss_difference-d.two_e_a-d.a_squared
    shifts=d.groupby('parent_id').agg(swap_shift_min=('a','min'),swap_shift_max=('a','max'));shifts['within_parent_swap_range']=shifts.swap_shift_max-shifts.swap_shift_min
    chk=chk.merge(shifts.reset_index(),on='parent_id');csv(chk,diag/'old_model_dose_and_swap_identity_checks.csv');csv(d,diag/'old_swap_squared_loss_decomposition.csv')
    csv([{k:float(q@d[k].to_numpy()) for k in ['two_e_a','a_squared','actual_squared_loss_difference','identity_residual']}],diag/'old_swap_squared_loss_decomposition_summary.csv')
    extended=pd.concat([old,*newpred],ignore_index=True);extended['error']=extended.prediction-extended.y;extended['absolute_error']=abs(extended.error);extended['squared_error']=extended.error**2;csv(extended,diag/'old_OOF_with_meanZ_and_median.csv')
    mets=[]
    for route,g in extended.groupby('route_model'):
        mets.append(dict(route_model=route,scope='pooled_OOF',**metrics(g.parent_id,g.y,g.prediction)))
        for f,h in g.groupby('outer_fold'):mets.append(dict(route_model=route,scope=f'fold_{int(f)}',**metrics(h.parent_id,h.y,h.prediction)))
    met=pd.DataFrame(mets);csv(met,diag/'old_metrics_with_ablations.csv')
    prior=pd.read_csv(OLD/'pilot/compound_equal_weight_metrics.csv');recon=prior.merge(met,on=['route_model','scope'],suffixes=('_saved','_recalculated'))
    for k in ['MAE_cmp','RMSE_cmp','R2_cmp']:recon[f'{k}_difference']=recon[k+'_recalculated']-recon[k+'_saved'];assert abs(recon[f'{k}_difference']).max()<1e-12
    csv(recon,diag/'old_metrics_reconciliation.csv')
    pr=pr.merge(pd.DataFrame(meanrefs),on=['parent_id','outer_fold','geneset'],validate='one_to_one');csv(pr,diag/'proxy_OOF_with_training_mean_reference.csv');pm=[]
    for name,g in pr.groupby('geneset'):
        y=g.measured.to_numpy();a=g.proxy.to_numpy();b=g.training_mean.to_numpy();v=np.var(y);sse=np.mean((a-y)**2);bse=np.mean((b-y)**2)
        pm.append(dict(geneset=name,n_compounds=len(g),MAE=np.mean(abs(a-y)),R2=1-sse/v,training_mean_MAE=np.mean(abs(b-y)),training_mean_R2=1-bse/v,MAE_minus_training_mean=np.mean(abs(a-y))-np.mean(abs(b-y)),MSE_skill_vs_training_mean=1-sse/bse))
    pm=pd.DataFrame(pm);csv(pm,diag/'proxy_all50_metrics_and_training_mean.csv')
    chkpm=pm.merge(proxy_metrics,on='geneset',suffixes=('_now','_saved'));assert abs(chkpm.MAE_now-chkpm.MAE_saved).max()<1e-12 and abs(chkpm.R2_now-chkpm.R2_saved).max()<1e-12
    assert chk.max_dose_slope_residual.max()<1e-10 and chk.within_parent_swap_range.max()<1e-10
    assert abs(ca.reconstruction_error).max()<1e-12 and abs(d.identity_residual).max()<1e-12
    dump(dict(completed_utc=now(),dev_only=True,source_parent_gene_files=len(ids),source_arrays=len(source),legacy_Ridge_objects_read=12,new_research_fits=0,holdout_diagnostics=0,**summary),diag/'completed.json')
    log('diagnostics_completed',**summary);print(summary);print(met.loc[met.scope.eq('pooled_OOF')].to_string(index=False))
