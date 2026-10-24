"""new_random=yes fixed-OOF group bootstrap; model_fit=no; model_query=no; used_evaluation_data=yes post-hoc O+N scores."""
from common import *
from scipy.stats import t as student_t
from sklearn.metrics import average_precision_score,roc_auc_score
E=P/'pooled_cv';ROLE='post-hoc used-domain grouped CV; not independent N validation'
CONTRASTS=[('M1','M0'),('M1','Mperm'),('M1','Mchem')]

def weighted_ap_many(y,score,weights):
    """Exactly the weighted step AP, including ties; rows are original cluster draws."""
    order=np.argsort(-score,kind='stable');s=score[order];yy=y[order];w=weights[:,order]
    ends=np.r_[np.flatnonzero(np.diff(s)),len(s)-1];tp=np.cumsum(w*yy,axis=1)[:,ends];den=np.cumsum(w,axis=1)[:,ends]
    inc=np.diff(np.column_stack([np.zeros(len(w)),tp]),axis=1);positive=tp[:,-1]
    precision=np.divide(tp,den,out=np.zeros_like(tp,dtype=float),where=den>0)
    result=np.divide(np.sum(inc*precision,axis=1),positive,out=np.full(len(w),np.nan),where=positive>0)
    total=w.sum(axis=1);result[(positive==0)|(positive==total)]=np.nan
    return result

def holm(p):
    p=np.asarray(p,float);out=np.full(len(p),np.nan);ok=np.flatnonzero(np.isfinite(p));order=ok[np.argsort(p[ok])];current=0.
    for rank,i in enumerate(order):current=max(current,min(1.,(len(p)-rank)*p[i]));out[i]=current
    return out

def run_statistics(m,cohorts,completed_b):
    d=readcsv(E/'pooled_oof_predictions.csv');models=['M0','M1','Mperm','Mchem','L2LR'];perrep=[];perfold=[];comparisons=[]
    for (cohort,rep,model),g in d.groupby(['cohort','repeat','model']):
        expected=set(m.member_id.iloc[cohorts[cohort]['rows']]);assert g.member_id.is_unique and set(g.member_id)==expected
        y=g.y.to_numpy(int);p=g.raw_score.to_numpy();c=g.calibrated_score.to_numpy();row=dict(cohort=cohort,repeat=rep,model=model,n=len(g),positive=int(y.sum()),AP=float(average_precision_score(y,p)),AP_minus_prevalence=float(average_precision_score(y,p)-y.mean()),ROC_AUC=float(roc_auc_score(y,p)),Brier_raw=float(np.mean((y-p)**2)),Brier_calibrated=float(np.mean((y-c)**2)) if np.isfinite(c).all() else np.nan,evidence_role=ROLE);perrep.append(row)
        for f,h in g.groupby('outer_fold'):
            yy=h.y.to_numpy(int);pp=h.raw_score.to_numpy();cc=h.calibrated_score.to_numpy();perfold.append(dict(cohort=cohort,repeat=rep,outer_fold=f,model=model,n_test=len(h),n_train=len(g)-len(h),query_groups=h.fold_group.nunique(),training_groups=g.loc[g.outer_fold.ne(f),'fold_group'].nunique(),AP=float(average_precision_score(yy,pp)),ROC_AUC=float(roc_auc_score(yy,pp)),Brier_raw=float(np.mean((yy-pp)**2)),Brier_calibrated=float(np.mean((yy-cc)**2)) if np.isfinite(cc).all() else np.nan))
    rep=pd.DataFrame(perrep);fold=pd.DataFrame(perfold);savecsv(rep,E/'repeat_pooled_metrics.csv',ROLE);savecsv(fold,E/'pooled_cv_results.csv',ROLE)
    for cohort,rr in rep.groupby('cohort'):
        for a,b in CONTRASTS:
            aa=rr.loc[rr.model.eq(a)].set_index('repeat');bb=rr.loc[rr.model.eq(b)].set_index('repeat').loc[aa.index]
            for r in aa.index:comparisons.append(dict(cohort=cohort,repeat=r,comparison=f'{a}-{b}',delta_AP=aa.loc[r,'AP']-bb.loc[r,'AP'],delta_AUC=aa.loc[r,'ROC_AUC']-bb.loc[r,'ROC_AUC'],delta_Brier_raw=aa.loc[r,'Brier_raw']-bb.loc[r,'Brier_raw'],delta_Brier_calibrated=aa.loc[r,'Brier_calibrated']-bb.loc[r,'Brier_calibrated']))
    comps=pd.DataFrame(comparisons);savecsv(comps,E/'repeat_comparisons.csv',ROLE)
    summary=rep.groupby(['cohort','model']).agg(mean_AP=('AP','mean'),min_AP=('AP','min'),max_AP=('AP','max'),mean_AUC=('ROC_AUC','mean'),mean_Brier_raw=('Brier_raw','mean'),mean_Brier_calibrated=('Brier_calibrated','mean')).reset_index();savecsv(summary,E/'repeat_metric_summary.csv','ten correlated split repetitions; range not ten independent datasets')
    nb=[]
    for cohort,ff in fold.groupby('cohort'):
        local=[]
        for a,b in CONTRASTS:
            aa=ff.loc[ff.model.eq(a)].set_index(['repeat','outer_fold']);bb=ff.loc[ff.model.eq(b)].set_index(['repeat','outer_fold']).loc[aa.index];delta=(aa.AP-bb.AP).to_numpy();n=len(delta);assert n==50
            mean=delta.mean();variance=delta.var(ddof=1);ratio=(aa.n_test/aa.n_train).mean();se=np.sqrt((1/n+ratio)*variance)
            tstat=mean/se if se>0 else np.nan;p=float(2*student_t.sf(abs(tstat),n-1)) if se>0 else np.nan;critical=student_t.ppf(.975,n-1)
            local.append(dict(cohort=cohort,comparison=f'{a}-{b}',estimand='arithmetic mean of 50 fold-specific AP differences',mean_fold_delta_AP=mean,fold_variance=variance,mean_n_test_over_n_train=ratio,m=n,df=n-1,corrected_SE=se,t=tstat,p_approximate=p,low_95=mean-critical*se if se>0 else np.nan,high_95=mean+critical*se if se>0 else np.nan,NA_reason='' if se>0 else 'zero fold-difference variance; all-zero is not equivalence',limitations='Nadeau-Bengio approximation; unequal scaffold groups, stratification and nonlinear AP do not ensure exact type-I control'))
        adjusted=holm([x['p_approximate'] for x in local])
        for row,p in zip(local,adjusted):row['p_Holm_three_AP_contrasts']=p
        nb.extend(local)
    savecsv(nb,E/'Nadeau_Bengio_corrected_AP.csv',ROLE)
    intervals=[];bootstrap_support=[];all_draws=[]
    for cohort in cohorts:
        dd=d.loc[d.cohort.eq(cohort)];ids=sorted(dd.member_id.unique());base=dd.loc[(dd.repeat==0)&(dd.model=='M0')].set_index('member_id').loc[ids];groups=sorted(base.fold_group.unique());group_index={g:i for i,g in enumerate(groups)};codes=np.array([group_index[g] for g in base.fold_group]);y=base.y.to_numpy(int);drawpath=E/f'{cohort}_group_bootstrap_indices.npz'
        if drawpath.exists():
            z=np.load(drawpath);assert np.array_equal(z['groups'],groups);indices=z['indices'];event('P3','group_bootstrap_indices_reused',cohort=cohort)
        else:
            indices=np.random.Generator(np.random.PCG64(61034001)).integers(0,len(groups),size=(2000,len(groups)),dtype=np.int32);np.savez_compressed(drawpath,groups=np.array(groups),indices=indices,seed=61034001);event('P3','new_group_bootstrap_draws',cohort=cohort,draws=2000,groups_per_draw=len(groups),seed=61034001,models_refit=0)
        counts=np.zeros((2000,len(groups)),dtype=np.int32)
        for i in range(2000):counts[i]=np.bincount(indices[i],minlength=len(groups))
        weights=counts[:,codes];tot=weights.sum(1);pi=(weights@y)/tot;valid=(pi>0)&(pi<1);group_sizes=np.bincount(codes,minlength=len(groups));maximum=(counts*group_sizes).max(1)/tot
        for i in range(2000):bootstrap_support.append(dict(cohort=cohort,draw=i,total_identity_weight=int(tot[i]),positive_fraction=float(pi[i]),maximum_group_fraction=float(maximum[i]),valid=bool(valid[i]),invalid_reason='' if valid[i] else 'single-class draw; no replacement draw'))
        metric_arrays={model:np.zeros(2000) for model in models}
        for model in models:
            arrays=[]
            for r in range(10):
                h=dd.loc[dd.model.eq(model)&dd.repeat.eq(r)].set_index('member_id').loc[ids];assert np.array_equal(h.y,y)
                arrays.append(weighted_ap_many(y,h.raw_score.to_numpy(),weights))
            metric_arrays[model]=np.mean(arrays,axis=0)
        for a,b in CONTRASTS:
            delta=metric_arrays[a]-metric_arrays[b];vals=delta[np.isfinite(delta)];point=float(comps.loc[comps.cohort.eq(cohort)&comps.comparison.eq(f'{a}-{b}'),'delta_AP'].mean())
            lo,hi=np.quantile(vals,[.025,.975],method='linear') if len(vals) else [np.nan,np.nan]
            intervals.append(dict(cohort=cohort,comparison=f'{a}-{b}',estimand='mean of ten repeat-specific pooled OOF AP differences',point=point,low_95=lo,high_95=hi,valid_draws=len(vals),invalid_draws=2000-len(vals),unit='whole Murcko group; all identities and ten predictions retained',conditions='fixed observed OOF and groups; no retraining uncertainty',evidence_role=ROLE))
            all_draws.extend(dict(cohort=cohort,comparison=f'{a}-{b}',draw=i,delta_AP=float(v)) for i,v in enumerate(delta))
    savecsv(intervals,E/'fixed_OOF_group_bootstrap_intervals.csv',ROLE);savecsv(bootstrap_support,E/'group_bootstrap_support.csv',ROLE);savecsv(all_draws,E/'group_bootstrap_differences.csv',ROLE)
    perm=[];base=rep.loc[rep.cohort.eq('P_mixed_full')&rep.model.eq('M0')].set_index('repeat');real=float(comps.loc[comps.cohort.eq('P_mixed_full')&comps.comparison.eq('M1-M0'),'delta_AP'].mean())
    for number in completed_b:
        h=readcsv(E/f'permutation_{number:03d}_oof.csv');delta=[]
        for r,g in h.groupby('repeat'):
            assert g.member_id.is_unique and len(g)==750;delta.append(float(average_precision_score(g.y,g.raw_score)-base.loc[r,'AP']))
        assert len(delta)==10;perm.append(dict(b=number,T_b=float(np.mean(delta)),repeat_delta_AP=json.dumps(delta),evidence_role='Monte Carlo random auxiliary supervision reference; not confirmatory clinical increment null'))
    savecsv(pd.DataFrame(perm,columns=['b','T_b','repeat_delta_AP','evidence_role']),E/'permutation_null.csv',ROLE)
    vals=np.array([r['T_b'] for r in perm]);B=len(vals);tail=int(np.sum(vals>=real));ref=dict(T_real=real,B_complete=B,B_expected=100,tail_count_ge=tail,reference_tail_fraction=(1+tail)/(B+1) if B else None,reference_fraction_less=float(np.mean(vals<real)) if B else None,minimum_attainable_for_B100=1/101,status='complete reference' if B==100 else 'incomplete reference; do not present as 100',interpretation='Only a Monte Carlo reference under task-wise A-label exchangeability; chemical/source and cross-task dependencies limit an incremental-null interpretation; not a bootstrap sign P value.')
    ref.update(reference_mean=float(vals.mean()) if B else None,reference_sample_SD=float(vals.std(ddof=1)) if B>1 else None,
               empirical_quantiles={str(q):float(np.quantile(vals,q,method='linear')) for q in [0,.025,.25,.5,.75,.975,1]} if B else {},
               empirical_quantile_method='linear; descriptive finite reference, not effect confidence limits')
    dump(ref,E/'permutation_reference_summary.json');plot(rep,comps,vals,real,E)
    log=json.loads((E/'run_log.json').read_text('utf8'));log.update(status='complete' if B==100 else 'partial_main_CV_and_statistics_complete',statistics_completed_utc=now(),bootstrap_draws_per_cohort=2000,reference=ref);dump(log,E/'run_log.json')
    with (E/'run_report.md').open('a',encoding='utf8') as f:f.write(f'\n三集合10次拼接OOF、50折近似校正t（同队列三AP比较Holm）和原结构整簇2000次固定OOF抽样已生成。置换完整数{B}/100，尾部计数{tail}，参考比例{ref["reference_tail_fraction"]}。校正t目标是平均逐折差；整簇区间目标是平均repeat拼接AP差，不互换估计对象，也不算独立证据。\n')

def plot(rep,comps,vals,real,E):
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(1,3,figsize=(170/25.4,96/25.4));models=['M0','M1','Mperm','Mchem','L2LR'];h=rep.loc[rep.cohort.eq('P_mixed_full')]
    axes[0].boxplot([h.loc[h.model.eq(n),'AP'] for n in models],tick_labels=models);axes[0].tick_params(axis='x',rotation=60);axes[0].set_ylabel('Repeat-pooled AP');axes[0].set_title('A  Post-hoc\nmixed cohort',fontweight='bold',fontsize=8)
    if len(vals):axes[1].hist(vals,bins=15,color='#6f7780',alpha=.7)
    axes[1].axvline(real,color='#21618c',lw=1.5,label='Real ΔAP');axes[1].set(xlabel='Mean repeat\nΔAP versus M0',ylabel='Reference sets');axes[1].set_title(f'B  Auxiliary permutations\n({len(vals)}/100 joint sets)',fontweight='bold',fontsize=8)
    axes[1].legend(loc='upper left',fontsize=7,frameon=False,handlelength=1)
    names=['P_mixed_full','P_mixed_common','P_v2_common']
    for i,name in enumerate(names):
        v=comps.loc[comps.cohort.eq(name)&comps.comparison.eq('M1-M0'),'delta_AP'];axes[2].scatter(np.full(len(v),i),v,s=12,color='#21618c');axes[2].plot([i-.2,i+.2],[v.mean()]*2,color='#b55d24')
    axes[2].set(xticks=range(3),xticklabels=['Full mixed','Common mixed','Common V2'],ylabel='Repeat-pooled ΔAP M1−M0');axes[2].set_title('C  Label/member\nsensitivity',fontweight='bold',fontsize=8);axes[2].tick_params(axis='x',rotation=55);axes[2].axhline(0,color='.5',lw=.6)
    fig.tight_layout()
    for ext in ['pdf','svg','png']:fig.savefig(E/f'Fig_S4_posthoc_pooled.{ext}',dpi=300)
    plt.close(fig)
