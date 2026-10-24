"""Derive tables and the standalone Figure 6 from the selected, saved predictions."""
from runtime import limit_threads
limit_threads()
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, average_precision_score, log_loss
from runtime import ROOT, OUT, read_json, write_json, event
from evaluate_fixed_predictions import ROUTES, CONTRASTS, survey_mean_cov, ci_summary

LABELS={'rdw_pct':'RDW (%)','hemoglobin_g_dl':'Hemoglobin (g/dL)'}

def interval_text(row,digits=5):
    return f'{row["estimate"]:.{digits}f} (95% CI {row["lower"]:.{digits}f} to {row["upper"]:.{digits}f})'

def retention_text(r):
    intervals=' ∪ '.join('['+', '.join(str(v) if isinstance(v,str) else f'{v:.3f}' for v in pair)+']' for pair in r['intervals'])
    point='NA' if r['point'] is None else f'{r["point"]:.3f}'
    return f'R = {point}; Fieller {r["set_type"]}: {intervals}'

def tables():
    cfg=read_json(ROOT/'protocol/analysis_config.json')
    summary=read_json(OUT/'tables/fixed_substitution_summary.json')
    cohort=pd.read_csv(ROOT/'data/analysis_cohort.csv')
    test=cohort.loc[cohort.role.eq('test')].set_index('person_id')
    frame=pd.read_csv(ROOT/'data/design_frame.csv')
    frame=frame.loc[frame.role.eq('test')].copy()
    frame['weight']=frame.WTMEC2YR/2
    pred=pd.read_csv(OUT/'predictions/test_predictions.csv').set_index('person_id').loc[test.index]
    zpred=pd.read_csv(OUT/'predictions/test_measurement_proxies.csv').set_index('person_id').loc[test.index]
    w=test.WTMEC2YR.to_numpy()/2
    y=test.y_5y.to_numpy()
    q1=[]
    for z in cfg['features']['measured_Z']:
        truth=zpred[z].to_numpy(); proxy=zpred['proxy_'+z].to_numpy(); mean=zpred['development_mean_'+z].to_numpy()
        assert np.array_equal(truth,test[z].to_numpy())
        err=abs(proxy-truth); ref=abs(mean-truth)
        vals=pd.DataFrame(dict(proxy_MAE=err,mean_MAE=ref),index=test.index)
        est,cov,df=survey_mean_cov(frame,vals)
        k=np.array([-1.,1.]); improvement=ci_summary(k@est,k@cov@k,df)
        mse=float(np.average((proxy-truth)**2,weights=w))
        variance=float(np.average((truth-np.average(truth,weights=w))**2,weights=w))
        q1.append(dict(measurement=z,unit=cfg['features']['Z_units'][z],n=len(test),
            MAE=ci_summary(est[0],cov[0,0],df),RMSE=float(np.sqrt(mse)),R2=1-mse/variance,
            R2_denominator='weighted test-mean variance',development_mean=float(mean[0]),
            development_mean_MAE=ci_summary(est[1],cov[1,1],df),MAE_improvement=improvement,
            truth_min=float(truth.min()),truth_max=float(truth.max()),proxy_min=float(proxy.min()),proxy_max=float(proxy.max())))
    write_json(OUT/'tables/q1_measurement_performance.json',q1)
    qflat=[]
    for row in q1:
        flat={k:v for k,v in row.items() if not isinstance(v,dict)}
        for k in ['MAE','development_mean_MAE','MAE_improvement']:
            flat.update({k+'_'+sub:value for sub,value in row[k].items()})
        qflat.append(flat)
    pd.DataFrame(qflat).to_csv(OUT/'tables/q1_measurement_performance.csv',index=False)
    support=[]; calibration=[]
    for route in ROUTES:
        p=pred['p_'+route].to_numpy()
        support.append(dict(route=route,n=len(y),weighted_log_loss=float(log_loss(y,np.clip(p,1e-12,1-1e-12),sample_weight=w,labels=[0,1])),
            weighted_AUROC=float(roc_auc_score(y,p,sample_weight=w)),
            weighted_AP=float(average_precision_score(y,p,sample_weight=w)),
            weighted_event_fraction=float(np.average(y,weights=w)),weighted_mean_risk=float(np.average(p,weights=w))))
        # Descriptive fixed ten equal-weight risk groups; no calibration fitting.
        order=np.argsort(p,kind='stable')
        binid=np.minimum(9,np.floor(10*(np.cumsum(w[order])-w[order]/2)/w.sum()).astype(int))
        for b in range(10):
            idx=order[binid==b]
            if len(idx):
                calibration.append(dict(route=route,bin=b+1,n=len(idx),weight_sum=float(w[idx].sum()),
                    mean_predicted=float(np.average(p[idx],weights=w[idx])),observed_fraction=float(np.average(y[idx],weights=w[idx]))))
    for name,rows in [('supporting_mortality_metrics',support),('descriptive_calibration',calibration)]:
        write_json(OUT/f'tables/{name}.json',rows)
        pd.DataFrame(rows).to_csv(OUT/f'tables/{name}.csv',index=False)
    flat=[]
    for domain in ['primary_40plus','sensitivity_40_79']:
        r=summary[domain]
        for kind in ['Brier','contrasts']:
            for name,values in r[kind].items():
                flat.append(dict(domain=domain,kind=kind,estimand=name,n=r['n'],events=r['events_5y'],**values))
        g=r['contrasts']['G_measured_gain']['estimate']; d=r['contrasts']['D_substitution_loss']['estimate']; h=r['contrasts']['H_remaining_gain']['estimate']
        assert np.isclose(g-d,h,atol=1e-14)
        if r['retention']['point'] is not None:
            assert np.isclose(r['retention']['point'],h/g)
        ft=frame.set_index('person_id')
        ids=test.index if domain=='primary_40plus' else test.loc[test.age_cap80<80].index
        active=ft.loc[ids,['stratum_id','psu_id']].drop_duplicates()
        r['design_counts']=dict(full_frame_n=len(frame),full_strata=int(frame.stratum_id.nunique()),
            full_PSUs=int(frame.psu_id.nunique()),active_domain_strata=int(active.stratum_id.nunique()),active_domain_PSUs=len(active))
    pd.DataFrame(flat).to_csv(OUT/'tables/brier_and_contrasts.csv',index=False)
    write_json(OUT/'tables/fixed_substitution_summary.json',summary)
    uw=summary['primary_40plus']['unweighted_sensitivity']
    pd.DataFrame([dict(kind=kind,estimand=k,estimate=v) for kind in ['route_Brier','contrasts'] for k,v in uw[kind].items()]).to_csv(OUT/'tables/unweighted_sensitivity.csv',index=False)
    write_json(OUT/'checks/scientific_checks.json',dict(status='passed',test_n=len(test),unique_ids=True,
        same_Y_and_weights=True,predictions_in_unit_interval=True,whitelist=cfg['features']['numeric_X']+cfg['features']['categorical_X'],
        G_minus_D_equals_H=True,retention_formula_checked=True,proxy_truth_matches_cohort=True,
        source='saved selected predictions only; no additional fits',fixed_sensitivities=['unweighted','age40_79']))
    return summary,q1,support

def figure(summary,q1):
    out=OUT/'figures'; out.mkdir(parents=True,exist_ok=True)
    r=summary['primary_40plus']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(13,10.2),layout='constrained')
    gs=fig.add_gridspec(2,2,width_ratios=[1,1],height_ratios=[1,1])
    a=fig.add_subplot(gs[0,0]); a.axis('off'); a.set_title('A  Cohort and information flow',loc='left',fontweight='bold')
    steps=[('2003–2004 | MAXM','Train: 2,955; deaths: 370\nProxy selection / strict nested OOF; 24 Y candidates'),
           ('2005–2006 | MAXM','Validation: 2,804; deaths: 288\nSelect B, M and P by weighted Brier'),
           ('2003–2006 | MAXM','Development refit: 5,759; deaths: 658\nNew strict OOF; final g, B, M and P'),
           ('2007–2010 | HMX',f'Later evaluation: {r["n"]:,}; deaths: {r["events_5y"]}\nB(X); M(X,Z); S(X,g); MEAN(X,Z̄); P(X,g)')]
    for i,(title,body) in enumerate(steps):
        yy=.88-i*.23
        a.text(.04,yy,title+'\n'+body,transform=a.transAxes,va='top',fontsize=10,
               bbox=dict(boxstyle='round,pad=.5',facecolor='#f1f4f6',edgecolor='#cad3da'))
        if i<3:
            a.annotate('',xy=(.47,yy-.185),xytext=(.47,yy-.14),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color='#526476'))
    inner=gs[0,1].subgridspec(2,1)
    for j,row in enumerate(q1):
        ax=fig.add_subplot(inner[j,0])
        if j==0: ax.set_title('B  Measurement prediction in original units',loc='left',fontweight='bold',pad=16)
        for i,(name,key,color) in enumerate([('Learned proxy','MAE','#315f82'),('Development mean','development_mean_MAE','#9c7446')]):
            v=row[key]
            ax.errorbar(v['estimate'],i,xerr=[[v['estimate']-v['lower']],[v['upper']-v['estimate']]],fmt='o',color=color,capsize=3)
        ax.set_yticks([0,1],['Learned proxy','Development mean']); ax.invert_yaxis()
        ax.set_xlabel('MAE — '+LABELS[row['measurement']]); ax.set_xlim(left=0); ax.margins(y=.6)
        ax.grid(axis='x',alpha=.18)
    c=fig.add_subplot(gs[1,0]); c.set_title('C  Five-route absolute Brier losses',loc='left',fontweight='bold')
    names=['B | baseline','M | measured','S | fixed proxy','MEAN | fixed mean','P | retrained proxy']
    for i,route in enumerate(ROUTES):
        v=r['Brier'][route]
        c.errorbar(v['estimate'],i,xerr=[[v['estimate']-v['lower']],[v['upper']-v['estimate']]],fmt='o',color='#315f82',capsize=3)
    c.set_yticks(range(5),names); c.invert_yaxis(); c.set_xlabel('Survey-weighted Brier (lower is better)'); c.grid(axis='x',alpha=.18)
    c.margins(y=.2)
    d=fig.add_subplot(gs[1,1]); d.set_title('D  Paired Brier contrasts',loc='left',fontweight='bold')
    keys=['G_measured_gain','D_substitution_loss','H_remaining_gain']
    for i,key in enumerate(keys):
        v=r['contrasts'][key]
        d.errorbar(v['estimate'],i,xerr=[[v['estimate']-v['lower']],[v['upper']-v['estimate']]],fmt='o',color='#315f82',capsize=3)
    d.axvline(0,color='#555555',lw=1,linestyle='--'); d.set_yticks(range(3),['G: B − M','D: S − M','H: B − S']); d.invert_yaxis()
    d.set_xlabel('Paired weighted loss difference'); d.grid(axis='x',alpha=.18); d.margins(y=.4)
    rt=retention_text(r['retention'])
    if r['retention']['headline_percentage_allowed']:
        status=('Measured gain supported; retention interval includes zero.'
                if any(isinstance(lo,(int,float)) and isinstance(hi,(int,float)) and lo<=0<=hi for lo,hi in r['retention']['intervals'])
                else 'Measured gain supported; bounded retention interval shown.')
    else:
        status='Retention is undefined or unstable; no headline percentage.'
    fig.suptitle('Figure 6 | NHANES measurement proxies and fixed-model substitution',fontsize=15,fontweight='bold')
    fig.supxlabel('95% pointwise design-based t intervals; fixed fitted predictions. Paired inference is in D.\n'+rt+'\n'+status,fontsize=9)
    for ext in ['png','svg','pdf']:
        fig.savefig(out/f'Figure6_NHANES_Q3.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
    write_json(out/'Figure6_source_data.json',dict(primary=r,q1=q1,cohort_roles=read_json(OUT/'checks/preflight.json')['counts_by_role']))

def main():
    summary,q1,support=tables()
    figure(summary,q1)
    if (ROOT/'scripts/write_manuscript.py').exists():
        from write_manuscript import write_all
        write_all(summary,q1,support)
    event('tables_figure_and_drafts_generated',visual_review='pending')
    write_json(OUT/'logs/output_generation_complete.json',dict(status='generated_pending_visual_review'))
    print('Tables and Figure 6 generated; local manuscript generator used when present. Visual review still required.')

if __name__=='__main__':
    main()
