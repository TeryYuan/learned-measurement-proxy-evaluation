"""Revised six-figure set from immutable saved aggregates and licensed SVGs.

Usage: python source/render_figures.py --figures 1 4 6 2 3 5
No raw records, estimators, bootstrap or training code are used.
"""
import argparse,csv,json
from render_common import *

def fig1():
    fig,a=canvas(540)
    icon(fig,'liver',14,13,30,23);icon(fig,'cell_culture_dish',49,15,24,20)
    txt(a,81,25,'Liver-safety applications',7.5)
    icon(fig,'adult_group',274,13,25,25);icon(fig,'blood_sample_tube',304,13,18,26)
    txt(a,330,25,'Clinical measurements',7.5)
    line(a,12,48,481,48)
    for x,label in [(14,'QUESTION'),(108,'COMPARISON'),(398,'EVIDENCE')]:
        txt(a,x,65,label,7,color=GRAY,weight='bold')
    for y in (207,348,507): line(a,12,y,481,y)
    for y1,y2 in [(82,197),(221,338),(362,497)]:line(a,389,y1,389,y2)

    txt(a,14,98,'Q1',10,weight='bold')
    txt(a,14,130,'Predict the\nmeasurement',8,weight='bold')
    txt(a,108,99,'Measured Z',7.5,BLUE)
    txt(a,182,99,'training only',7,GRAY)
    a.annotate('',xy=(199,137),xytext=(150,105),arrowprops=dict(
        arrowstyle='-|>',ls='--',color=BLUE,lw=.8,mutation_scale=7,
        connectionstyle='angle,angleA=0,angleB=-90,rad=0'))
    txt(a,108,154,'Allowed X',7.5)
    arrow(a,(151,154),(165,154))
    card(a,168,135,65,39,'Learned proxy\ng(X)',ORANGE,PALEORANGE)
    arrow(a,(234,154),(248,154))
    card(a,250,135,91,28,'Proxy Ẑ',ORANGE,PALEORANGE,color=ORANGE)
    card(a,250,175,91,22,'Observed Z',BLUE,PALEBLUE,color=BLUE)
    line(a,347,149,354,149,GRAY);line(a,354,149,354,186,GRAY)
    line(a,347,186,354,186,GRAY);txt(a,361,168,'vs',7,GRAY)
    txt(a,108,188,'BSEP; RDW / Hb',7,GRAY)
    txt(a,398,130,'Prediction error\n/ discrimination',7.5,weight='bold')
    txt(a,398,173,'Proxy versus\nobserved value',7,GRAY)

    txt(a,14,231,'Q2',10,weight='bold')
    txt(a,14,269,'Improve\ndownstream\nprediction',8,weight='bold')
    txt(a,108,226,'Two separately fitted models',7.5,weight='bold')
    for y,inp,label,col in [(250,'X','Baseline B',GRAY),(290,'X + Ẑ','Proxy-input P',ORANGE)]:
        txt(a,108,y+13,inp,7.5,color=col)
        arrow(a,(151,y+13),(170,y+13),col)
        card(a,173,y,88,27,label,col)
        arrow(a,(263,y+13),(282,y+13),col)
        txt(a,288,y+13,'Prediction',7.5,col)
    line(a,344,263,351,263,GRAY);line(a,351,263,351,303,GRAY)
    line(a,344,303,351,303,GRAY);txt(a,359,283,'vs',7,GRAY)
    txt(a,108,334,'DILI concern; NHANES mortality',7,GRAY)
    txt(a,398,259,'Increment\nover baseline',7.5,weight='bold')
    txt(a,398,305,'Same cases\nand outcome Y',7,GRAY)

    txt(a,14,374,'Q3',10,weight='bold')
    txt(a,14,416,'Evaluate\ninput\nreplacement',8,weight='bold')
    txt(a,108,369,'First establish measured-input use and benefit',7.5,weight='bold')
    txt(a,108,391,'Two queries of the same fitted model',7.5)
    card(a,108,411,120,63,edge=SEP)
    txt(a,114,424,'Query state',7,color=GRAY)
    txt(a,114,443,'X + measured Z',7.5,BLUE)
    txt(a,114,462,'or X + proxy Ẑ',7.5,ORANGE)
    arrow(a,(230,445),(245,445))
    card(a,247,419,73,53,'Fitted M\nunchanged',BLUE,PALEBLUE,weight='bold')
    arrow(a,(322,445),(338,445))
    txt(a,344,443,'Predictions:\nM / S',7.5)
    txt(a,108,492,'PHH–rat ALT / constructed control; NHANES',7,GRAY)
    txt(a,398,394,'Measured gain',7.5,BLUE,weight='bold')
    txt(a,398,409,r'$G = L_B - L_M$',7)
    txt(a,398,435,'Substitution loss',7.5,ORANGE,weight='bold')
    txt(a,398,450,r'$D = L_S - L_M$',7)
    txt(a,398,476,'Remaining gain',7.5,TEAL,weight='bold')
    txt(a,398,491,r'$H = L_B - L_S$',7)
    txt(a,W/2,525,'Establish measured-input use and benefit before interpreting replacement.',7.5,
        weight='bold',ha='center')
    save(fig,'Fig1_three_questions')

def saved():
    return json.loads((ROOT/'source/data/plotted_saved_values.json').read_text(encoding='utf-8'))

def fig4():
    data=saved()['Figure4']
    fig,a=canvas(539)
    panel(a,'A','Input removal makes substitution loss zero',12,17)
    # Candidate inputs and retained query path. Muted branches stop before f(c).
    card(a,15,37,22,20,'c',GRAY,size=8)
    txt(a,43,43,'Experimental conditions',7.5,weight='bold')
    txt(a,43,56,'PHH concentration, rat dose / route',7,GRAY)
    arrow(a,(188,47),(228,47))
    card(a,232,35,102,58,edge=GRAY,fill='#F7F9FA')
    txt(a,283,48,'Selected function',7,GRAY,ha='center')
    txt(a,283,68,r'$f(c)=b+\beta_C^T C(c)$',8,ha='center')
    txt(a,283,85,'Conditions only',7.5,weight='bold',ha='center')
    card(a,15,69,22,18,'X',SEP,color=MUTED)
    txt(a,43,78,'Molecular structure',7.5,MUTED)
    line(a,144,78,175,78,SEP,.9,ls='--');line(a,175,73,175,83,GRAY,.7)
    txt(a,179,78,'not selected',7,MUTED)
    icon(fig,'cell_culture_dish',14,96,24,18)
    txt(a,43,100,'Measured PHH responses',7.5,MUTED)
    line(a,164,102,175,102,SEP,.9,ls='--');line(a,175,97,175,107,GRAY,.7)
    txt(a,179,103,'not selected',7,MUTED)
    txt(a,43,114,'24 h; 50 Hallmark member means',7,GRAY)
    arrow(a,(337,63),(350,63))
    txt(a,356,43,'Same predictions',7.5,weight='bold')
    txt(a,356,61,'B = M = S = f(c)',8)
    txt(a,356,83,'G = H = D = 0',8,weight='bold')
    txt(a,356,103,'R = H/G: undefined',7.5)
    txt(a,15,138,'Matched compound / condition records',7,GRAY)
    icon(fig,'rat',221,117,47,29)
    txt(a,273,129,'Public rat ALT data',7.5,weight='bold')
    txt(a,273,142,'Existing single-dose, 24 h biochemistry',7,GRAY)
    line(a,12,151,481,151)
    txt(a,12,164,'Holdout: 24 parent compounds / 72 conditions',7.5,weight='bold')
    txt(a,12,176,'Development: 97 parent compounds',7,GRAY)

    panel(a,'B','Paired ALT comparisons',12,192)
    b=chart(fig,205,215,273,101)
    keys=['M1_measured','M2_predicted','M3_substitution']
    for i,key in enumerate(keys):
        for cov,lw in [(.9833333333333333,.8),(.95,2.1)]:
            row=next(r for r in data['paired'] if r['comparison']==key and float(r['confidence'])==cov)
            interval(b,row,i,('estimate','ci_low','ci_high'),BLUE if i==0 else ORANGE,
                     ['o','D','s'][i],lw=lw,tag=f'Figure4A:{key}:{cov}',open_marker=i==2)
    b.set(yticks=range(3),yticklabels=['Nominal measured − reference','Independent proxy − reference',
          'Fixed swap − measured'],ylim=(2.6,-.6),xlim=(-.013,.009),
          xlabel='Paired MAE difference',xticks=[-.01,-.005,0,.005])
    b.axvline(0,color=GRAY,ls=':',lw=.8)
    for y in (0,2):
        b.annotate('= 0 by construction',xy=(0,y),xytext=(7,0),textcoords='offset points',
                   fontsize=7,va='center',ha='left',color=INK)
    txt(a,205,204,'95% thick; 98.33% thin',7,GRAY)
    panel(a,'C','Absolute ALT error',12,362)
    c=chart(fig,205,384,273,119)
    labs=['Condition-only (6 identical routes)','Independent predicted input','Permutation-trained input',
          'Zero response','Training compound mean','Training weighted median']
    for i,row in enumerate(data['absolute']):
        interval(c,row,i,('estimate','ci_low','ci_high'),ORANGE if i==1 else GRAY,
                 ['s','D','x','v','^','<'][i],tag=f"Figure4B:{row['route_model']}",open_marker=i in (4,5))
    c.set(yticks=range(6),yticklabels=labs,ylim=(5.6,-.6),xlim=(.11,.27),
          xlabel='Parent-equal MAE (ALT log2-response units)',xticks=[.12,.16,.20,.24])
    grid(c);txt(a,205,374,'Single-model 95% intervals',7,GRAY)
    save(fig,'Fig4_NAM_case')

def fig6():
    source=json.loads((ROOT/'source/data/Figure6_source_data.json').read_text(encoding='utf-8'))
    r=source['primary'];q=source['q1'];roles=source['cohort_roles']
    fig,a=canvas(590)
    panel(a,'A','Temporal roles in NHANES',12,16)
    icon(fig,'adult_group',14,34,24,25)
    txt(a,44,46,'Adults ≥40 years',7.5,weight='bold')
    icon(fig,'questionnaire',175,34,21,25)
    txt(a,203,41,'X: 16 baseline variables',7.5)
    txt(a,203,54,'Survey / examination',7,GRAY)
    icon(fig,'blood_sample_tube',320,30,18,28)
    txt(a,345,41,'Measured Z',7.5,BLUE)
    txt(a,345,54,'RDW / hemoglobin',7,GRAY)

    for x,w,key,period,label in [(15,123,'train','2003–2004','Fit'),
                                 (156,123,'validation','2005–2006','Select'),
                                 (349,130,'test','2007–2010','Evaluate')]:
        rr=roles[key]
        card(a,x,77,w,64,edge=SEP)
        txt(a,x+10,91,period,8,weight='bold')
        txt(a,x+10,109,label,8,weight='bold',color=BLUE if key!='test' else INK)
        txt(a,x+10,128,f"n={rr['n']:,}; events={rr['events']:,}",7)
    txt(a,306,67,'MAXM → HMX instrument change',7,GRAY,ha='center')
    line(a,306,77,306,153,GRAY,.75,ls=':')
    line(a,306,191,306,209,GRAY,.75,ls=':')
    # A pooling bracket, not an additional population or reverse time arrow.
    line(a,15,145,15,150,GRAY);line(a,15,150,279,150,GRAY);line(a,279,145,279,150,GRAY)
    line(a,147,150,147,157,GRAY)
    card(a,15,158,264,48,edge=GRAY,fill='#F7F9FA')
    txt(a,25,170,'2003–2006   Refit / final development',8,weight='bold')
    n=roles['train']['n']+roles['validation']['n'];ev=roles['train']['events']+roles['validation']['events']
    txt(a,25,186,f'n={n:,}; events={ev:,}',7)
    txt(a,25,198,'Pooled first two periods',7,GRAY)
    card(a,349,158,130,48,edge=SEP)
    arrow(a,(279,183),(349,183),GRAY)
    txt(a,314,174,'Fitted models',7,GRAY,ha='center')
    txt(a,359,172,'Later-cycle evaluation',7.5,weight='bold')
    txt(a,359,191,'Same cases for all five routes',7,GRAY)
    line(a,12,217,481,217)

    panel(a,'B','Measurement prediction',12,232)
    for j,row in enumerate(q):
        x=104 if j==0 else 352
        title='RDW' if j==0 else 'Hemoglobin'
        txt(a,x,252,f"{title}  ·  R²={row['R2']:.3f}",7.5,weight='bold')
        ax=chart(fig,x,265,126,57)
        for i,(key,col,mark) in enumerate([('MAE',ORANGE,'o'),('development_mean_MAE',GRAY,'^')]):
            interval(ax,row[key],i,('estimate','lower','upper'),col,mark,
                     tag=f"Figure6B:{row['measurement']}:{key}",open_marker=i==1)
        ax.set(yticks=[0,1],yticklabels=['Proxy','Measurement mean'],ylim=(1.65,-.65),
               xlim=(0,.85) if j==0 else (0,1.25),xticks=[0,.4,.8] if j==0 else [0,.5,1.0],
               xlabel='MAE (percentage points)' if j==0 else 'MAE (g/dL)')
        grid(ax)

    panel(a,'C','Five-route Brier losses',12,370)
    panel(a,'D','Paired Brier contrasts',261,370)
    c=chart(fig,112,407,120,120);d=chart(fig,367,407,111,120)
    for i,key in enumerate(['B','M','S','MEAN','P']):
        interval(c,r['Brier'][key],i,('estimate','lower','upper'),
                 [GRAY,BLUE,ORANGE,GRAY,ORANGE][i],['s','o','s','^','D'][i],
                 tag=f'Figure6C:{key}',open_marker=key in ['S','MEAN'])
    c.set(yticks=range(5),yticklabels=['B · baseline','M · measured','S · fixed proxy',
          'MEAN · mean query','P · separate fit'],ylim=(4.65,-.65),xlim=(.041,.054),
          xticks=[.042,.048,.054],xlabel='Weighted Brier loss')
    grid(c)
    txt(a,12,390,'M / S / MEAN: one fitted M',7.5,weight='bold')
    for i,key in enumerate(['G_measured_gain','D_substitution_loss','H_remaining_gain']):
        interval(d,r['contrasts'][key],i,('estimate','lower','upper'),[BLUE,ORANGE,TEAL][i],
                 ['o','s','D'][i],tag=f'Figure6D:{key}',open_marker=i==1)
    d.set(yticks=range(3),yticklabels=['G · B − M\nMeasured gain','D · S − M\nSubstitution loss',
          'H · B − S\nRemaining gain'],ylim=(2.55,-.55),xlim=(-.0004,.0032),
          xticks=[0,.0015,.003],xlabel='Paired Brier difference')
    from matplotlib.ticker import FuncFormatter
    d.xaxis.set_major_formatter(FuncFormatter(lambda v,pos:f'{v:.4f}'.replace('-', '−')))
    c.xaxis.set_major_formatter(FuncFormatter(lambda v,pos:f'{v:.3f}'.replace('-', '−')))
    d.axvline(0,color=GRAY,ls=':',lw=.85);grid(d)
    txt(a,261,390,'G / D / H use paired losses',7.5,weight='bold')
    txt(a,12,566,'B / P separately fitted; P uses OOF training proxies.\nB–D: 95% survey-design t intervals; fixed fits.',7,GRAY)
    rr=r['retention'];lo,hi=rr['intervals'][0]
    txt(a,261,566,f"R = H/G: {rr['point']:.1%}\n95% Fieller: {lo:.1%} to {hi:.1%}".replace('-', '−'),7.5)
    txt(a,12,585,'Primary ≥40 years: positive remaining gain not established. Age 40–79 sensitivity: Table S30.',7,GRAY)
    save(fig,'Figure6_NHANES_Q3')

def fig2():
    rows=list(csv.DictReader((ROOT/'source/data/Fig3_BSEP_bridge_BSEP_intervals.csv').open(encoding='utf-8-sig')))
    fig,a=canvas(355)
    txt(a,12,19,'Public human BSEP vesicle assay',8,weight='bold')
    txt(a,12,34,r'Source-defined inhibition: $\mathrm{IC}_{50}$ <50 μM',7,GRAY)
    txt(a,12,56,'E201: 201 identities (37+/164−)   •   E172: nested subset (34+/138−)',7.5)
    line(a,12,70,481,70)
    panel(a,'A','Measurement prediction',12,88)
    txt(a,478,88,'Single-model 95% intervals',7,GRAY,ha='right')
    b=chart(fig,157,108,321,104)
    for i,(cohort,target) in enumerate([('E_BSEP','g_BSEPi_real'),('E_BSEP','g_BSEPi_perm'),
                                        ('E_no_anyA','g_BSEPi_real'),('E_no_anyA','g_BSEPi_perm')]):
        row=next(r for r in rows if r['cohort']==cohort and r['target']==target and r['metric']=='AP')
        interval(b,row,i,('point','low','high'),GRAY if i%2 else ORANGE,'x' if i%2 else 'o',
                 tag=f'Figure2A:{cohort}:{target}')
    b.set(yticks=range(4),yticklabels=['E201 · real','E201 · permuted','E172 (nested) · real','E172 (nested) · permuted'],
          ylim=(3.6,-.6),xlim=(0,.9),xticks=[0,.2,.4,.6,.8],xlabel='BSEP average precision (AP)')
    grid(b)
    # Local prevalence references span only the two matching cohort rows.
    prevalences=json.loads((ROOT/'source/data/display_references.json').read_text(encoding='utf-8'))['prevalence']
    for row,(ylo,yhi) in zip(prevalences,[(-.35,1.35),(1.65,3.35)]):
        name,num,den=row['cohort'],row['positive'],row['total']
        value=num/den
        b.plot([value,value],[ylo,yhi],color=GRAY,ls='--',lw=.8,zorder=1)
        b.text(value+.018,(ylo+yhi)/2,f'{value:.3f}',fontsize=7,color=GRAY,va='center')
        LEDGER.append(dict(tag=f'Added:Figure2A:prevalence:{name}',point=value,numerator=num,denominator=den,
                           unit='proportion',source='source/data/display_references.json#prevalence'))
    txt(a,157,247,'Dashed local marks: class-prevalence reference',7,GRAY)
    panel(a,'B','Matched supervision control',12,263)
    txt(a,478,263,'Paired 95% intervals',7,GRAY,ha='right')
    c=chart(fig,157,281,321,37)
    for i,cohort in enumerate(['E_BSEP','E_no_anyA']):
        row=next(r for r in rows if r['cohort']==cohort and r['target']=='real-perm' and r['metric']=='AP')
        interval(c,row,i,('point','low','high'),ORANGE,tag=f'Figure2B:{cohort}')
    c.set(yticks=[0,1],yticklabels=['E201','E172 (nested)'],ylim=(1.6,-.6),xlim=(0,.65),
          xticks=[0,.2,.4,.6],xlabel='Paired AP difference (real − permuted)')
    c.axvline(0,color=GRAY,ls=':',lw=.8);grid(c)
    # Definitions of absolute vs paired nominal 95% intervals stay in the legend.
    save(fig,'Fig2_BSEP_endpoint')

def fig3():
    data=saved()['Figure3'];fig,a=canvas(422)
    txt(a,12,17,'Drug-level DILI concern',8,weight='bold')
    txt(a,12,34,'125 V2 evaluation identities  ·  43 positive / 82 negative',7.5,GRAY)
    line(a,12,49,481,49)
    panel(a,'A','Downstream increment',12,67)
    b=chart(fig,153,89,325,59)
    for i,row in enumerate(data['paired']):
        for cov,lw in [('97_5',.8),('95',2.1)]:
            interval(b,row,i,('delta_AP',f'CI_{cov}_low',f'CI_{cov}_high'),ORANGE,
                     'D' if i==0 else 'p',lw=lw,tag=f'Figure3A:{i}:{cov}')
    b.set(yticks=[0,1],yticklabels=[r'$M_1 − M_0$',r'$R_{\mathrm{gate}} − M_0$'],ylim=(1.5,-.5),
          xlim=(-.06,.05),xlabel='Paired AP difference',xticks=[-.06,-.03,0,.03])
    b.axvline(0,color=GRAY,ls=':',lw=.8)
    txt(a,153,192,'95% thick; 97.5% thin paired intervals',7,GRAY)
    panel(a,'B','Average precision',153,218);panel(a,'C','ROC AUC',343,218)
    routes=['M0','M1','R_gate','Mperm','R_gate_perm','R_gate_struct']
    labs=[r'Structure ($M_0$)',r'Proxy concat. ($M_1$)',r'Residual ($R_{\mathrm{gate}}$)',
          'Concat. permutation','Residual permutation','Structural kernel']
    for j,(metric,x,w) in enumerate([('AP',153,144),('ROC_AUC',343,135)]):
        ax=chart(fig,x,239,w,129)
        for i,key in enumerate(routes):
            row=next(r for r in data['absolute'] if r['target']==key and r['metric']==metric)
            interval(ax,row,i,('point_estimate','low','high'),GRAY if i==0 or i>=3 else ORANGE,
                     ['s','D','p','x','x','h'][i],tag=f'Figure3:{key}:{metric}',open_marker=i==5)
        ax.set(yticks=range(6),yticklabels=labs if j==0 else [],ylim=(5.6,-.6),
               xlabel='AP' if j==0 else 'AUC',xlim=(.27,.59) if j==0 else (.39,.75),
               xticks=[.3,.4,.5] if j==0 else [.4,.5,.6,.7])
        ax.axvline(.344 if j==0 else .5,color=GRAY,ls=':',lw=.8);grid(ax)
    txt(a,12,413,'B–C: 95% intervals. Dotted references: prevalence 0.344 / chance 0.5.',7,GRAY)
    save(fig,'Fig3_DILI_N')

def fig5():
    data=saved()['Figure5'];fig,a=canvas(547)
    panel(a,'A','DILI seed sensitivity',12,17);panel(a,'B','Real ALT sensitivity',267,17)
    b=chart(fig,50,44,151,119);c=chart(fig,347,44,131,119)
    # Group is a descriptive selection category. Preserve the saved jitter geometry.
    for x,leaf,col,mark in [(0,1,ORANGE,'o'),(1,3,ORANGE,'D')]:
        rows=[r for r in data['A_all_seeds'] if int(r['M1_min_samples_leaf'])==leaf]
        xx=[x-.17+.34*i/(len(rows)-1) for i in range(len(rows))]
        vals=[float(r['DeltaAP']) for r in rows]
        b.scatter(xx,vals,color=col,marker=mark,s=17,edgecolors='white',lw=.4,zorder=3)
        LEDGER.extend(dict(tag=f"Figure5A:seed{r['root_seed']}",point=v) for r,v in zip(rows,vals))
    mean=data['A_stored_mean'];b.axhline(0,color=GRAY,lw=.7);b.axhline(mean,color=INK,ls='--',lw=.9)
    LEDGER.append(dict(tag='Figure5A:mean',point=mean))
    original=next(r for r in csv.DictReader((ROOT/'source/data/source_tables/N_primary_expected.csv').open(encoding='utf-8-sig'))
                  if r['comparison_id']=='original')
    primary=float(original['delta_AP'])
    b.axhline(primary,color=BLUE,ls='-.',lw=.8,zorder=1)
    LEDGER.append(dict(tag='Added:Figure5A:original_fixed_fit',point=primary,unit='AP difference',
                       source='source/data/source_tables/N_primary_expected.csv#comparison_id=original'))
    b.set(xticks=[0,1],xticklabels=['Leaf 1\n26 seeds','Leaf 3\n4 seeds'],
          ylabel=r'ΔAP ($M_1 − M_0$)',xlim=(-.45,1.45),ylim=(-.055,.112),yticks=[-.05,0,.05,.1])
    b.tick_params(axis='x',length=0)
    labs=['Zero response','No-Z','Forced-Z',r'$\lambda_Z=10$',r'$\lambda_Z=0.1$',r'$\lambda_Z=0.001$']
    for i,row in enumerate(data['B_saved_summary']):
        interval(c,row,i,('MAE_cmp_mean','MAE_cmp_min','MAE_cmp_max'),GRAY if i<2 else BLUE,
                 ['v','s','o','o','o','o'][i],tag=f"Figure5B:{row['route_model']}")
    c.set(yticks=range(6),yticklabels=labs,ylim=(5.65,-.65),xlim=(.201,.267),
          xlabel='Parent-equal MAE',xticks=[.21,.23,.25]);grid(c)
    txt(a,12,202,'30 paired refits; 18/30 positive; descriptive leaf groups.',7,GRAY)
    line(a,12,215,29,215,INK,.9,ls='--');txt(a,34,215,'Mean of 30 refits: +0.0125',7)
    line(a,12,228,29,228,BLUE,.8,ls='-.')
    txt(a,34,228,'Original fixed fit (Fig. 3A): −0.0093',7,BLUE)
    txt(a,267,213,'121 parents; mean of five repeat-pooled MAEs.\nWhiskers: min–max across repeats (not CI).',7,GRAY)
    line(a,12,240,481,240)
    panel(a,'C','Semisynthetic positive control',12,259)
    card(a,52,279,104,22,'PHH PC1 signal',BLUE,PALEBLUE,size=7)
    arrow(a,(158,290),(181,290))
    card(a,183,279,119,22,'Constructed target',GRAY,size=7)
    txt(a,316,290,'Defined within outer training',7,GRAY)
    d=chart(fig,61,343,417,156)
    for metric,label,col,mark,ls,ms in [('measured_gain','G: measured-input gain',BLUE,'o','-',4.3),
                                      ('substitution_loss','D: substitution loss',ORANGE,'s','--',6)]:
        rows=[r for r in data['C_saved_curves'] if r['metric']==metric]
        x=[float(r['gamma']) for r in rows];y=[float(r['mean']) for r in rows];se=[float(r['MCSE']) for r in rows]
        d.errorbar(x,y,yerr=se,color=col,marker=mark,ls=ls,ms=ms,lw=1.2,capsize=2.4,
                   mfc='none' if mark=='s' else col,mew=.9,label=label)
        LEDGER.extend(dict(tag=f"Figure5C:{metric}:gamma{r['gamma']}",point=v,MCSE=e) for r,v,e in zip(rows,y,se))
    d.set(xlabel='Constructed signal strength γ',ylabel='MAE gain or loss\n(synthetic target units)',
          xlim=(-.03,1.04),ylim=(-.035,.64),xticks=[0,.1,.25,.5,1],yticks=[0,.2,.4,.6])
    d.axhline(0,color=GRAY,lw=.7)
    d.legend(loc='lower left',bbox_to_anchor=(0,1.06),frameon=False,handlelength=2.1,ncol=2,borderaxespad=0)
    d.text(.985,.09,'γ = 1 · original noise\n≈0.6% of constructed measured-input\nMAE benefit retained',
           transform=d.transAxes,ha='right',va='bottom',fontsize=7)
    # In-panel R axis uses saved repeat-wise means, never a ratio of aggregate G/D.
    rows=[r for r in csv.DictReader((ROOT/'source/data/source_tables/S24_curves_and_MCSE.csv').open(encoding='utf-8-sig'))
          if r['metric']=='ratio_fixed' and float(r['sigma_multiplier'])==1]
    ri=chart(fig,101,370,136,64)
    ri.set_facecolor('white')
    valid=[r for r in rows if r['mean']!='']
    ri.errorbar([float(r['gamma']) for r in valid],[100*float(r['mean']) for r in valid],
                yerr=[100*float(r['MCSE']) for r in valid],color=TEAL,marker='o',ms=3.2,lw=.85,capsize=2)
    ri.set(xlim=(-.03,1.04),ylim=(-2.2,9.5),xticks=[0,.1,.25,.5,1],yticks=[-2,0,4,8],xlabel='')
    ri.text(.96,.06,'γ',transform=ri.transAxes,ha='right',va='bottom',fontsize=7)
    ri.set_xticklabels(['0\nNA','0.1','0.25','0.5','1'])
    ri.tick_params(length=2,pad=2)
    ri.axhline(0,color=GRAY,lw=.55)
    ri.set_title('Mean repeat-wise retention (%)\n± Monte Carlo SE',fontsize=7,fontweight='normal',pad=5)
    for row in rows:
        LEDGER.append(dict(tag=f"Added:Figure5C:R:gamma{row['gamma']}",
                           gamma=float(row['gamma']),point=100*float(row['mean']) if row['mean'] else None,
                           MCSE=100*float(row['MCSE']) if row['MCSE'] else None,
                           valid_repeats=int(row['valid_repeats']),unit='percent; MCSE in percentage points',
                           source='source/data/source_tables/S24_curves_and_MCSE.csv#ratio_fixed,sigma_multiplier=1'))
    txt(a,12,540,'C: means ± Monte Carlo SE; 20 repeats, original noise. G/D at identical γ; R at γ=0: NA (0/20).',7,GRAY)
    save(fig,'Fig5_sensitivity_positive_control')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--figures',nargs='+',type=int,default=[1,4,6,2,3,5]);args=p.parse_args()
    for n in args.figures: globals()[f'fig{n}']()
