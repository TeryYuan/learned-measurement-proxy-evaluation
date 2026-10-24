from b1_common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

BLUE='#27698B';GOLD='#B17824';GRAY='#6D737A';DARK='#27333D'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':DARK,'text.color':DARK,'axes.titleweight':'bold','pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white'})
def read(rel):return pd.read_csv(E/rel)
def table(df):
    def v(x):return f'{x:.6f}' if isinstance(x,(float,np.floating)) else str(x)
    return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(v(a) for a in row)+' |' for row in df.itertuples(index=False,name=None))
def savefig(fig,name,sources):
    for ext in ['png','pdf','svg']:fig.savefig(E/f'figures/{name}.{ext}',dpi=190,bbox_inches='tight')
    plt.close(fig)
    for key,d in sources.items():csv(d,E/f'figures/source_data/{name}_{key}.csv')
def figures(m,t):
    support=read('diagnostics/target_support_and_variance.csv').iloc[0];a=read('diagnostics/ALT_sampling_variance_and_shared_controls.csv');q=qweights(t.parent_id)
    fig,axs=plt.subplots(1,3,figsize=(14,4.1),layout='constrained');bins=np.linspace(t.y.min(),t.y.max(),25)
    axs[0].hist(t.y,bins=bins,weights=q,color=BLUE,alpha=.8);axs[0].axvline(0,color=DARK,lw=1);axs[0].set(xlabel='ALT response (log2 treated - control)',ylabel='Compound-weighted probability',title='A  Development target')
    axs[1].bar(['Between\nparents','Within\nparent'],[support.V_between,support.V_within],color=[BLUE,GOLD]);axs[1].set(ylabel='Weighted variance (log2 response squared)',title='B  Descriptive variance')
    axs[2].scatter(abs(a.y),a.estimated_SE,s=17,alpha=.6,color=BLUE,edgecolors='none');axs[2].set(xlabel='Absolute ALT response',ylabel='Estimated sampling SE',title='C  Individual ALT sampling')
    fig.suptitle('NAM_dev: 97 parents, 298 conditions | post-hoc diagnostics',fontsize=13)
    savefig(fig,'fig1_target_and_sampling',{'conditions':t,'sampling':a,'variance':pd.DataFrame([support])})
    rep=read('diagnostics/PHH_repeat_response_stability.csv').set_index('geneset');proxy=read('diagnostics/proxy_all50_metrics_and_training_mean.csv').set_index('geneset');can=read('diagnostics/Hallmark_cancellation_summary.csv').set_index('geneset');names=sorted(rep.index);ii=np.arange(50)
    fig,axs=plt.subplots(1,3,figsize=(15,17),sharey=True,layout='constrained');axs[0].scatter(rep.loc[names,'half_response_correlation'],ii,color=BLUE,s=28);axs[0].set_xlim(0,1);axs[0].set_yticks(ii,[n.replace('HALLMARK_','').replace('_',' ') for n in names],fontsize=9);axs[0].invert_yaxis();axs[0].set(xlabel='Half-response Pearson correlation',title='A  Source repeat stability')
    axs[1].scatter(proxy.loc[names,'MAE_minus_training_mean'],ii,color=BLUE,s=28);axs[1].axvline(0,color=GRAY,lw=1);axs[1].set(xlabel='Proxy MAE - outer-training-mean MAE',title='B  Proxy vs training mean')
    axs[2].scatter(can.loc[names,'median_cancellation_ratio'],ii,color=GOLD,s=28);axs[2].set_xlim(0,1);axs[2].set(xlabel='Median |mean response| / mean |response|',title='C  Within-set sign cancellation')
    for ax in axs:ax.grid(axis='y',alpha=.15)
    fig.suptitle('All 50 fixed Hallmark sets | 97 development parents | no set selection',fontsize=14)
    savefig(fig,'fig2_all50_source_and_proxy',{'repeat':rep.reset_index(),'proxy':proxy.reset_index(),'cancellation':can.reset_index()})
    props=read('diagnostics/actual_old_model_properties.csv');decomp=read('diagnostics/old_swap_squared_loss_decomposition_summary.csv');old=read('diagnostics/old_metrics_with_ablations.csv');oldp=old.loc[old.scope.eq('pooled_OOF')].set_index('route_model')
    fig,axs=plt.subplots(1,3,figsize=(14,4.4),layout='constrained')
    for k,(name,g) in enumerate(props.groupby('route_model')):axs[0].scatter(g.outer_fold+(k-1.5)*.04,g.dose_slope,marker=['o','s','^','D'][k],label=name,color=[BLUE,GOLD,GRAY,DARK][k],s=38)
    axs[0].set(xticks=[0,1,2],xlabel='Outer fold',ylabel='Old log10-dose coefficient',title='A  Additive-model slopes');axs[0].legend(fontsize=8)
    vals=decomp.iloc[0][['two_e_a','a_squared','actual_squared_loss_difference']].to_numpy(float);axs[1].bar(['2 e a','a squared','Total'],vals,color=[BLUE,GOLD,GRAY]);axs[1].axhline(0,color=DARK,lw=.8);axs[1].set(ylabel='Compound-weighted MSE contribution',title='B  Old swap loss identity');axs[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    order=['F_meas','F_swap','F_meanZ','F_struct','training_compound_mean'];axs[2].scatter(oldp.loc[order,'MAE_cmp'],np.arange(len(order)),color=BLUE);axs[2].set_yticks(np.arange(len(order)),['F_meas','F_swap','F_meanZ','F_struct','Training mean']);axs[2].invert_yaxis();axs[2].set(xlabel='Compound-equal MAE',title='C  Old model and input ablations');axs[2].set_xlim(.20,.24)
    savefig(fig,'fig3_old_model_diagnostics',{'slopes':props,'decomposition':decomp,'old_metrics':old})
    new=read('modeling/metrics.csv');newp=new.loc[new.scope.eq('pooled_OOF')].set_index('route_model');oldorder=['F_cond','F_struct','F_pred','F_meas','F_swap','F_meanZ','zero_response','training_compound_mean','training_weighted_median'];neworder=['D_cond','D_struct','D_meas_only','D_pred','D_meas','D_perm','D_swap','D_meanZ','zero_response','training_compound_mean','training_weighted_median']
    fig,axs=plt.subplots(2,2,figsize=(14,10),layout='constrained')
    for row,(data,order,title) in enumerate([(oldp,oldorder,'Phase4B0 preserved + stated ablations'),(newp,neworder,'Phase4B1 finite nested family')]):
        for col,(metric,label,lim) in enumerate([('MAE_cmp','Compound-equal MAE',(.19,.25)),('RMSE_cmp','Compound-equal RMSE',(.33,.38))]):
            ax=axs[row,col];ax.scatter(data.loc[order,metric],np.arange(len(order)),color=[GRAY if 'training_' in n or n=='zero_response' else BLUE for n in order],s=45);ax.set_yticks(np.arange(len(order)),[n.replace('training_compound_mean','Training mean').replace('training_weighted_median','Training median').replace('zero_response','Zero response') for n in order]);ax.invert_yaxis();ax.set(xlabel=label,xlim=lim,title=title);ax.grid(axis='x',alpha=.2)
    fig.suptitle('Same 97 parents and 298 OOF queries | lower error is better | no holdout results',fontsize=13)
    savefig(fig,'fig4_old_and_new_errors',{'old_metrics':old,'new_metrics':new})
    comp=read('modeling/comparisons.csv');o=read('modeling/OOF_predictions.csv');fig,axs=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    for ax,(name,title) in zip(axs[0], [('D_struct','A  Structure + condition'),('D_meas','B  Structure + measured PHH')]):
        d=o[o.route_model.eq(name)];ax.scatter(d.y,d.prediction,s=18,alpha=.65,color=BLUE,edgecolors='none');lo=min(d.y.min(),d.prediction.min());hi=max(d.y.max(),d.prediction.max());ax.plot([lo,hi],[lo,hi],color=GRAY,lw=1);ax.set(xlabel='Observed ALT log2 response',ylabel='Outer OOF prediction',title=title,xlim=(lo-.1,hi+.1),ylim=(lo-.1,hi+.1))
    for ax,(names,title) in zip(axs[1], [(['measured','predicted','substitution'],'C  Prespecified development differences'),(['meanZ_ablation','swap_vs_meanZ','permutation'],'D  Supporting controls')]):
        for k,name in enumerate(names):
            g=comp.loc[comp.comparison.eq(name)];pool=g.loc[g.scope.eq('pooled_OOF'),'delta_MAE'].iloc[0];ax.scatter(pool,k,s=55,color=BLUE,marker='D')
            for j in range(3):ax.scatter(g.loc[g.scope.eq(f'fold_{j}'),'delta_MAE'],k+(.14*(j-1)),s=25,color=GRAY,marker=['o','s','^'][j],alpha=.65)
        ax.set_yticks(range(len(names)),names);ax.invert_yaxis();ax.axvline(0,color=DARK,lw=1);ax.set(xlabel='Delta compound-equal MAE (signed as named)',title=title);ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    fig.suptitle('Development only | diamonds: pooled OOF; small marks: descriptive folds, not intervals',fontsize=12)
    savefig(fig,'fig5_predictions_and_differences',{'OOF':o,'comparisons':comp})
    sel=read('modeling/selected_lambdas_and_drop.csv');routes=['D_cond','D_struct','D_meas_only','D_pred','D_meas','D_perm'];fig,axs=plt.subplots(1,2,figsize=(11,5),layout='constrained');mapping={'0.001':0,'0.1':1,'10.0':2,'DROP':3}
    from matplotlib.colors import ListedColormap
    cmap=ListedColormap(['#214F68','#598FAC','#C5DDE8','#E9E9E9'])
    for ax,col,title in zip(axs,['lambdaX','lambdaZ'],['Structure block penalty','Biological block penalty']):
        mat=sel.pivot(index='route_model',columns='outer_fold',values=col).loc[routes];values=np.array([[mapping[str(v)] for v in row] for row in mat.to_numpy()]);ax.imshow(values,cmap=cmap,vmin=0,vmax=3,aspect='auto');ax.set_xticks([0,1,2],['Outer 0','Outer 1','Outer 2']);ax.set_yticks(range(len(routes)),routes);ax.set_title(title)
        for i in range(len(routes)):
            for j in range(3):ax.text(j,i,str(mat.iloc[i,j]),ha='center',va='center',color='white' if values[i,j]<2 else DARK,fontsize=11)
    fig.suptitle('Selections from training-only pooled inner MAE | lambdaC fixed at 0.01',fontsize=12)
    savefig(fig,'fig6_selected_penalties',{'selections':sel})

def run():
    done=json.loads((E/'modeling/completed.json').read_text(encoding='utf8'));assert done['ridge_fits']==531
    m,t,z,x=data();figures(m,t)
    metrics_new=read('modeling/metrics.csv');new=metrics_new.query("scope=='pooled_OOF'").set_index('route_model');metrics_old=read('diagnostics/old_metrics_with_ablations.csv');old=metrics_old.query("scope=='pooled_OOF'").set_index('route_model');diff=read('modeling/comparisons.csv');sel=read('modeling/selected_lambdas_and_drop.csv');a=read('diagnostics/ALT_sampling_variance_and_shared_controls.csv');rep=read('diagnostics/PHH_repeat_response_stability.csv');p=read('diagnostics/proxy_all50_metrics_and_training_mean.csv');ca=read('diagnostics/Hallmark_cancellation_per_parent.csv');sup=read('diagnostics/target_support_and_variance.csv').iloc[0];dist=read('diagnostics/target_distribution.csv');decomp=read('diagnostics/old_swap_squared_loss_decomposition_summary.csv').iloc[0];edge=read('diagnostics/shared_control_edges.csv');check=read('diagnostics/old_model_dose_and_swap_identity_checks.csv');block=read('modeling/training_block_support.csv');q=qweights(t.parent_id)
    summary=dict(mean_vhat_compound_equal=float(q@a.vhat.to_numpy()),median_condition_sampling_SE=float(a.estimated_SE.median()),min_condition_sampling_SE=float(a.estimated_SE.min()),max_condition_sampling_SE=float(a.estimated_SE.max()),shared_control_edges=len(edge),partial_overlap_edges=int((~edge.complete_same_control).sum()),PHH_correlation_min=rep.half_response_correlation.min(),PHH_correlation_median=rep.half_response_correlation.median(),PHH_correlation_max=rep.half_response_correlation.max(),proxy_R2_positive=int((p.R2>0).sum()),proxy_MAE_better_than_training_mean=int((p.MAE_minus_training_mean<0).sum()),cancellation_ratio_median=ca.cancellation_ratio.median(),legacy_dose_max_residual=check.max_dose_slope_residual.max(),legacy_swap_within_parent_max_range=check.within_parent_swap_range.max())
    dump(summary,E/'diagnostics/summary.json')
    distributions=[];oof=read('modeling/OOF_predictions.csv')
    for name,g in oof.groupby('route_model'):
        for scope,h in [('pooled_OOF',g)]+[(f'fold_{int(f)}',h) for f,h in g.groupby('outer_fold')]:
            for weighting,w in [('compound_equal',qweights(h.parent_id)),('condition_equal',np.ones(len(h))/len(h))]:
                for variable in ['y','prediction']:
                    v=h[variable].to_numpy();order=np.argsort(v,kind='stable');cum=np.cumsum(w[order]);row=dict(route_model=name,scope=scope,variable=variable,weighting=weighting,n_parents=h.parent_id.nunique(),n_conditions=len(h),mean=float(w@v),SD=float(np.sqrt(w@(v-w@v)**2)))
                    row.update({f'p{p}':float(v[order[min(np.searchsorted(cum,p/100-1e-14),len(v)-1)]]) for p in [0,5,25,50,75,95,100]});distributions.append(row)
    csv(distributions,E/'modeling/prediction_and_target_distribution.csv')
    executed=json.loads((E/'config/execution.json').read_text(encoding='utf8'));executed.update(status=done['status'],actual_run=done,final_configuration_copied_from_prefit_rules=True)
    dump(executed,E/'config/execution_completed.json')
    snapshot=pd.concat([read('metadata/readonly_input_snapshot.csv'),read('metadata/PHH_actual_read_inventory.csv')],ignore_index=True).drop_duplicates('path');pres=[]
    for r in snapshot.itertuples():
        s=Path(r.path).stat();pres.append(dict(path=r.path,unchanged=s.st_size==r.bytes and s.st_mtime_ns==r.mtime_ns,bytes_before=r.bytes,bytes_after=s.st_size,mtime_ns_before=r.mtime_ns,mtime_ns_after=s.st_mtime_ns))
    assert all(v['unchanged'] for v in pres);csv(pres,E/'metadata/legacy_preservation_check.csv')
    ts=ET.parse(E/'logs/final_tests.xml').getroot();suite=ts if ts.tag=='testsuite' else ts.find('testsuite');tests={k:int(suite.attrib[k]) for k in ['tests','failures','errors','skipped']};assert tests['failures']==tests['errors']==0;dump(tests,E/'logs/final_test_summary.json')
    env={k:__import__(mod).__version__ for k,mod in [('numpy','numpy'),('pandas','pandas'),('scipy','scipy'),('sklearn','sklearn'),('matplotlib','matplotlib'),('joblib','joblib')]};dump(dict(python=sys.version,executable=sys.executable,packages=env,cpu_workers=4,blas_threads=1),E/'logs/environment.json')
    counts=[]
    for f in range(3):
        s=t.loc[t.outer_fold.ne(f)];qq=t.loc[t.outer_fold.eq(f)];obj=joblib.load(E/f'models/outer{f}_D_meas.joblib');cats=obj.transformer.categories
        counts.append(dict(outer_fold=f,train_parents=s.parent_id.nunique(),query_parents=qq.parent_id.nunique(),train_conditions=len(s),query_conditions=len(qq),query_unknown_route_conditions=int((~qq.route.isin(cats)).sum()),known_training_routes=';'.join(cats),identity_overlap=len(set(s.parent_id)&set(qq.parent_id))))
    csv(counts,E/'modeling/outer_fold_support.csv')
    order=['zero_response','training_compound_mean','training_weighted_median','D_cond','D_struct','D_meas_only','D_pred','D_meas','D_perm','D_swap','D_meanZ'];oldorder=['zero_response','training_compound_mean','training_weighted_median','F_cond','F_struct','F_pred','F_meas','F_swap','F_meanZ']
    newtable=new.loc[order,['MAE_cmp','RMSE_cmp','R2_cmp']].reset_index();oldtable=old.loc[oldorder,['MAE_cmp','RMSE_cmp','R2_cmp']].reset_index();main=diff.loc[diff.scope.eq('pooled_OOF')&diff.comparison.isin(['measured','predicted','substitution']),['comparison','model_a','model_b','delta_MAE','delta_RMSE']]
    controls=diff.loc[diff.scope.eq('pooled_OOF')&~diff.comparison.isin(['measured','predicted','substitution']),['comparison','delta_MAE','delta_RMSE']]
    captions='''# 真实图稿说明\n\n全部图仅使用NAM_dev；PNG、PDF、SVG与CSV源数据位于figures。图中的折是开发描述，不是独立生物实验或置信区间。\n\n1. **fig1_target_and_sampling**：目标直方图按母体等权；方差精确分解为药物间和药物内，药物内包括剂量变化。每条件SE来自处理与对照个体logALT样本方差，不是模型误差上限。\n2. **fig2_all50_source_and_proxy**：50集合逐项源重复半份相关、代理相对于对应外训练均值的MAE差、集合内符号抵消比。每母体半份为1处理+1对照，供者身份未证实；不筛集合。\n3. **fig3_old_model_diagnostics**：旧模型原尺度log剂量斜率、2ea+a²分解、旧输入消融。斜率来自实际对象，不据此称程序错误。\n4. **fig4_old_and_new_errors**：旧结果与本轮新族分栏，横轴为明确标出的局部误差范围，点图不以柱长夸大差异；所有路线同97母体/298条件。旧图加入本轮明确标识的F_meanZ和中位数，不覆盖旧先导。\n5. **fig5_predictions_and_differences**：新OOF实际散点和规定差值；菱形为pooled，其他符号为三折描述，不是区间。差值必须按CSV中model_a减model_b解释。\n6. **fig6_selected_penalties**：每条路线的训练内选参，DROP为严格不含该块；D_cond原本无X/Z，D_struct原本无Z，D_meas_only原本无X，这些不应误算成可选块被丢弃的频率。\n'''
    (E/'reports/figure_captions_zh.md').write_text(captions,encoding='utf8')
    txt=f'''# Phase4B1：NAM信号诊断与有限剂量交互开发

交付时点：{now()}。**已完成开发集诊断、一次限定嵌套开发、保存对象核验及真实图稿。实测PHH输入在本轮仍未优于结构和简单参照，因此当前没有建立值得由proxy保留的实测增量。建议停止把当前ALT任务作为“预测替代实验”的主模型，不开展主动补测；如需完整限制性结论，可另行授权一次固定留出检验。** 本轮已经停止，不评价NAM_holdout。

本方法是在观察Phase4B0负向结果后提出的事后开发，不是原始预登记、原创算法证明或独立确认性验证。原alpha10结果、DILI Phase0—3和已用旧N均保留。二维身份、ALT终点、最高PHH浓度、50集合和原折未改变；没有删除大响应、难预测药物或集合。

## 范围、输入与真实完成

实际合格集合仍为121母体：**NAM_dev 97母体、298条件；NAM_holdout 24身份仅用于排除**。原失败NAM0014与NAM0067没有补抽。所有新目标分布、响应关联、诊断、变换、拟合、选参和预测均限开发集。全体原文件在读取后立即按开发身份子集；不计算留出结局统计。

在拟合前保存了 [方法修订](../config/method_amendment.md) 和 [执行配置](../config/execution.json)。实际读入Phase4B0报告、next_decisions、配置、分区/折、ALT动物/条件/对照、50集合和逐数组归一化基因矩阵、模型清单、12个Ridge及12个g对象、旧OOF。逐数组分析实际读取97个开发母体、388个源数组的保存矩阵；未读CEL作新转换。[输入快照](../metadata/readonly_input_snapshot.csv)、[逐母体表达读取清单](../metadata/PHH_actual_read_inventory.csv)、[开发成员](../metadata/dev_members.csv)、[共享折](../metadata/shared_development_folds.csv)、[开发共同条件](../derived/dev_queries.csv)可定位来源。

实际运行 **513个内层候选 + 18个外层选定重拟合 = 531次下游拟合**；补充27个U内部g，严格复用核对通过的12个旧g，所有lambda共享这些proxy。保存18个外层分块Ridge和27个新g，旧g只引用原对象。D_swap、D_meanZ额外拟合为0；全97生产拟合、留出预测、CEL重处理和新增下载均为0。完整状态见 [completed.json](../modeling/completed.json)。

## A. 目标与测量支持

目标仍是处理动物mean(log2 ALT)减匹配对照mean(log2 ALT)，单位为ALT log2响应；不是坏死、DILI或临床风险。当前开发集包括 **{int(sup.n_unique_treated)}只唯一合格处理动物、{int(sup.n_unique_controls)}只唯一对照**。同药多剂量及共享对照不增加独立药物数。

{table(dist)}

分位数为对应权重经验CDF的观测值，非插值；总体SD按相应权重计算。[完整目标统计](../diagnostics/target_distribution.csv)。母体等权方差：总 {sup.V_total:.9f} = 药物间 {sup.V_between:.9f} + 药物内 {sup.V_within:.9f}，代数残差 {sup.identity_residual:.3g}。药物内占约{100*sup.V_within/sup.V_total:.1f}%，含真实剂量变化，不能全部命名为测量噪声。[方差表](../diagnostics/target_support_and_variance.csv)、[逐条件完整贡献](../diagnostics/all_target_contributions_descriptive_not_exclusions.csv)保留所有尾部，没有winsorize或删除。

298个条件的有效处理动物为2—5只、对照均5只。按组内独立动物抽样等假设，vhat=s²T/nT+s²C/nC；其母体等权均值为 {summary['mean_vhat_compound_equal']:.6f}。条件SE中位数 {summary['median_condition_sampling_SE']:.6f}、范围 {summary['min_condition_sampling_SE']:.6f}—{summary['max_condition_sampling_SE']:.6f}，相较目标SD {dist.loc[dist.weighting.eq('compound_equal'),'SD'].iloc[0]:.6f}不可忽略，但尺度比较不能建立确定的信号/噪声分界。小样本方差不稳定，也没有将vhat从总方差相减以宣称可靠性或R²上限。

实际发现{len(edge)}条处理条件对共享对照边，全部为完整相同控制集合，部分重叠边{summary['partial_overlap_edges']}；没有跨母体共享边。每条完整共享边记录s²C/nC协方差项，原始动物ID、源文件和行号见 [ALT逐条件抽样表](../diagnostics/ALT_sampling_variance_and_shared_controls.csv)、[动物映射](../diagnostics/ALT_condition_animal_source_mapping.csv)、[对照边](../diagnostics/shared_control_edges.csv)、[复用模式](../diagnostics/control_reuse_patterns.csv)。不把条件误差视为相互独立药物证据。

PHH源重复按处理/控制各自源ID排序，2+2划为两个不重叠1+1半份；使用原冻结fRMA、原Entrez映射和全部Hallmark 2025.1.Hs H固定成员。50集合跨97母体的半份响应Pearson相关范围 **{rep.half_response_correlation.min():.3f}—{rep.half_response_correlation.max():.3f}，中位数{rep.half_response_correlation.median():.3f}**。这是当前源重复的响应稳定性，不是已证实跨供者验证；单个半份噪声较大。全部集合保留。

原g的OOF仍为{int((p.R2>0).sum())}/50集合R²>0；相对于各自外训练均值，{int((p.MAE_minus_training_mean<0).sum())}/50集合MAE更低。两种参照分母不同，不能混用，也不能据此筛选23或32个集合。[逐项代理表现及训练均值](../diagnostics/proxy_all50_metrics_and_training_mean.csv)。4850个母体×集合条目中，abs(mean r)/mean(abs r)中位数 {ca.cancellation_ratio.median():.3f}；符号抵消只说明均值摘要压缩，并不证明被抵消部分含ALT预测信息。[半份源映射](../diagnostics/PHH_repeat_half_source_mapping.csv)、[半份响应](../diagnostics/PHH_repeat_half_responses.csv)、[重复稳定性](../diagnostics/PHH_repeat_response_stability.csv)、[抵消逐药表](../diagnostics/Hallmark_cancellation_per_parent.csv)。

![目标与抽样诊断](../figures/fig1_target_and_sampling.png)

## B. 原加法模型的实际构造检查

实际模型和特征顺序与Phase4B0一致。同外折、同给药途径时，各母体共享原尺度log鼠剂量斜率，z仅改变药物截距。验证的最大斜率残差为 {check.max_dose_slope_residual.max():.3g}；每药swap−meas偏移的最大组内范围 {check.within_parent_swap_range.max():.3g}。这是已核实的表示约束，不是程序错误。[逐母体验证](../diagnostics/old_model_dose_and_swap_identity_checks.csv)、[真实对象属性](../diagnostics/actual_old_model_properties.csv)。

旧swap平方误差恒等式中，母体等权2ea={decomp.two_e_a:.9f}，a²={decomp.a_squared:.9f}，相加={decomp.actual_squared_loss_difference:.9f}；平均恒等式残差{decomp.identity_residual:.3g}。负平均损失差是当前误差与输入偏移共同作用的描述，不等于实测/代理等价，也不能据此宣称消除了实验噪声。[逐查询分解](../diagnostics/old_swap_squared_loss_decomposition.csv)。

新增旧对象F_meanZ只将外查询z替换成相应训练实测z加权均值，未重新拟合。新增训练中位数基线则按训练条件1/m权重的累计0.5下界观测值计算，不是每药均值的中位数。下表保留原七路线，另明确增加F_meanZ和中位数参照；原指标复算差小于1e-12。

{table(oldtable)}

旧F_swap MAE {old.loc['F_swap','MAE_cmp']:.6f}，F_meanZ {old.loc['F_meanZ','MAE_cmp']:.6f}。替换改善与仅消除个体z变化的数值很接近，不能独立支持“proxy保留了实测用途”。[旧OOF及本轮消融](../diagnostics/old_OOF_with_meanZ_and_median.csv)、[原指标核对](../diagnostics/old_metrics_reconciliation.csv)、[训练中位数](../diagnostics/training_median_baseline.csv)。

## C. 限定分块模型：实际方法和边界

六拟合路线严格为D_cond、D_struct、D_meas_only、D_pred、D_meas、D_perm；候选分别1/4/4/16/16/16。每当前训练集合使用q=1/(n母体×每药条件数)，总和1，预处理和损失同权。C含训练标准化PHH logμM、鼠log mg/kg及仅训练途径one-hot，未知途径全零。X的2048位只中心化，8描述符标准化；Z的50维按当前实测、内部proxy或训练置换输入各自标准化。X/Z均加入训练标准化鼠log剂量一阶交互，然后按训练加权block RMS能量缩放。常量置0，未按ALT筛列。没有μM/mg/kg换算、dose²、route×gene或PCA。

目标函数为Σq(y−ŷ)²+0.01||βC||²+λX||βX||²+λZ||βZ||²，截距不罚。两块各λ∈{{0.001,0.1,10,DROP}}，DROP精确不进入求解。它是归一化损失λ，不能与旧sklearn alpha10直接比较数值。实际用训练行空间的加权中心化Gram与Cholesky解，合成数据已对照独立primal解；所有候选复用变换/Gram和合法proxy，无数千维显式逆。

{table(pd.DataFrame(counts))}

外折沿用原分组；内层U/V逐项继承旧g成员。对U训练proxy，在U内部再三折产生；V用只在U训练的旧g_U预测。没有把旧S拼接proxy切成U训练输入。内层候选预测拼接后以母体等权MAE选择；1e-10内并列依少块、强Z、强X，不能以折MAE简单平均。各外层S再拟合所选模型，统一Q查询。D_perm仅在当前训练唯一母体间整体置换50维向量，种子2718+100×outer+inner（外重拟合inner=3），查询使用真实测量；没有改y、剂量或查询成员。

外层训练C中心Gram秩为3或4；结构交互4112列的中心Gram秩为126/128/128，实测Z交互100列均秩100。高维列数不等于有效自由度，也不据p/n直接诊断过拟合。全部训练上下文的常量、缩放能量、实际数值秩与容差见 [block支持](../modeling/training_block_support.csv)。所有X/Z非零block的训练加权能量均为1，没有零能量排药。

[proxy上下文与核对](../modeling/proxy_training_contexts.csv)、[下游边界](../modeling/downstream_training_contexts.csv)、[整向量置换映射](../modeling/permutation_maps.csv)、[全部内层候选指标](../modeling/inner_candidates.csv)、[候选逐条件预测](../modeling/inner_candidate_predictions.csv)、[实际拟合记录](../modeling/actual_fit_contexts.csv)均已保存。D_swap和D_meanZ确实复用同一D_meas及其训练变换，0次额外拟合。

## D. 新开发结果及合理参照

每条路线均覆盖97母体、298个相同外查询。MAE先药内平均再母体平均；RMSE为母体等权均方误差平方根；加权R²分母用同查询真实结局母体等权总体变异，不能把各折R²简单平均作pooled。条件等权结果只作补充，见CSV。

{table(newtable)}

三项规定主要差值（负MAE差仅表示当前开发误差数值更小）：

{table(main)}

完整支持性比较：

{table(controls)}

实测D_meas比旧F_meas的MAE下降{old.loc['F_meas','MAE_cmp']-new.loc['D_meas','MAE_cmp']:.6f}，但本轮同时改变结构缩放、交互和正则，不能把下降单独归因于其中一项。新实测相对结构差 **{main.loc[main.comparison.eq('measured'),'delta_MAE'].iloc[0]:+.9f}**；相对训练均值差 **{new.loc['D_meas','MAE_cmp']-new.loc['training_compound_mean','MAE_cmp']:+.6f}**。它仍未显示可用的实测增量；D_meas_only也高于D_cond及常数参照，置换路线没有更差。不能只比较两个弱模型并称实测有效。

D_pred−D_struct={main.loc[main.comparison.eq('predicted'),'delta_MAE'].iloc[0]:+.9f}，数值极小，而且对应路线仍弱于训练均值。D_swap−D_meas={main.loc[main.comparison.eq('substitution'),'delta_MAE'].iloc[0]:+.9f}，而D_meanZ−D_meas={new.loc['D_meanZ','MAE_cmp']-new.loc['D_meas','MAE_cmp']:+.9f}；swap与meanZ仅差{new.loc['D_swap','MAE_cmp']-new.loc['D_meanZ','MAE_cmp']:+.9f}。这些开发结果没有建立超出训练均值输入消融解释的替代用途，也不构成等效检验。未建立实测用途之前，不能宣称代理保留了该用途。

实测−结构逐折MAE差为0、−0.000001409、+0.000294154；没有一致改善。折数仅3，不做折t检验、确认性P值或把这些差异称统计显著。所有新拟合/消融路线加权R²为负；常数训练中位数pooled R²约0.001889，并不证明可部署预测效用。所有逐药与逐条件结果保留，不定义“获益亚组”。

所选参数如下（表中不可用块的DROP是路线定义；不能当成搜索后主动丢弃）：

{table(sel[['outer_fold','route_model','lambdaX','lambdaZ','selected_inner_MAE','n_candidates']])}

D_meas生物块一折DROP、两折λZ=10；D_pred两折生物块DROP，另折λZ=10。D_struct两折DROP结构，剩余一折λX=0.1。选择由训练内MAE决定，保留DROP分支并不保证外折表现更好。[选参及DROP](../modeling/selected_lambdas_and_drop.csv)、[频率](../modeling/lambda_DROP_frequency.csv)、[全OOF](../modeling/OOF_predictions.csv)、[全部指标](../modeling/metrics.csv)、[每药误差](../modeling/per_compound_errors.csv)、[每药差值](../modeling/per_parent_differences.csv)、[每条件差值](../modeling/per_condition_differences.csv)。

![旧结果与新开发分栏](../figures/fig4_old_and_new_errors.png)

## 解释限制与后续选择

必须区分四件事：当前目标分布接近零且有长尾、抽样与共享对照限制；原加法模型确实限制剂量斜率；新有限表示降低了旧模型误差；但这并没有建立实测输入在结构和常数参照上的增量。PHH源重复相关较高、部分集合proxy可预测，都不能替代ALT用途证据。均值摘要、单最高PHH浓度、人鼠体系/条件不等效、供者未知、有限药物和母体随机折仍限制推广。原因解释只是假设，没有做ALT相关选基因或新终点来验证它们。

历史PHH浓度和鼠剂量选择已使用活性、溶解度或先行毒性资料；结构+条件并非零成本结构输入。实际实测路线也不等于动物替代成功；目标只是单次24h ALT相对响应。

**推荐停止将当前ALT设定作为替代实测的主模型，不开展主动补测，不追加算法/集合/时点。** 若研究者希望形成完整的限制性结论，可另行确认一次固定24母体留出评价；这不是因为本轮建立了正向收益。最少待确认事项仅为：停止该主应用，或按 [未来固定留出协议草案](next_holdout_protocol_draft.md)做一次终局检验。草案明确未执行，参数冻结/全dev生产拟合/留出推断须有新授权。旧DILI N已经使用，绝不因新任务名刷新。

## 交付、测试及真实失败状态

- **完成**：A全部开发诊断、旧对象消融、六路线限定嵌套开发、两同对象查询消融、三基线、全部OOF/差值/选参/置换和模型上下文、六组真实图及源数据。
- **测试通过**：拟合前11项；含实际模型/选择/保留检查的最终{tests['tests']}项全部通过。验证包括独立primal/dual、权重、常量/未知途径、内部proxy边界、真实对象重构、训练内汇总选择及旧文件保留。[最终测试XML](../logs/final_tests.xml)、[测试摘要](../logs/final_test_summary.json)。
- **计算失败**：本轮诊断与限定训练无失败候选，无技术删药或填分。开始时两个文件定位尝试使用了不存在的路径，已定位实际共享技能说明和旧scripts/pilot.py；不涉及科学规则或数据修复。[定位记录](../logs/read_path_resolution_notes.json)。原两例CEL失败保持，无新补抽。
- **保留**：{len(pres)}份本轮记录的旧输入/对象/日志及开发表达文件大小和mtime未变；这是普通文件保留检查，不是哈希链。未重复数万项历史审查。[保留表](../metadata/legacy_preservation_check.csv)。
- **停止边界**：0留出预测/关联/结局诊断，0全97生产拟合，0主动补测/DrugMatrix/临床库/旧DILI评价或投稿。未上传推送或重建环境。

[18个外层模型清单](../modeling/outer_model_manifest.csv)、[27个新g清单](../modeling/new_proxy_model_manifest.csv)、[39个g新建/复用上下文](../modeling/proxy_training_contexts.csv)、[命令与资源日志](../logs/commands.jsonl)、[事件日志](../logs/events.jsonl)、[环境](../logs/environment.json)、[README](../README.md)。CPU最多4 worker，BLAS单线程，下游候选串行；缓存合法上下文复用。实际耗时及峰值RSS见命令日志，不能把写好代码当作计算完成。

全部六组图含PNG/PDF/SVG：目标/抽样、全部50集合、原模型诊断、旧新误差、预测与差值、选参。对应 [图注](figure_captions_zh.md) 和 [图稿索引](../figures/figure_index.csv) 提供逐图源表路径。
'''
    txt=txt.replace('[全部指标](../modeling/metrics.csv)','[全部指标](../modeling/metrics.csv)、[所有路线预测/目标分布](../modeling/prediction_and_target_distribution.csv)')
    txt=txt.replace('[执行配置](../config/execution.json)','[拟合前配置](../config/execution.json) 和 [实际完成配置](../config/execution_completed.json)')
    (E/'reports/phase4b1_NAM_signal_report_zh.md').write_text(txt,encoding='utf8')
    entry=txt.replace('](../','](../extensions/nam_phase4b1/').replace('](next_holdout_protocol_draft.md)','](../extensions/nam_phase4b1/reports/next_holdout_protocol_draft.md)').replace('](figure_captions_zh.md)','](../extensions/nam_phase4b1/reports/figure_captions_zh.md)')
    (ROOT/'reports/phase4b1_NAM_signal_report_zh.md').write_text(entry,encoding='utf8')
    protocol='''# 下一轮固定留出评价协议草案——本轮未执行

当前建议：停止当前ALT任务作为实测替代主应用；不做主动补测。如果研究者为完整限制性结论另行授权一次终局检验，才执行下述冻结设计。不能将本轮开发结论或本草案视为留出评价授权。

1. 保持最终121母体的原97/24分区、原ALT定义、精确条件、50固定集合和技术失败规则。24留出既不补抽也不用于选参、输入缩放、代理学习或终点诊断；旧DILI N仍为已使用数据。
2. 只沿用本轮六条有限路线及D_swap/D_meanZ同对象消融、三常数参照，不按开发最高分改名主方法。所有路线评价同一留出条件/母体。结构+条件仍包含先行实验选择成本。
3. 为全97训练选择参数时，使用原三折作为训练内U/V，严格按本轮网格、母体等权pooled MAE和并列规则。U训练proxy使用U内部交叉拟合，V由g_U产生；有效旧上下文可核对复用。各路线在97上只进行一次选定生产拟合。最终g训练只用97；D_pred训练用97内部合法crossfit，留出由g_97给出。当前文件没有这些全dev最终模型，也没有假定三个外折参数的多数表决可替代这一步。
4. D_perm全97的唯一训练向量置换须在下一轮拟合前固定单一context编号/种子并记录，不能挑坏的负控制；不改变留出测量。D_swap/D_meanZ复用全97最终D_meas。所有转换、block均值/能量、阈值（连续模型无阈值）不根据留出调整。
5. 下一轮先保存协议核对、模型/参数/上下文清单和不含结局的固定预测，再关联已保存留出ALT。技术失败明确报告，不按误差删药补值。
6. 主要仍为MAE(D_meas)-MAE(D_struct)、MAE(D_pred)-MAE(D_struct)、MAE(D_swap)-MAE(D_meas)，并完整报告常数/条件/仅实测/置换/meanZ控制、母体等权MAE/RMSE/R²。若需要不确定性，建议一次固定2000次配对母体bootstrap、PCG64(4109)，同抽样索引覆盖所有路线及差值，95%百分位区间明确仅涵盖固定模型在小留出样本上的抽样不确定性；这是待下一轮确认的统计口径，不是当前计算产物。不做确认性P值或等效/动物替代声明。
7. 一次评价后停止；无论方向都不据留出改算法再称其未触碰。新数据/终点/主动补测/DrugMatrix不属于此草案。若选择停止，无需先消耗留出以“证明阴性”。
'''
    (E/'reports/next_holdout_protocol_draft.md').write_text(protocol,encoding='utf8')
    readme='''# Phase4B1 实际交付

入口：reports/phase4b1_NAM_signal_report_zh.md。全部新产物与本轮代码位于本目录；旧Phase4B0输入和对象只读引用，无CEL复制/重处理、无留出预测、无全dev生产模型。

使用项目既有 `.venv/Scripts/python.exe`，不安装依赖。实际执行顺序：

```powershell
.\\.venv\\Scripts\\python.exe extensions\\nam_phase4b1\\scripts\\run.py initialize
.\\.venv\\Scripts\\python.exe extensions\\nam_phase4b1\\scripts\\run.py diagnostics
.\\.venv\\Scripts\\python.exe -m pytest extensions/nam_phase4b1/tests/test_rules.py -q
# pre_model_test_gate.json 实际通过后才运行；不可无依据覆盖
.\\.venv\\Scripts\\python.exe extensions\\nam_phase4b1\\scripts\\run.py train
.\\.venv\\Scripts\\python.exe -m pytest extensions/nam_phase4b1/tests -q
.\\.venv\\Scripts\\python.exe extensions\\nam_phase4b1\\scripts\\run.py report
```

本轮已完成，**不要重复研究拟合**。initialize和train检测已有完成状态会拒绝重复。report只读取已保存结果并重建图/报告，不触发模型训练。轻量验证可以只运行tests；没有CEL或R依赖。

加载新模型需将本目录src加入sys.path，类在banded.py/proxy_cache.py。旧模型定义来自nam_phase4b0/src/modeling.py；本轮只读底层sklearn对象，不调用会写旧日志的旧fit/predict包装。joblib仅限可信本地对象。

复现输入：metadata/dev_members.csv和原共享折；derived/dev_queries.csv、dev_PHH50.csv、dev_structure.npz、dev_ALT_individuals.csv；PHH源重复只引用旧dev基因矩阵。源身份/记录映射、全部候选及OOF、上下文和45个新模型都可追溯。要独立重现计算，请在另行授权的复制目录修改输出根，不覆盖本轮完成状态，不打包CEL或环境。

所有输出为事后开发证据，不是留出验证。Hallmark沿用2025.1.Hs H和原项目许可/来源，目标只是rat single24h ALT相对响应。新Ridge为已有分块正则思想的限定应用，不主张算法原创或动物替代。
'''
    readme=readme.replace('.\\.venv\\Scripts\\python.exe -m pytest extensions/nam_phase4b1/tests/test_rules.py -q','.\\.venv\\Scripts\\python.exe extensions\\nam_phase4b1\\scripts\\run.py pre_tests').replace('.\\.venv\\Scripts\\python.exe -m pytest extensions/nam_phase4b1/tests -q','.\\.venv\\Scripts\\python.exe extensions\\nam_phase4b1\\scripts\\run.py final_tests')
    (E/'README.md').write_text(readme,encoding='utf8')
    index=[]
    for path in sorted((E/'figures').glob('*.png')):
        index.append(dict(figure=path.stem,png=str(path.relative_to(E)),pdf=str(path.with_suffix('.pdf').relative_to(E)),svg=str(path.with_suffix('.svg').relative_to(E)),source_files=';'.join(str(a.relative_to(E)) for a in sorted((E/'figures/source_data').glob(path.stem+'_*.csv'))),scope='NAM_dev_only'))
    csv(index,E/'figures/figure_index.csv')
    dump([dict(path='nam_phase4b0/src/run_pilot.py',status='initial_nonexistent_path_resolved',actual='nam_phase4b0/scripts/pilot.py',scientific_effect='none')],E/'logs/read_path_resolution_notes.json')
    events=[json.loads(a) for a in (E/'logs/events.jsonl').read_text(encoding='utf8').splitlines()];fail=[a for a in events if 'failure' in a['event']];dump(dict(stage_failures=fail,failed_ridge_candidates=0,model_scientific_repairs=0,excluded_for_new_results=0),E/'logs/failure_summary.json')
    assets=[dict(path=str(p.relative_to(E)),bytes=p.stat().st_size) for p in sorted(E.rglob('*')) if p.is_file() and '__pycache__' not in str(p)];csv(assets,E/'reports/artifact_index.csv')
    dump(dict(**done,tests=tests,figures=len(index),model_objects=len(list((E/'models').glob('*.joblib'))),legacy_files_preserved=len(pres),output_MiB=sum(a['bytes'] for a in assets)/2**20,decision='recommend_stop_ALT_as_main_replacement_application; optional separately authorized fixed holdout for limiting conclusion'),E/'reports/actual_delivery_summary.json')
    log('report_completed',figures=len(index),tests=tests,model_objects=45);print('Report and six figure groups generated; actual artifacts indexed.')
