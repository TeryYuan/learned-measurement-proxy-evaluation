"""Read saved P6/P4B data only. No fit, prediction, PCA, or random-number generation."""
from common7 import *
from helpers7 import selection_crossing, gain_and_ratios, generator_variance, zero_event_upper_bound
import argparse
CHECKS=[]
def ck(v,label):
 assert bool(v),label
 CHECKS.append(label)
def rr(p):return read(track(p))
def run(out):
 E=OLD/'simulation'; L=ROOT/'extensions/nam_phase4b1'; out.mkdir(parents=True,exist_ok=True);track(OLD/'p2_simulation.py')
 d=rr(E/'sim_results.csv'); old_summary=rr(E/'sim_summary.csv'); sel=rr(E/'selected_parameters.csv'); cells=rr(E/'seed_and_scenario_table.csv')
 protocol=json.loads(track(E/'protocol.json').read_text('utf8'));eps=protocol['ratio_denominator_threshold'];sigma0=protocol['sigma0']
 t=rr(L/'derived/dev_queries.csv');m=rr(L/'metadata/dev_members.csv'); ck(len(m)==97 and len(t)==298,'dev member/condition shape')
 ck(not d.duplicated(['outer_fold','cell_id','route']).any() and len(d)==4200,'4200 unique saved losses')
 noise=np.load(track(E/'cache/noise.npz'));ck(np.array_equal(noise['condition_ids'],t.condition_id),'noise condition order')
 contexts=[];var=[];audit=[];losscheck=[];weights=[];blockaudit=[]
 for f in range(3):
  gen=np.load(track(E/f'cache/generator_outer{f}.npz'));pr=np.load(track(E/f'outer{f}_predictions.npz'))
  ck(np.array_equal(gen['condition_ids'],t.condition_id),f'generator {f} order')
  T=gen['training_ids'];Q=gen['query_ids'];ck(set(T).isdisjoint(Q) and set(T)|set(Q)==set(m.parent_id),f'generator {f} parents')
  ck(set(Q)==set(m.loc[m.outer_fold.eq(f),'parent_id']),f'generator {f} shared folds')
  contexts.append(dict(context='simulation_outer_refit',outer_fold=f,n_train=len(T),n_query=len(Q),inner_folds=3,noise_repeats=20,generator_scope='outer-training locked; no simulation_full97_selection',source_file=f'phase6/simulation/cache/generator_outer{f}.npz'))
  # Verify selections against saved coefficient objects, rather than only counting a report flag.
  for p in (E/'models').glob(f'outer{f}_batch*_D_meas_*.npz'):
   a=np.load(track(p));ids=a['cell_ids'];h=sel.loc[sel.outer_fold.eq(f)&sel.route.eq('D_meas')&sel.cell_id.isin(ids)]
   ck(len(h)==len(ids) and (h.lambdaZ.astype(str)==str(a['lambdaZ'])).all(),f'selection/coefficient identity {p.name}')
   blockaudit.extend(dict(outer_fold=f,cell_id=int(i),lambdaZ=str(a['lambdaZ']),coefficient_file=p.relative_to(ROOT).as_posix(),Z_input_columns=50,Z_interaction_columns=50,Z_selected=str(a['lambdaZ'])!='DROP',source_implementation='phase6/p2_simulation.py: BlockTransform and BatchDual; non-DROP Z uses actual measured 50+interaction block') for i in ids)
  for scope,ids in [('T',T),('Q',Q)]:
   idx=np.flatnonzero(t.parent_id.isin(ids));p=t.iloc[idx].parent_id;w=(1/p.map(p.value_counts())/len(ids)).to_numpy();ck(abs(w.sum()-1)<1e-12,f'weights {f} {scope}')
   weights.extend(dict(outer_fold=f,scope=scope,condition_id=t.iloc[i].condition_id,parent_id=t.iloc[i].parent_id,weight=float(wj)) for i,wj in zip(idx,w))
   fc=gen['f_C'][idx];s=gen['signal'][idx]
   for (g,sm),cellg in cells.groupby(['gamma','sigma_multiplier']):
    sg=float(cellg.sigma.iloc[0]);v=generator_variance(fc,s,w,g,sg)
    var.append(dict(outer_fold=f,generator_scope=scope,n_parents=len(ids),gamma=g,sigma_multiplier=sm,sigma=sg,weight_sum_squared=float(w@w),variance_definition='population mixture and expected finite fixed design kept distinct',**v))
    for cell in cellg.itertuples():
     e=noise['noise'][idx,cell.repeat];vv=generator_variance(fc,s,w,g,sg,e);y=gen['targets'][idx,cell.cell_id];vy=float(w@((y-w@y)**2));target_res=float(np.max(abs(y-(fc+g*s+sg*e))))
     ck(target_res<1e-12,f'target {f} {scope} {cell.cell_id}')
     ck(abs(vy-vv['finite_decomposition'])<1e-11,f'variance {f} {scope} {cell.cell_id}')
     cov=lambda a,b:float(w@((a-w@a)*(b-w@b)))
     audit.append(dict(outer_fold=f,generator_scope=scope,cell_id=cell.cell_id,gamma=g,sigma_multiplier=sm,repeat=cell.repeat,variance_noise_saved=cov(e,e),cov_C_noise=cov(fc,e),cov_s_noise=cov(s,e),finite_variance_from_saved_y=vy,finite_decomposition=vv['finite_decomposition'],residual=vy-vv['finite_decomposition'],target_reconstruction_max_abs=target_res))
  tq=t[t.parent_id.isin(Q)].set_index('condition_id').loc[pr['condition_ids']];ck(np.array_equal(tq.parent_id,pr['parent_ids']),f'prediction {f} parent order');w=(1/tq.parent_id.map(tq.parent_id.value_counts())/len(Q)).to_numpy()
  for route in d.route.unique():
   er=pr[route]-pr['synthetic_target'];mae=w@abs(er);mse=w@(er*er);h=d[d.outer_fold.eq(f)&d.route.eq(route)].set_index('cell_id').loc[pr['cell_ids']]
   ck(np.allclose(mae,h.MAE,atol=1e-12,rtol=1e-12),f'actual saved errors {f} {route} MAE');ck(np.allclose(mse,h.MSE,atol=1e-12,rtol=1e-12),f'actual saved errors {f} {route} MSE')
   losscheck.extend(dict(outer_fold=f,cell_id=int(c),route=route,MAE=float(mae[j]),MSE=float(mse[j]),n_query_parents=len(Q),absolute_MAE_residual=abs(float(mae[j]-h.MAE.iloc[j]))) for j,c in enumerate(pr['cell_ids']))
 save(blockaudit,out/'selection_object_audit.csv');save(weights,out/'generator_condition_weights.csv');save(var,out/'variance_decomposition.csv');save(audit,out/'finite_variance_identity.csv');save(losscheck,out/'saved_prediction_loss_replay.csv')
 b=pd.DataFrame(blockaudit);ck(len(b)==600 and not b.duplicated(['outer_fold','cell_id']).any(),'600 actual D_meas selected coefficient contexts')
 z=d.loc[d.route.eq('D_meas')].merge(b[['outer_fold','cell_id','Z_selected']],on=['outer_fold','cell_id'],validate='one_to_one');ck((z.Z_selected==z.Z_block_retained).all(),'Z presence versus saved loss flag')
 save(z,out/'D_meas_outer_selections.csv')
 loss=[]
 for keys,h in d.groupby(['gamma','sigma_multiplier','repeat','route']):
  w=h.n_query_parents.to_numpy()/len(m);ck(abs(w.sum()-1)<1e-12,f'outer weight {keys}')
  loss.append(dict(zip(['gamma','sigma_multiplier','repeat','route'],keys),MAE=float(w@h.MAE),MSE=float(w@h.MSE)))
 loss=pd.DataFrame(loss);save(loss,out/'repeat_losses.csv');comp=[]
 for (g,sm,r),h in loss.groupby(['gamma','sigma_multiplier','repeat']):
  x=h.set_index('route').MAE;zz=z[z.gamma.eq(g)&z.sigma_multiplier.eq(sm)&z.repeat.eq(r)];f=float((zz.n_query_parents/97)@zz.Z_selected.astype(float))
  comp.append(dict(gamma=g,sigma_multiplier=sm,repeat=r,weighted_selection=f,**gain_and_ratios(x.D_struct,x.D_meas,x.D_pred,x.D_swap,eps)))
 comp=pd.DataFrame(comp);save(comp,out/'repeat_gains_ratios_selection.csv');summ=[]
 metrics=['measured_gain','predicted_gain','substitution_loss','ratio_retrained','ratio_fixed','weighted_selection']
 for (g,sm),h in comp.groupby(['gamma','sigma_multiplier']):
  for metric in metrics:
   v=h[metric].dropna();summ.append(dict(gamma=g,sigma_multiplier=sm,metric=metric,mean=v.mean(),MCSE=v.std(ddof=1)/np.sqrt(len(v)) if len(v)>1 else np.nan,valid_repeats=len(v),undefined_fraction=1-len(v)/20,MC_unit='noise repeat after outer weighting'))
 summary=pd.DataFrame(summ)
 mapping={'measured_gain':'measured_increment','predicted_gain':'predicted_increment','substitution_loss':'substitution_loss','ratio_retrained':'relative_retrained_increment_ratio','ratio_fixed':'fixed_model_retention','weighted_selection':'D_meas_Z_retention'}
 for r in summary.itertuples():
  old=old_summary[old_summary.gamma.eq(r.gamma)&old_summary.sigma_multiplier.eq(r.sigma_multiplier)&old_summary.metric.eq(mapping[r.metric])].iloc[0]
  ck(np.allclose([r.mean,r.MCSE],[old['mean'],old.MCSE],atol=1e-12,rtol=1e-12,equal_nan=True),f'original summary {r.gamma} {r.sigma_multiplier} {r.metric}')
 save(summary,out/'curves_and_MCSE.csv');rows=[]
 for (g,sm),h in comp.groupby(['gamma','sigma_multiplier']):
  v=h[h.ratio_na_reason.isna()];rows.append(dict(gamma=g,sigma_multiplier=sm,valid_repeats=len(v),nonpositive=int(h.ratio_na_reason.eq('NONPOSITIVE_MEASURED_GAIN').sum()),small_positive=int(h.ratio_na_reason.eq('SMALL_POSITIVE_DENOMINATOR').sum()),mean_of_ratios_retrained=v.ratio_retrained.mean(),ratio_of_means_retrained=v.predicted_gain.mean()/v.measured_gain.mean() if len(v) else np.nan,mean_of_ratios_fixed=v.ratio_fixed.mean(),ratio_of_means_fixed=1-v.substitution_loss.mean()/v.measured_gain.mean() if len(v) else np.nan,empirical_positive_gain_frequency=float(h.measured_gain.gt(0).mean())))
 save(rows,out/'ratio_aggregation_and_NA.csv')
 save(z.groupby(['outer_fold','gamma','sigma_multiplier'],as_index=False).agg(selection_frequency=('Z_selected','mean'),n_repeats=('repeat','count')),out/'per_outer_selection_curves.csv')
 curve=summary[summary.metric.eq('weighted_selection')].copy();threshold=[]
 for sm,h in curve.groupby('sigma_multiplier'):
  h=h.sort_values('gamma')
  for p in [.5,.9]:
   a=selection_crossing(h.gamma,h['mean'],p);br=a.pop('bracket');a['upcrossing_brackets']=json.dumps(a['upcrossing_brackets']);gg=h.gamma.to_numpy();ff=h['mean'].to_numpy();a['all_crossings_with_direction']=json.dumps([dict(lower=float(gg[i-1]),upper=float(gg[i]),direction='up' if ff[i]>=p else 'down') for i in range(1,len(h)) if (ff[i-1]<p<=ff[i]) or (ff[i]<p<=ff[i-1])]);a.update(sigma_multiplier=sm,bracket_lower=br[0] if br else np.nan,bracket_upper=br[1] if br else np.nan,endpoint_MCSE_lower=float(h.loc[h.gamma.eq(br[0]),'MCSE'].iloc[0]) if br else np.nan,endpoint_MCSE_upper=float(h.loc[h.gamma.eq(br[1]),'MCSE'].iloc[0]) if br else np.nan,estimand='empirical_selection_crossing; not detection limit or confidence interval');threshold.append(a)
 save(threshold,out/'gamma_crossings.csv');thr=pd.DataFrame(threshold);wide=[]
 for sm,h in thr.groupby('sigma_multiplier'):
  a=h[h.level.eq(.5)].iloc[0];b=h[h.level.eq(.9)].iloc[0];wide.append(dict(sigma_multiplier=sm,gamma_star_50=a.interpolated_gamma,gamma_star_90=b.interpolated_gamma,status_50=a.status,status_90=b.status,label='empirical_selection_crossing; not a detection limit or real ALT confidence bound'))
 save(wide,out/'gamma_threshold_summary.csv');curve['q_sim']=1-curve['mean'];curve['empirical_q_ge_0_05']=curve.q_sim.ge(.05);curve['not_a_real_data_test']=True;save(curve,out/'scenario_nonretention_curve.csv')
 actual=rr(L/'figures/source_data/fig6_selected_penalties_selections.csv');actual=actual[actual.route_model.eq('D_meas')]
 for row in actual.itertuples():
  n=int(m.outer_fold.eq(row.outer_fold).sum());contexts.append(dict(context='actual_block_outer_refit',outer_fold=row.outer_fold,n_train=97-n,n_query=n,inner_folds=3,noise_repeats=1,lambdaX=row.lambdaX,lambdaZ=row.lambdaZ,actual_Z_retained=row.active_Z,source_file='extensions/nam_phase4b1/figures/source_data/fig6_selected_penalties_selections.csv'))
 final=rr(ROOT/'extensions/nam_phase4c/figures/source_data/fig4_final_parameters_selection.csv');f=final[final.route_model.eq('D_meas')].iloc[0]
 contexts.append(dict(context='actual_full97_selection',n_train=97,n_query=np.nan,inner_folds=3,noise_repeats=1,lambdaX=f.lambdaX,lambdaZ=f.lambdaZ,actual_Z_retained=f.actual_Z_retained,source_file='extensions/nam_phase4c/figures/source_data/fig4_final_parameters_selection.csv',query_values_read=False));save(contexts,out/'selection_context_comparison.csv')
 diag=json.loads(track(L/'diagnostics/summary.json').read_text('utf8'));v=rr(L/'diagnostics/target_support_and_variance.csv');old=rr(ROOT/'extensions/nam_phase4b0/pilot/compound_equal_weight_metrics.csv');old=old[old.scope.eq('pooled_OOF')].set_index('route_model');cmp=rr(L/'modeling/comparisons.csv');delta=cmp.loc[cmp.comparison.eq('measured')&cmp.scope.eq('pooled_OOF'),'delta_MAE'].iloc[0]
 anchors=[('compound_equal_target_SD',np.sqrt(v.V_total.iloc[0]),'4B1','diagnostics/target_support_and_variance.csv','97 observed dev; no model','scale; not regression residual SD'),('condition_SE_median',sigma0,'6P2','phase6/simulation/protocol.json','97 dev; historical exact sigma','simulation noise scale; not residual distribution'),('PHH_source_repeat_correlation_median',diag['PHH_correlation_median'],'4B1','diagnostics/summary.json','dev source repeats','repeat stability; not predictive utility'),('proxy_positive_R2_sets',diag['proxy_R2_positive'],'4B1','diagnostics/summary.json','dev 50-set OOF','23/50; not explained ALT variance'),('additive_meas_minus_struct_MAE',old.loc['F_meas','MAE_cmp']-old.loc['F_struct','MAE_cmp'],'4B0','extensions/nam_phase4b0/pilot/compound_equal_weight_metrics.csv','original alpha10 dev OOF','downstream increment in that context'),('block_meas_minus_struct_MAE',delta,'4B1','modeling/comparisons.csv','original three outer refits','context-specific dev increment')]
 save([dict(statistic=a,value=b,source_stage=c,source_file=d,population='NAM_dev97',training_context=e,what_it_calibrates=f,what_it_does_not_identify='real gamma, biological mechanism, holdout behavior') for a,b,c,d,e,f in anchors],out/'real_data_anchors.csv')
 save([dict(gamma_lower=np.nan,gamma_upper=np.nan,real_variance_share=np.nan,status='NOT_IDENTIFIABLE_FROM_EXISTING_DESIGN',reason=r) for r in ['outer64/65 selection differs from final97 event','real ALT is not established to follow fC+gamma*PC1; residual can absorb any coefficient change','sigma0 is sampling SE calibration, not known regression residual SD; heteroscedastic and shared-control real errors','20 correlated-context noise repeats cannot give precise 5% real-data inversion']],out/'real_data_position.csv')
 dump(dict(ideal_iid_zero_event_n=20,one_sided95_upper=zero_event_upper_bound(20),min_n_strictly_below_05=next(n for n in range(1,1000) if zero_event_upper_bound(n)<.05),not_applied_to_weighted_outer_frequency=True,rounded_example_fixed_retention=1-.566123/.569687,rounded_example_minus_0062=1-.566123/.569687-.0062,epsilon=eps,new_fit=0,new_model_query=0,new_random=0,NAM_holdout_values_read=0),out/'definitions_numeric_checks.json')
 save([dict(quantity=k,definition=v,source_script='phase6/p2_simulation.py + phase7/replay_phase7B.py',input='saved generator/predictions/losses/selection objects',order='condition -> equal parent -> nQ/97 within repeat -> 20 repeats') for k,v in [('G_meas','Lstruct-Lmeas'),('G_pred','Lstruct-Lpred'),('D_sub','Lswap-Lmeas'),('ratio_retrained','G_pred/G_meas when G_meas>epsilon'),('ratio_fixed','1-D_sub/G_meas when G_meas>epsilon'),('F_repeat','sum(nQ/97 * actual selected Z)'),('variance','population mixture and finite fixed-design formulas in protocols/04')]],out/'definition_source_map.csv')
 save([dict(check=x,passed=True) for x in CHECKS],out/'replay_checks.csv')
 dump(dict(status='complete',checks=len(CHECKS),rows=4200,actual_D_meas_contexts=600,new_fit=0,new_query=0,new_random=0,simulation_only=True),out/'completion.json')
 print('7B deterministic checks',len(CHECKS));return summary,threshold

def plot(s,out):
 import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
 fig,axes=plt.subplots(2,3,figsize=(7.1,4.8));items=[('weighted_selection','Measured Z selection frequency'),('measured_gain','Measured MAE gain'),('substitution_loss','Fixed-input substitution loss'),('ratio_retrained','Retrained gain ratio'),('ratio_fixed','Fixed-model retention ratio')]
 for ax,(key,label) in zip(axes.flat,items):
  for sm,col in [(1.,'#23669a'),(.5,'#b05b24')]:
   h=s[s.metric.eq(key)&s.sigma_multiplier.eq(sm)].sort_values('gamma');ax.errorbar(h.gamma,h['mean'],h.MCSE,color=col,marker='o',markersize=3,capsize=2,label=f'noise {sm:g} x sigma0')
  ax.set(xlabel='Injected gamma',ylabel=label);ax.axhline(0,c='.6',lw=.5)
 axes[0,0].axhline(.5,c='.5',ls=':',lw=.6);axes[0,0].axhline(.9,c='.5',ls=':',lw=.6);axes[0,0].legend(fontsize=6)
 for sm,col in [(1.,'#23669a'),(.5,'#b05b24')]:
  h=s[s.metric.eq('ratio_fixed')&s.sigma_multiplier.eq(sm)].sort_values('gamma');axes[1,2].plot(h.gamma,h.undefined_fraction,color=col,marker='o',markersize=3)
 axes[1,2].set(xlabel='Injected gamma',ylabel='Undefined-ratio fraction')
 for tag,ax in zip('ABCDEF',axes.flat):ax.spines[['top','right']].set_visible(False);ax.text(-.16,1.04,tag,transform=ax.transAxes,fontweight='bold')
 fig.suptitle('Post-hoc simulation replay: outer-training selection contexts',fontsize=10)
 fig.text(.5,.012,'20 shared noise repeats; bars: Monte Carlo SE. Real ALT gamma is not identified.',ha='center',fontsize=8);fig.tight_layout(rect=[0,.055,1,.96])
 for ext in ['pdf','svg','png']:
  fig.savefig(out/f'Fig_threshold.{ext}',dpi=400);meta(out/f'Fig_threshold.{ext}',source='curves_and_MCSE.csv; saved results only',uncertainty='Monte Carlo SE, not real-data confidence limits',raster_dpi=400 if ext=='png' else None)
 plt.close(fig)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=B/'replay' if args.verify else B
 s,th=run(out)
 if not args.verify:plot(s,out);event('7B','complete',checks=len(CHECKS),new_fit=0,new_query=0,new_random=0)
 else:
  for p in out.glob('*.csv'):
   pd.testing.assert_frame_equal(read(p),read(B/p.name),check_exact=False,atol=1e-12,rtol=1e-12)
  print('All deterministic output tables agree with original Phase7B tables')
