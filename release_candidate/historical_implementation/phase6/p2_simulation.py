"""new_random=yes synthetic noise; model_fit=yes PCA/transforms/dual Ridge; model_query=yes dev only; used_evaluation_data=no NAM holdout values."""
from common import *
import itertools,psutil,joblib
from scipy.linalg import cho_factor,cho_solve
LEG=ROOT/'extensions/nam_phase4b1';OLD=ROOT/'extensions/nam_phase4b0'
sys.path.insert(0,str(LEG/'src'))
from banded import BlockTransform,WeightedScale,candidates,select,tie_key
from b1_common import qweights
ROLE='conditional semisynthetic development-only simulation; not real ALT or new animal validation'
ROUTES=['D_cond','D_struct','D_meas','D_pred','D_perm'];EXTRA=['D_swap','D_meanZ']

def inputs():
    m=readcsv(LEG/'metadata/dev_members.csv');ids=set(m.parent_id)
    # Explicit columns exclude real ALT response and animal measurements.
    t=readcsv(LEG/'derived/dev_queries.csv',usecols=['condition_id','parent_id','route','log10_PHH_uM','log10_rat_mg_kg','outer_fold'])
    z=readcsv(LEG/'derived/dev_PHH50.csv').set_index('parent_id')
    a=np.load(LEG/'derived/dev_structure.npz',allow_pickle=False);x=pd.DataFrame(a['X'],index=a['parent_id'],columns=a['columns'])
    excluded=readcsv(LEG/'metadata/excluded_ids_only.csv',usecols=['parent_id'])
    assert len(ids)==97 and len(t)==298 and set(t.parent_id)==ids==set(z.index)==set(x.index) and not ids&set(excluded.parent_id)
    assert z.shape[1]==50 and x.shape[1]==2057 and t.condition_id.is_unique
    se=readcsv(LEG/'diagnostics/ALT_sampling_variance_and_shared_controls.csv',usecols=['condition_id','parent_id','estimated_SE']);assert set(se.parent_id)==ids and set(se.condition_id)==set(t.condition_id)
    sigma=float(se.estimated_SE.median());prior=json.loads((LEG/'diagnostics/summary.json').read_text('utf8'));assert abs(sigma-prior['median_condition_sampling_SE'])<1e-15
    dump(dict(actual_saved_table_median=sigma,historical_report_median=prior['median_condition_sampling_SE'],difference=sigma-prior['median_condition_sampling_SE'],reason='round-trip CSV parser preserves source decimals; earlier report used default parser',selected_source='actual dev SE table; no holdout values'),P/'simulation/sigma_source_check.json')
    return m,t,z,x,sigma

def proxies(m,z,x):
    ctx=readcsv(LEG/'modeling/proxy_training_contexts.csv');assert len(ctx)==39
    a=readcsv(OLD/'pilot/proxy_unique_parent_OOF.csv',usecols=['parent_id','geneset','proxy']);assert set(a.parent_id)==set(m.parent_id)
    outer=a.pivot(index='parent_id',columns='geneset',values='proxy').loc[:,z.columns]
    inner={f:readcsv(OLD/f'pilot/outer{f}_training_proxy_crossfit.csv').set_index('parent_id').loc[:,z.columns] for f in range(3)}
    innerinner={(f,j):readcsv(LEG/f'modeling/outer{f}_inner{j}_training_proxy_crossfit.csv').set_index('parent_id').loc[:,z.columns] for f in range(3) for j in range(3)}
    checked=[]
    for row in ctx.itertuples():
        tr=row.training_ids.split(';');va=row.query_ids.split(';');assert not set(tr)&set(va) and set(tr+va)<=set(m.parent_id)
        f=int(row.context[5]);saved=json.loads((LEG/f'modeling/outer{f}_context.json').read_text('utf8'))
        assert saved['geneset_order']==z.columns.tolist() and saved['proxy_input_columns']==x.columns.tolist()
        # Validate exact source prediction rows and their partition, without querying old models.
        if '_crossfit' in row.context:
            direct=readcsv(LEG/f'modeling/{row.context}_proxy_predictions.csv').set_index('parent_id');j=int(row.context.split('_inner')[1][0]);aggregate=innerinner[f,j]
        elif '_inner' in row.context:
            j=int(row.context.split('_inner')[1]);direct=inner[f].loc[va];aggregate=inner[f]
        else:direct=outer.loc[va];aggregate=outer
        assert set(direct.index)==set(va) and np.array_equal(direct.loc[va,z.columns],aggregate.loc[va,z.columns])
        assert np.isfinite(direct.loc[va,z.columns]).all().all()
        checked.append(dict(context=row.context,training_ids=row.training_ids,query_ids=row.query_ids,status='reused existing exact context; original B1 dev matrices and feature/target orders retained',new_proxy_fit=False,new_proxy_query=False,source_validation=row.checks))
    savecsv(checked,P/'simulation/proxy_cache_audit.csv',ROLE)
    return outer,inner,innerinner

class BatchDual:
    def __init__(self,blocks,q,context):
        self.q=q;self.sqrtq=np.sqrt(q);self.means={k:q@v for k,v in blocks.items()};self.c={k:v-self.means[k] for k,v in blocks.items()};self.grams={k:(v*self.sqrtq[:,None])@(v*self.sqrtq[:,None]).T for k,v in self.c.items()};self.factors={};self.context=context
    def solve(self,Y,lx,lz,query):
        pair=(lx,lz);pen={'C':.01}
        if lx!='DROP':pen['X']=float(lx)
        if lz!='DROP':pen['Z']=float(lz)
        if pair not in self.factors:
            G=np.eye(len(self.q))
            for k,lam in pen.items():G+=self.grams[k]/lam
            self.factors[pair]=cho_factor(G,lower=True)
            event('P2','factorization',context=self.context,lambdaX=lx,lambdaZ=lz,n=len(self.q))
        mu=self.q@Y;rhs=self.sqrtq[:,None]*(Y-mu);u=cho_solve(self.factors[pair],rhs);alpha=self.sqrtq[:,None]*u
        return self.query(alpha,mu,lx,lz,query),alpha,mu
    def query(self,alpha,mu,lx,lz,query):
        pen={'C':.01}
        if lx!='DROP':pen['X']=float(lx)
        if lz!='DROP':pen['Z']=float(lz)
        pred=np.tile(mu,(len(query['C']),1))
        for k,lam in pen.items():pred+=(query[k]-self.means[k])@self.c[k].T@alpha/lam
        assert np.isfinite(pred).all();return pred

def raw(t,x,z):return t[['log10_PHH_uM','log10_rat_mg_kg']].to_numpy(),x.loc[t.parent_id].to_numpy()[:,:-1],z.loc[t.parent_id].to_numpy(),t.route.to_numpy()

def run(resume=False,estimate=False):
    E=P/'simulation';E.mkdir(exist_ok=True);(E/'cache').mkdir(exist_ok=True);(E/'models').mkdir(exist_ok=True)
    cfg=spec()['P2'];started=time.time();budget=spec()['resources']['P2_wall_budget_hours']*3600
    m,t,z,x,sigma=inputs();ids=set(m.parent_id);sigmaeps=max(1e-12,.001*sigma)
    cells=[dict(cell_id=i,gamma=g,sigma_multiplier=s,sigma=s*sigma,repeat=r,noise_seed=f'SeedSequence([20261010,200,{r}])') for i,(g,s,r) in enumerate(itertools.product(cfg['gammas'],cfg['noise_scale_multipliers'],range(20)))]
    protocol=dict(locked_utc=now(),sigma0=sigma,sigma_source='nam_phase4b1 dev estimated_SE median',simulation_role=ROLE,generator='outer-training unique-parent standardized PHH50 full SVD PC1; sign and training SD fixed per outer',condition_mean='0.1 standardized PHH log10 uM + 0.1 standardized rat log10 mg/kg',noise='condition iid; common vector for each repeat over all gamma/sigma/outer',seed_rule='SeedSequence([20261010,200,repeat]) PCG64',order='outer 0,1,2 then gamma ascending; sigma0 then sigma0/2; repeat 0..19',inner_target='same locked outer-generator values for all inner contexts',ratio_denominator_threshold=sigmaeps,folds='original 3 outer, original inner3; not 97 LOPO',worker_policy='one process; BLAS1; no parallel heavy P3',budget_hours=24,cache_reuse='39 existing exact proxy contexts; no real ALT downstream model reused')
    if not (E/'protocol.json').exists():dump(protocol,E/'protocol.json');savecsv(cells,E/'seed_and_scenario_table.csv','protocol fixed before simulation outcomes')
    else:
        old=json.loads((E/'protocol.json').read_text('utf8'));assert old['sigma0']==sigma and old['seed_rule']==protocol['seed_rule']
    outproxy,innerproxy,innerinner=proxies(m,z,x)
    noise_path=E/'cache/noise.npz'
    if noise_path.exists():
        a=np.load(noise_path);assert np.array_equal(a['condition_ids'],t.condition_id);noise=a['noise'];event('P2','noise_reused',vectors=20)
    else:
        noise=np.column_stack([np.random.Generator(np.random.PCG64(np.random.SeedSequence([20261010,200,r]))).standard_normal(len(t)) for r in range(20)])
        np.savez_compressed(noise_path,noise=noise,condition_ids=t.condition_id.to_numpy(str));event('P2','new_random',distribution='standard_normal',vectors=20,values=int(noise.size),seed_rule=protocol['seed_rule'])
    contexts=readcsv(LEG/'metadata/fold_support_and_proxy_inner_members.csv');perm=readcsv(LEG/'modeling/permutation_maps.csv')
    allres=[];oracle=[];selections=[];gdiag=[];completed=[]
    for f in range(3):
        done=E/f'outer{f}_completed.json'
        if done.exists():
            allres.extend(readcsv(E/f'outer{f}_results.csv').to_dict('records'));oracle.extend(readcsv(E/f'outer{f}_oracle.csv').to_dict('records'));selections.extend(readcsv(E/f'outer{f}_selections.csv').to_dict('records'));gdiag.extend(readcsv(E/f'outer{f}_generator_diagnostics.csv').to_dict('records'));completed.append(f);continue
        begin=time.time();S=sorted(m.loc[m.outer_fold.ne(f),'parent_id']);Q=sorted(m.loc[m.outer_fold.eq(f),'parent_id']);it=np.flatnonzero(t.parent_id.isin(S));iq=np.flatnonzero(t.parent_id.isin(Q));ts=t.iloc[it];tq=t.iloc[iq]
        genpath=E/f'cache/generator_outer{f}.npz'
        if genpath.exists():
            a=np.load(genpath);signal=a['signal'];fc=a['f_C'];Y=a['targets'];event('P2','generator_reused',outer=f)
        else:
            Z=z.loc[S].to_numpy();zm=Z.mean(0);zs=Z.std(0);zs[zs==0]=1;Zs=(Z-zm)/zs
            _,sing,V=np.linalg.svd(Zs,full_matrices=False);v=V[0];v*=1 if v[np.argmax(abs(v))]>=0 else -1
            a=Zs@v;am=a.mean();sd=a.std();assert sd>0
            signals=pd.Series(((z-zm)/zs).to_numpy()@v,index=z.index);signals=(signals-am)/sd;signal=signals.loc[t.parent_id].to_numpy()
            w=qweights(ts.parent_id);cs=WeightedScale().fit(ts[['log10_PHH_uM','log10_rat_mg_kg']].to_numpy(),w);cstd=cs.transform(t[['log10_PHH_uM','log10_rat_mg_kg']].to_numpy());fc=.1*cstd[:,0]+.1*cstd[:,1]
            Y=np.column_stack([fc+c['gamma']*signal+c['sigma']*noise[:,c['repeat']] for c in cells]);assert np.isfinite(Y).all()
            np.savez_compressed(genpath,signal=signal,f_C=fc,targets=Y,z_mean=zm,z_scale=zs,loading=v,score_mean=am,score_scale=sd,condition_mean=cs.mean,condition_scale=cs.scale,training_ids=np.array(S),query_ids=np.array(Q),condition_ids=t.condition_id.to_numpy(str))
            event('P2','PCA_generator_fit',outer=f,training_parents=len(S),query_parents=len(Q),svd_solver='full',new_fit=1,new_random=False)
        diagnostics=dict(outer_fold=f,n_train_parents=len(S),n_query_parents=len(Q),fC_variance=float(np.var(fc[it])),signal_variance_conditions=float(np.var(signal[it])),fC_signal_correlation=float(np.corrcoef(fc[it],signal[it])[0,1]),role=ROLE)
        builders=[]
        for j in range(4):
            if j<3:
                row=contexts.loc[contexts.outer_fold.eq(f)&contexts.inner_fold.eq(j)].iloc[0];U=row.training_ids.split(';');V=row.validation_ids.split(';');predtrain=innerinner[f,j];predquery=innerproxy[f].loc[V];name=f'outer{f}_inner{j}'
            else:U=S;V=Q;predtrain=innerproxy[f];predquery=outproxy.loc[Q];name=f'outer{f}_refit'
            assert not set(U)&set(V) and set(U+V)<=ids
            if j<3:assert not set(U+V)&set(Q) and set(U+V)==set(S)
            iu=np.flatnonzero(t.parent_id.isin(U));iv=np.flatnonzero(t.parent_id.isin(V));tu=t.iloc[iu];tv=t.iloc[iv];w=qweights(tu.parent_id)
            pm=perm.loc[perm.context.eq(name)].set_index('parent_id').loc[U];assert set(pm.assigned_z_parent)==set(U)
            zp=pd.DataFrame(z.loc[pm.assigned_z_parent].to_numpy(),index=U,columns=z.columns)
            branches={}
            for branch,zt,zq in [('meas',z,z),('pred',predtrain,predquery),('perm',zp,z)]:
                tr=BlockTransform().fit(*raw(tu,x,zt),w);dual=BatchDual(tr.transform(*raw(tu,x,zt)),w,name+'_'+branch);query=tr.transform(*raw(tv,x,zq));branches[branch]=(tr,dual,query)
                event('P2','block_transform_fit',context=name,branch=branch,training_ids=U,query_ids=V,outer_query_ids=Q,columns_X=2056,columns_Z=50)
            extra={}
            if j==3:
                tr=branches['meas'][0];zmean=w@z.loc[tu.parent_id].to_numpy();meanframe=pd.DataFrame(np.tile(zmean,(len(V),1)),index=V,columns=z.columns)
                extra={'D_swap':tr.transform(*raw(tv,x,predquery)),'D_meanZ':tr.transform(*raw(tv,x,meanframe))}
            builders.append(dict(iu=iu,iv=iv,branches=branches,extra=extra,weights=qweights(tv.parent_id),U=U,V=V,name=name))
        context_build_seconds=time.time()-begin;results=[];ors=[];selected=[];candidate_rows=[]
        predictions={r:np.full((len(iq),len(cells)),np.nan) for r in ROUTES+EXTRA}
        # First complete outer/cell is the technical and timing context; then fixed remaining cells.
        batches=[[0],list(range(1,len(cells)))] if f==0 else [list(range(len(cells)))]
        for batch_index,batch in enumerate(batches):
            batch_start=time.time();inner_loss={};logical=0
            for route in ROUTES:
                branch='pred' if route=='D_pred' else 'perm' if route=='D_perm' else 'meas'
                for pair in candidates(route):
                    loss=np.zeros(len(batch))
                    for j,bld in enumerate(builders[:3]):
                        tr,dual,query=bld['branches'][branch];yp,_,_=dual.solve(Y[np.ix_(bld['iu'],batch)],*pair,query);delta=abs(yp-Y[np.ix_(bld['iv'],batch)])
                        loss+=len(bld['V'])/len(S)*(bld['weights']@delta);logical+=len(batch)
                    inner_loss[route,pair]=loss
                    candidate_rows.extend(dict(outer_fold=f,cell_id=cidx,route=route,lambdaX=pair[0],lambdaZ=pair[1],inner_MAE=float(loss[pos])) for pos,cidx in enumerate(batch))
            bld=builders[3]
            for route in ROUTES:
                branch='pred' if route=='D_pred' else 'perm' if route=='D_perm' else 'meas';tr,dual,query=bld['branches'][branch]
                picks=[select({pair:inner_loss[route,pair][c] for pair in candidates(route)}) for c in range(len(batch))]
                for pair in sorted(set(picks),key=tie_key):
                    loc=[i for i,p in enumerate(picks) if p==pair];cols=[batch[i] for i in loc];yy=Y[np.ix_(it,cols)];yp,alpha,mu=dual.solve(yy,*pair,query);logical+=len(cols)
                    predictions[route][:,cols]=yp
                    np.savez_compressed(E/f'models/outer{f}_batch{batch_index}_{route}_{str(pair).replace(" ","")}.npz',alpha=alpha,intercept=mu,cell_ids=np.array(cols),train_conditions=t.iloc[it].condition_id.to_numpy(str),lambdaX=str(pair[0]),lambdaZ=str(pair[1]))
                    for k,cidx in enumerate(cols):
                        c=cells[cidx];err=yp[:,k]-Y[iq,cidx];mae=float(bld['weights']@abs(err));mse=float(bld['weights']@(err*err))
                        results.append(dict(outer_fold=f,n_query_parents=len(Q),**c,route=route,lambdaX=pair[0],lambdaZ=pair[1],Z_block_retained=pair[1]!='DROP',MAE=mae,MSE=mse,RMSE=np.sqrt(mse),evidence_role=ROLE))
                        selected.append(dict(outer_fold=f,**c,route=route,lambdaX=pair[0],lambdaZ=pair[1],inner_MAE=inner_loss[route,pair][batch.index(cidx)],n_candidates=len(candidates(route))))
                    if route=='D_meas':
                        for ablation,qb in bld['extra'].items():
                            ap=dual.query(alpha,mu,*pair,qb);predictions[ablation][:,cols]=ap
                            for k,cidx in enumerate(cols):
                                err=ap[:,k]-Y[iq,cidx];mse=float(bld['weights']@(err*err));results.append(dict(outer_fold=f,n_query_parents=len(Q),**cells[cidx],route=ablation,lambdaX=pair[0],lambdaZ=pair[1],Z_block_retained=pair[1]!='DROP',MAE=float(bld['weights']@abs(err)),MSE=mse,RMSE=np.sqrt(mse),evidence_role=ROLE))
            for cidx in batch:
                c=cells[cidx];err=c['sigma']*noise[iq,c['repeat']];mse=float(bld['weights']@(err*err));ors.append(dict(outer_fold=f,**c,n_query_parents=len(Q),MAE=float(bld['weights']@abs(err)),MSE=mse,RMSE=np.sqrt(mse),role='oracle known noiseless generating mean; not deployable model'))
            event('P2','logical_solves_completed',outer=f,cell_ids=batch,count=logical,ablation_additional_fits=0,elapsed_seconds=time.time()-batch_start,rss_gib=psutil.Process().memory_info().rss/2**30)
            if f==0 and batch_index==0:
                dump(dict(status='completed_one_full_nested_outer_context',outer=0,cell=cells[0],logical_solves=logical,context_build_seconds=context_build_seconds,solve_seconds=time.time()-batch_start,rss_gib=psutil.Process().memory_info().rss/2**30,batch_factor_reuse='subsequent 199 RHS share current-context transforms/factors',direction_used_for_decision=False),E/'technical_timing_context.json')
                savecsv(results,E/'technical_context_results.csv',ROLE)
                if estimate:return
            savecsv(results,E/f'outer{f}_checkpoint_results.csv',ROLE)
        # Save input transformations; weights are in small dual-alpha files rather than repeated large primal vectors.
        joblib.dump(dict(outer=f,train_ids=S,query_ids=Q,context_blocks=[dict(name=b['name'],U=b['U'],V=b['V'],transforms={k:v[0] for k,v in b['branches'].items()},block_means={k:v[1].means for k,v in b['branches'].items()}) for b in builders],feature_order=x.columns.tolist(),geneset_order=z.columns.tolist()),E/f'models/outer{f}_transform_context.joblib',compress=3)
        assert all(np.isfinite(v).all() for v in predictions.values())
        np.savez_compressed(E/f'outer{f}_predictions.npz',**predictions,condition_ids=t.iloc[iq].condition_id.to_numpy(str),parent_ids=t.iloc[iq].parent_id.to_numpy(str),synthetic_target=Y[iq],cell_ids=np.arange(len(cells)))
        dump(dict(evidence_role=ROLE,outer=f,routes=ROUTES+EXTRA,dimensions='query condition by cell; target is outer-specific synthetic outcome',no_holdout=True),E/f'outer{f}_predictions.npz.metadata.json')
        savecsv(candidate_rows,E/f'outer{f}_inner_candidates.csv',ROLE)
        savecsv(results,E/f'outer{f}_results.csv',ROLE);savecsv(ors,E/f'outer{f}_oracle.csv',ROLE);savecsv(selected,E/f'outer{f}_selections.csv',ROLE);savecsv([diagnostics],E/f'outer{f}_generator_diagnostics.csv',ROLE)
        dump(dict(status='complete',outer=f,primary_rows=len(results),logical_solves=164*len(cells),factorizations=sum(len(d.factors) for b in builders for _,d,_ in b['branches'].values()),elapsed_seconds=time.time()-begin),done)
        allres.extend(results);oracle.extend(ors);selections.extend(selected);gdiag.append(diagnostics);completed.append(f);print(f'P2 outer {f} completed: {len(results)} rows',flush=True)
        if time.time()-started>=budget:break
    savecsv(allres,E/'sim_results.csv',ROLE);savecsv(oracle,E/'oracle_results.csv',ROLE);savecsv(selections,E/'selected_parameters.csv',ROLE);savecsv(gdiag,E/'generator_diagnostics.csv',ROLE)
    if len(completed)==3:
        d=pd.DataFrame(allres);assert len(d)==4200
        agg=[]
        for (g,sm,r,route),h in d.groupby(['gamma','sigma_multiplier','repeat','route']):
            w=h.n_query_parents.to_numpy()/97;assert abs(w.sum()-1)<1e-12
            agg.append(dict(gamma=g,sigma_multiplier=sm,repeat=r,route=route,MAE=float(w@h.MAE),MSE=float(w@h.MSE),RMSE=float(np.sqrt(w@h.MSE)),Z_retention=float(w@h.Z_block_retained.astype(float))))
        agg=pd.DataFrame(agg);savecsv(agg,E/'repeat_conditional_losses.csv',ROLE);rr=[]
        for (g,sm,r),h in agg.groupby(['gamma','sigma_multiplier','repeat']):
            h=h.set_index('route');gm=h.loc['D_struct','MAE']-h.loc['D_meas','MAE'];gp=h.loc['D_struct','MAE']-h.loc['D_pred','MAE'];sub=h.loc['D_swap','MAE']-h.loc['D_meas','MAE'];valid=gm>sigmaeps
            rr.append(dict(gamma=g,sigma_multiplier=sm,repeat=r,measured_increment=gm,predicted_increment=gp,substitution_loss=sub,relative_retrained_increment_ratio=gp/gm if valid else np.nan,fixed_model_retention=1-sub/gm if valid else np.nan,ratio_NA_reason='' if valid else 'nonpositive denominator' if gm<=0 else 'positive denominator at or below threshold',D_meas_Z_retention=h.loc['D_meas','Z_retention']))
        rr=pd.DataFrame(rr);savecsv(rr,E/'repeat_comparisons.csv',ROLE);summ=[]
        for (g,sm),h in rr.groupby(['gamma','sigma_multiplier']):
            for metric in ['measured_increment','predicted_increment','substitution_loss','relative_retrained_increment_ratio','fixed_model_retention','D_meas_Z_retention']:
                v=h[metric].dropna();summ.append(dict(gamma=g,sigma_multiplier=sm,sigma=sigma*sm,metric=metric,mean=v.mean(),MCSE=v.std(ddof=1)/np.sqrt(len(v)) if len(v)>1 else np.nan,valid_repeats=len(v),undefined_fraction=1-len(v)/20,MCSE_role='20-noise-repeat Monte Carlo SE; ratios conditional on valid denominator'))
        savecsv(summ,E/'sim_summary.csv',ROLE);plot(pd.DataFrame(summ),E)
    savecsv(pd.DataFrame(columns=['context','error']),E/'failures.csv',ROLE)
    finish('P2','complete' if len(completed)==3 else 'resource_checkpoint_incomplete',dict(completed_outer=completed,expected_outer=3,actual_primary_rows=len(allres),expected_primary_rows=4200,logical_solves=len(completed)*32800,new_PCA_fits=3,proxy_contexts_reused=39,holdout_value_reads=0,real_ALT_downstream_reused=False,elapsed_seconds=time.time()-started,actual_counts_source='events.jsonl; cache/reuse events separate from new operations'),f'# P2 实际执行\n\n完成外情境{completed}；主结果{len(allres)}/4200行。σ0={sigma:.17g}来自开发SE表。39类proxy上下文复用，模拟目标重新选参/求解。只有97开发母体数据；无24留出表达或ALT值读取。3外×3内替代97次LOPO已事前明示。\n\n每外训练PCA生成器锁定；内层不重造目标，20噪声重复跨γ/σ共用。γ=0不以必须DROP验收，γ>0不自动代表结构以外条件增量。比值只对正且超过{sigmaeps:.17g}的分母计算，不截断。Monte Carlo SE不是真实药物泛化区间。真实ALT的γ未知。\n')

def plot(s,E):
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axs=plt.subplots(2,3,figsize=(170/25.4,112/25.4));metrics=['D_meas_Z_retention','measured_increment','substitution_loss','relative_retrained_increment_ratio','fixed_model_retention']
    labels=['Measured Z retention','Measured increment (MAE)','Fixed substitution loss (MAE)','Retrained increment ratio','Fixed-model retention']
    for ax,metric,label in zip(axs.flat,metrics,labels):
        for sm,color in [(1.,'#21618c'),(.5,'#b55d24')]:
            h=s.loc[s.metric.eq(metric)&s.sigma_multiplier.eq(sm)].sort_values('gamma');ax.errorbar(h.gamma,h['mean'],yerr=h.MCSE,color=color,marker='o',markersize=3,label=f'σ = {sm:g} σ₀',capsize=2)
        ax.set(xlabel='Injected γ',ylabel=label);ax.axhline(0,color='.6',lw=.5);ax.spines[['top','right']].set_visible(False)
    ax=axs.flat[-1]
    for sm,color in [(1.,'#21618c'),(.5,'#b55d24')]:
        h=s.loc[s.metric.eq('fixed_model_retention')&s.sigma_multiplier.eq(sm)].sort_values('gamma');ax.plot(h.gamma,h.undefined_fraction,color=color,marker='o',markersize=3,label=f'σ = {sm:g} σ₀')
    ax.set(xlabel='Injected γ',ylabel='Undefined ratio fraction');axs.flat[0].legend(fontsize=7);fig.suptitle('Synthetic targets: conditional out-of-parent simulation',fontsize=9)
    for tag,ax in zip('ABCDEF',axs.flat):ax.text(-.12,1.04,tag,transform=ax.transAxes,fontweight='bold')
    fig.text(.5,.012,'Real ALT: γ unknown; historical DROP-Z identity is a separate observation, not γ = 0.',ha='center',fontsize=7)
    fig.tight_layout(rect=[0,.045,1,.96])
    for ext in ['pdf','svg','png']:fig.savefig(E/f'Fig_S3_simulation.{ext}',dpi=300)
    plt.close(fig)
