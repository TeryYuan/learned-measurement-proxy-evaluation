from b1_common import *
from banded import *
from proxy_cache import ProxyCache

def raw(t,x,z):
    return (t[['log10_PHH_uM','log10_rat_mg_kg']].to_numpy(),x.loc[t.parent_id].to_numpy()[:,:-1],z.loc[t.parent_id].to_numpy(),t.route.to_numpy())
def output(t,p,route):
    a=t[['condition_id','parent_id','source_name','y','dose_mg_kg','route','log10_PHH_uM','outer_fold']].copy();a['route_model']=route;a['prediction']=p;a['error']=p-a.y;a['absolute_error']=abs(a.error);a['squared_error']=a.error**2;return a
def run():
    if (E/'modeling/completed.json').exists():raise RuntimeError('Finite family already completed; refusing repeat training')
    if not (E/'diagnostics/completed.json').exists():raise RuntimeError('Diagnostics must precede modeling')
    test=json.loads((E/'logs/pre_model_test_gate.json').read_text(encoding='utf8'));assert test['passed']
    m,t,z,x=data();allowed=set(m.parent_id);pc=ProxyCache(m,z,x);inner_original=pd.read_csv(E/'metadata/fold_support_and_proxy_inner_members.csv')
    all_predictions=[];inner_predictions=[];inner_scores=[];selected=[];perms=[];supports=[];manifest=[];fitrows=[];boundaries=[]
    for f in range(3):
        S=sorted(m.loc[m.outer_fold.ne(f),'parent_id']);Q=sorted(m.loc[m.outer_fold.eq(f),'parent_id']);ts=t.loc[t.parent_id.isin(S)].copy();tq=t.loc[t.parent_id.isin(Q)].copy();guard(S+Q,allowed)
        zscf,zq=pc.outer_pair(f,S,Q);pooled={(r,p):[] for r in ROUTES for p in candidates(r)}
        for j in range(3):
            ctx=inner_original.loc[inner_original.outer_fold.eq(f)&inner_original.inner_fold.eq(j)].iloc[0]
            U=ctx.training_ids.split(';');V=ctx.validation_ids.split(';');assert not set(U)&set(V) and set(U)|set(V)==set(S)
            tu=t.loc[t.parent_id.isin(U)].copy();tv=t.loc[t.parent_id.isin(V)].copy();zu,zv=pc.inner_pair(f,j,U,V)
            zp,maps=permute_vectors(U,z,f,j);perms.extend(maps);w=qweights(tu.parent_id);contexts={}
            for branch,zt,zvinput in [('meas',z,z),('pred',zu,zv),('perm',zp,z)]:
                rr=raw(tu,x,zt);transform=BlockTransform().fit(*rr,w);b=transform.transform(*rr);dual=DualContext(b,tu.y.to_numpy(),w);qc=dual.query_cache(transform.transform(*raw(tv,x,zvinput)));contexts[branch]=(dual,qc)
                for row in transform.support(b,w):supports.append(dict(context=f'outer{f}_inner{j}',branch=branch,**row))
            boundaries.append(dict(context=f'outer{f}_inner{j}',role='inner_candidate_training',training_ids=';'.join(U),query_ids=';'.join(V),outer_query_ids=';'.join(Q),training_conditions=len(tu),query_conditions=len(tv),identity_overlap=0,weight_sum=w.sum(),min_parent_weight=pd.Series(w).groupby(tu.parent_id.to_numpy()).sum().min(),max_parent_weight=pd.Series(w).groupby(tu.parent_id.to_numpy()).sum().max()))
            for route in ROUTES:
                branch='pred' if route=='D_pred' else 'perm' if route=='D_perm' else 'meas';dual,qcache=contexts[branch]
                for lx,lz in candidates(route):
                    start=time.time();pred,_=dual.solve(lx,lz,qcache)
                    rec=dict(context=f'outer{f}_inner{j}_{route}_{keypair(lx,lz)}',outer_fold=f,inner_fold=j,route_model=route,lambdaX=lx,lambdaZ=lz,role='inner_candidate',status='completed',training_parents=len(U),query_parents=len(V),elapsed_seconds=time.time()-start)
                    fitrows.append(rec);log('banded_ridge_fit',**rec)
                    frame=output(tv,pred,route);frame['inner_fold']=j;frame['selection_outer_fold']=f;frame['lambdaX']=str(lx);frame['lambdaZ']=str(lz);inner_predictions.append(frame);pooled[(route,(lx,lz))].append(frame)
            pc.save();print(f'outer {f} inner {j}: legal proxies and 57 candidates completed',flush=True)
        chosen={}
        for route in ROUTES:
            scores={}
            for pair in candidates(route):
                d=pd.concat(pooled[(route,pair)],ignore_index=True);assert set(d.parent_id)==set(S) and d.condition_id.is_unique
                met=metrics(d.parent_id,d.y,d.prediction);scores[pair]=met['MAE_cmp'];inner_scores.append(dict(outer_fold=f,route_model=route,lambdaX=pair[0],lambdaZ=pair[1],candidate_id=keypair(*pair),**met))
            pick=select(scores);chosen[route]=pick;selected.append(dict(outer_fold=f,route_model=route,lambdaC=.01,lambdaX=pick[0],lambdaZ=pick[1],active_X=pick[0]!='DROP',active_Z=pick[1]!='DROP',selected_inner_MAE=scores[pick],n_candidates=len(scores),best_inner_MAE=min(scores.values()),n_tied=sum(v<=min(scores.values())+1e-10 for v in scores.values())))
        zp,maps=permute_vectors(S,z,f,3);perms.extend(maps);w=qweights(ts.parent_id);contexts={};meanz=w@z.loc[ts.parent_id].to_numpy()
        for branch,zt,zvinput in [('meas',z,z),('pred',zscf,zq),('perm',zp,z)]:
            rr=raw(ts,x,zt);transform=BlockTransform().fit(*rr,w);b=transform.transform(*rr);dual=DualContext(b,ts.y.to_numpy(),w);qc=dual.query_cache(transform.transform(*raw(tq,x,zvinput)));contexts[branch]=(transform,dual,qc,zvinput)
            for row in transform.support(b,w):supports.append(dict(context=f'outer{f}_refit',branch=branch,**row))
        boundaries.append(dict(context=f'outer{f}_refit',role='outer_selected_training',training_ids=';'.join(S),query_ids=';'.join(Q),outer_query_ids=';'.join(Q),training_conditions=len(ts),query_conditions=len(tq),identity_overlap=0,weight_sum=w.sum(),min_parent_weight=pd.Series(w).groupby(ts.parent_id.to_numpy()).sum().min(),max_parent_weight=pd.Series(w).groupby(ts.parent_id.to_numpy()).sum().max()))
        for route in ROUTES:
            branch='pred' if route=='D_pred' else 'perm' if route=='D_perm' else 'meas';transform,dual,qc,zvinput=contexts[branch];lx,lz=chosen[route]
            start=time.time();pred,coef=dual.solve(lx,lz,qc,True);key=f'outer{f}_{route}'
            model=BandedModel(transform,dual,coef,S,allowed,lx,lz,key,meanz,x.columns[:-1].tolist(),z.columns.tolist())
            direct=model.predict(tq.parent_id.tolist(),*raw(tq,x,zvinput));assert np.allclose(direct,pred,atol=1e-10,rtol=0)
            path=E/f'models/{key}.joblib';joblib.dump(model,path,compress=3);manifest.append(dict(context=key,type='selected_outer_banded_Ridge',path=str(path.relative_to(E)),training_ids=';'.join(S),query_ids=';'.join(Q),lambdaX=lx,lambdaZ=lz,production_model=False))
            rec=dict(context=key,outer_fold=f,inner_fold=None,route_model=route,lambdaX=lx,lambdaZ=lz,role='outer_selected',status='completed',training_parents=len(S),query_parents=len(Q),elapsed_seconds=time.time()-start);fitrows.append(rec);log('banded_ridge_fit',**rec);all_predictions.append(output(tq,pred,route))
            if route=='D_meas':
                for ablation,zinput in [('D_swap',zq),('D_meanZ',pd.DataFrame(np.tile(meanz,(len(Q),1)),index=Q,columns=z.columns))]:
                    ap=model.predict(tq.parent_id.tolist(),*raw(tq,x,zinput));all_predictions.append(output(tq,ap,ablation));log('same_object_query_ablation',context=key,ablation=ablation,model_path=str(path),query_ids=Q,additional_fits=0)
        for name,val in [('zero_response',0),('training_compound_mean',w@ts.y.to_numpy()),('training_weighted_median',median(ts.y.to_numpy(),w))]:all_predictions.append(output(tq,np.full(len(tq),val),name))
        csv(selected,E/'modeling/selected_lambdas_and_drop.csv');csv(fitrows,E/'modeling/actual_fit_contexts.csv');csv(manifest,E/'modeling/outer_model_manifest.csv');csv(inner_scores,E/'modeling/inner_candidates.csv');csv(perms,E/'modeling/permutation_maps.csv');csv(supports,E/'modeling/training_block_support.csv');csv(boundaries,E/'modeling/downstream_training_contexts.csv');csv(pd.concat(all_predictions,ignore_index=True),E/'modeling/OOF_predictions.csv');csv(pd.concat(inner_predictions,ignore_index=True),E/'modeling/inner_candidate_predictions.csv');pc.save()
        dump(dict(outer_fold=f,S=S,Q=Q,geneset_order=z.columns.tolist(),proxy_input_columns=x.columns.tolist(),selected={r:list(p) for r,p in chosen.items()},training_conditions=ts.condition_id.tolist(),query_conditions=tq.condition_id.tolist(),same_D_meas_object_for_swap_and_meanZ=True),E/f'modeling/outer{f}_context.json')
        print(f'outer {f} selected models and same-object ablations completed',flush=True)
    oo=pd.concat(all_predictions,ignore_index=True);mets=[]
    for route,g in oo.groupby('route_model'):
        assert len(g)==len(t) and set(g.parent_id)==allowed and g.condition_id.is_unique
        mets.append(dict(route_model=route,scope='pooled_OOF',**metrics(g.parent_id,g.y,g.prediction),condition_MAE=abs(g.error).mean(),condition_RMSE=np.sqrt(np.mean(g.error**2))))
        for f,h in g.groupby('outer_fold'):mets.append(dict(route_model=route,scope=f'fold_{int(f)}',**metrics(h.parent_id,h.y,h.prediction),condition_MAE=abs(h.error).mean(),condition_RMSE=np.sqrt(np.mean(h.error**2))))
    met=pd.DataFrame(mets);csv(met,E/'modeling/metrics.csv');perdrug=oo.groupby(['route_model','parent_id','source_name','outer_fold'],as_index=False).agg(n_conditions=('condition_id','size'),mean_absolute_error=('absolute_error','mean'),mean_squared_error=('squared_error','mean'),mean_prediction=('prediction','mean'),mean_y=('y','mean'));csv(perdrug,E/'modeling/per_compound_errors.csv')
    diffs=[];dc=[]
    comparisons=[('measured','D_meas','D_struct'),('predicted','D_pred','D_struct'),('substitution','D_swap','D_meas'),('measured_only','D_meas_only','D_cond'),('permutation','D_meas','D_perm'),('meanZ_ablation','D_meanZ','D_meas'),('swap_vs_meanZ','D_swap','D_meanZ')]+[(f'meas_vs_{b}','D_meas',b) for b in ['zero_response','training_compound_mean','training_weighted_median','D_cond','D_meas_only']]
    for name,a,b in comparisons:
        for scope,h in met.groupby('scope'):
            h=h.set_index('route_model');diffs.append(dict(comparison=name,model_a=a,model_b=b,scope=scope,delta_MAE=h.loc[a,'MAE_cmp']-h.loc[b,'MAE_cmp'],delta_RMSE=h.loc[a,'RMSE_cmp']-h.loc[b,'RMSE_cmp']))
        arows=oo.loc[oo.route_model.eq(a)].set_index('condition_id');brows=oo.loc[oo.route_model.eq(b)].set_index('condition_id').loc[arows.index]
        d=arows[['parent_id','source_name','outer_fold','y']].copy();d['comparison']=name;d['delta_absolute_error']=arows.absolute_error-brows.absolute_error;d['delta_squared_error']=arows.squared_error-brows.squared_error;dc.append(d.reset_index())
    csv(diffs,E/'modeling/comparisons.csv');dc=pd.concat(dc,ignore_index=True);csv(dc,E/'modeling/per_condition_differences.csv');csv(dc.groupby(['comparison','parent_id','source_name','outer_fold'],as_index=False).agg(delta_MAE=('delta_absolute_error','mean'),delta_MSE=('delta_squared_error','mean')),E/'modeling/per_parent_differences.csv')
    csv(pd.DataFrame(selected).groupby(['route_model','lambdaX','lambdaZ'],dropna=False).size().reset_index(name='n_outer_folds'),E/'modeling/lambda_DROP_frequency.csv')
    assert len(fitrows)==531 and len(manifest)==18 and len(pc.manifest)==27 and len(pc.contexts)==39
    complete=dict(status='completed_one_finite_family_stopped_before_holdout',completed_utc=now(),ridge_fits=len(fitrows),inner_ridge_fits=513,outer_ridge_fits=18,extra_proxy_fits=len(pc.manifest),legacy_proxy_reuses=12,ablation_fits=0,full_dev_production_fits=0,holdout_prediction_calls=0,n_dev_parents=len(allowed),n_dev_conditions=len(t),routes_scored=11)
    dump(complete,E/'modeling/completed.json');log('family_completed',**complete);print(met.loc[met.scope.eq('pooled_OOF')].to_string(index=False))
