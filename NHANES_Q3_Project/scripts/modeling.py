"""Prespecified CPU-only estimators with training-local transforms and nested OOF."""
from runtime import limit_threads
limit_threads()
import time
import warnings
import numpy as np
import pandas as pd
from sklearn.preprocessing import OneHotEncoder, StandardScaler, SplineTransformer
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.ensemble import HistGradientBoostingRegressor, HistGradientBoostingClassifier
from sklearn.exceptions import ConvergenceWarning

class FeatureMap:
    def __init__(self, cfg, z_present=False, spline=False, standardize=False):
        self.numeric = list(cfg['features']['numeric_X'])
        self.categories = list(cfg['features']['categorical_X'])
        self.zcols = list(cfg['features']['measured_Z']) if z_present else []
        self.use_spline = spline
        self.standardize = standardize
        self.spline_defaults = cfg['downstream']['spline_defaults']

    def _raw(self, frame):
        # Deliberate whitelist: metadata and mortality fields cannot enter X.
        n = frame[self.numeric].to_numpy(float)
        missing = ~np.isfinite(n)
        values = np.where(missing, self.medians[None,:], n)
        if self.zcols:
            z = frame[self.zcols].to_numpy(float)
            if not np.isfinite(z).all():
                raise ValueError('Z must be observed/predicted in original units; no Z imputation.')
            values = np.column_stack([values,z])
        cat = frame[self.categories].astype('string').fillna('__MISSING__').to_numpy(str)
        return values,missing.astype(float),cat

    def fit(self, frame):
        self.feature_names_in_ = self.numeric+self.categories+self.zcols
        n = frame[self.numeric].to_numpy(float)
        self.medians = np.nanmedian(np.where(np.isfinite(n),n,np.nan),axis=0)
        if not np.isfinite(self.medians).all():
            raise ValueError('A numeric training column has no finite median; diagnose rather than silently invent an imputation.')
        values,ind,cat=self._raw(frame)
        self.encoder=OneHotEncoder(handle_unknown='ignore',sparse_output=False,dtype=np.float64).fit(cat)
        if self.use_spline:
            d=self.spline_defaults
            self.spline=SplineTransformer(n_knots=d['n_knots'],degree=d['degree'],
                knots=d['knots'],include_bias=d['include_bias'],extrapolation=d['extrapolation']).fit(values)
        a=self._expanded(values,ind,cat)
        self.scaler=StandardScaler().fit(a) if self.standardize else None
        self.n_features_out_=a.shape[1]
        self.fit_n_=len(frame)
        return self

    def _expanded(self,values,ind,cat):
        numeric=self.spline.transform(values) if self.use_spline else values
        return np.column_stack([numeric,ind,self.encoder.transform(cat)])

    def transform(self,frame):
        a=self._expanded(*self._raw(frame))
        return self.scaler.transform(a) if self.scaler is not None else a

def weighted_targets(z,w,cfg):
    mean=np.average(z,axis=0,weights=w)
    var=np.average((z-mean)**2,axis=0,weights=w)
    if np.any(var<=cfg['proxy']['minimum_training_target_variance']):
        raise ValueError('Training target variance at/below prespecified numerical threshold.')
    return mean,np.sqrt(var)

class ProxyModel:
    def __init__(self,cfg,candidate):
        self.cfg=cfg
        self.candidate=candidate

    def fit(self,frame,z,weights):
        w=np.asarray(weights,float)
        w=w/w.mean()
        self.mean_,self.scale_=weighted_targets(z,w,self.cfg)
        self.features=FeatureMap(self.cfg,standardize=self.candidate['model']=='Ridge').fit(frame)
        a=self.features.transform(frame)
        target=(z-self.mean_)/self.scale_
        d=self.candidate
        if d['model']=='Ridge':
            self.estimator=Ridge(alpha=d['alpha'],solver=d['solver']).fit(a,target,sample_weight=w)
        else:
            kwargs={k:d[k] for k in ['loss','learning_rate','max_iter','max_leaf_nodes','min_samples_leaf','l2_regularization','early_stopping']}
            self.estimator=[HistGradientBoostingRegressor(**kwargs,categorical_features=None,
                random_state=self.cfg['random_state']).fit(a,target[:,j],sample_weight=w) for j in range(target.shape[1])]
        return self

    def predict(self,frame):
        a=self.features.transform(frame)
        p=np.column_stack([m.predict(a) for m in self.estimator]) if isinstance(self.estimator,list) else self.estimator.predict(a)
        return p*self.scale_+self.mean_

class OutcomeModel:
    def __init__(self,cfg,candidate,z_present):
        self.cfg=cfg
        self.candidate=candidate
        self.z_present=z_present

    def fit(self,frame,y,weights):
        d=self.candidate
        spline=d['family']=='spline_logistic'
        self.features=FeatureMap(self.cfg,z_present=self.z_present,spline=spline,standardize=spline).fit(frame)
        a=self.features.transform(frame)
        w=np.asarray(weights,float); w=w/w.mean()
        self.convergence_messages_=[]
        if spline:
            defaults=self.cfg['downstream']['spline_defaults']
            # Default l1_ratio=0 is L2 in sklearn >=1.8; default penalty=L2 in older supported versions.
            self.estimator=LogisticRegression(C=d['C'],solver=defaults['logistic_solver'],
                max_iter=defaults['logistic_max_iter'],class_weight=None,random_state=self.cfg['random_state'])
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter('always',ConvergenceWarning)
                self.estimator.fit(a,y,sample_weight=w)
            self.convergence_messages_=[str(v.message) for v in caught if issubclass(v.category,ConvergenceWarning)]
            if self.convergence_messages_:
                # Same estimator/penalty, only extend the authorized iteration limit.
                self.estimator.set_params(max_iter=12000)
                with warnings.catch_warnings():
                    warnings.simplefilter('error',ConvergenceWarning)
                    self.estimator.fit(a,y,sample_weight=w)
        else:
            kwargs=self.cfg['downstream']['hgb_defaults'].copy()
            self.estimator=HistGradientBoostingClassifier(**kwargs,max_leaf_nodes=d['max_leaf_nodes'],
                l2_regularization=d['l2_regularization'],random_state=self.cfg['random_state']).fit(a,y,sample_weight=w)
        return self

    def predict_proba(self,frame):
        return self.estimator.predict_proba(self.features.transform(frame))[:,1]

def group_folds(frame,seed,n=3):
    groups=np.array(sorted(frame.psu_id.astype(str).unique()),dtype=object)
    if len(groups)<n:
        raise ValueError('Insufficient PSU groups for the fixed fold count.')
    np.random.default_rng(seed).shuffle(groups)
    mapping={str(g):int(i) for i,part in enumerate(np.array_split(groups,n)) for g in part}
    fold=frame.psu_id.astype(str).map(mapping).to_numpy(int)
    return fold,mapping

def pick_candidate(scores,tolerance=1e-12):
    best=0
    for i in range(1,len(scores)):
        if scores[i]<scores[best]-tolerance:
            best=i
    return best

def timed_fit(model,frame,target,weights,label,callback=None):
    if callback:
        callback('start',label,frame,model,None)
    start=time.perf_counter()
    model.fit(frame,target,weights)
    elapsed=time.perf_counter()-start
    if callback:
        callback('finish',label,frame,model,elapsed)
    return model

def select_proxy(frame,cfg,pool_k,label,callback=None):
    fold,mapping=group_folds(frame,cfg['random_state'],cfg['proxy']['inner_cv_folds'])
    z=frame[cfg['features']['measured_Z']].to_numpy(float)
    w=frame.WTMEC2YR.to_numpy(float)/pool_k
    rows=[]; scores=[]
    for candidate in cfg['proxy']['candidates']:
        parts=[]
        for k in range(cfg['proxy']['inner_cv_folds']):
            train=fold!=k; val=~train
            model=timed_fit(ProxyModel(cfg,candidate),frame.loc[train],z[train],w[train],
                f'{label}/{candidate["id"]}/inner{k}',callback)
            pred=model.predict(frame.loc[val])
            mse=np.average((pred-z[val])**2,axis=0,weights=w[val])
            scaled=mse/model.scale_**2
            score=float(np.mean(scaled))
            parts.append((score,float(w[val].sum())))
            rows.append(dict(candidate=candidate['id'],fold=k,score=score,
                per_target_scaled_MSE=scaled.tolist(),training_target_variance=(model.scale_**2).tolist(),
                validation_weight_total=float(w[val].sum()),
                fit_psus=sorted(frame.loc[train,'psu_id'].unique().tolist()),
                heldout_psus=sorted(frame.loc[val,'psu_id'].unique().tolist())))
        scores.append(float(np.average([p[0] for p in parts],weights=[p[1] for p in parts])))
    chosen=pick_candidate(scores,cfg['proxy']['tie_absolute_tolerance'])
    record=dict(label=label,pool_K=pool_k,group_to_fold=mapping,fold_records=rows,
                candidates=[dict(candidate=c['id'],score=s) for c,s in zip(cfg['proxy']['candidates'],scores)],
                selected=cfg['proxy']['candidates'][chosen]['id'],selection_uses_Y=False)
    return cfg['proxy']['candidates'][chosen],record

def replace_z(frame,z,cfg):
    out=frame.copy()
    out.loc[:,cfg['features']['measured_Z']]=np.asarray(z)
    return out
