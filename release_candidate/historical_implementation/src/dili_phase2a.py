"""Fixed Phase 2A query inference and paired statistics; no training entry point."""
import os
for _key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
    os.environ[_key] = '1'
import pathlib, json, datetime, contextlib, collections, math
from unittest.mock import patch
import numpy as np
import pandas as pd
import joblib
from rdkit import DataStructs
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import average_precision_score, roc_auc_score, matthews_corrcoef
import scipy.optimize, scipy.linalg
import dili_phase1 as old
import dili_phase1p5 as fusion
from dili_phase1_log import append_jsonl

ROOT = pathlib.Path(__file__).resolve().parents[1]
LOG = ROOT / 'logs/phase2a_N'
MODELS = ['M0', 'M1', 'Mperm', 'Mchem', 'B_star', 'L2LR', *fusion.ROUTES]
TASKS = old.TASKS
QUERY_COLUMNS = ['member_id', 'model_input_smiles', 'fold_group']
METRICS = ['AP', 'ROC_AUC', 'Brier_raw', 'Brier_calibrated']


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def event(kind, **fields):
    append_jsonl(LOG / 'events.jsonl', {'utc': utc(), 'event': kind, **fields})


@contextlib.contextmanager
def no_learning(scope, record=True):
    """Fail closed at estimator fit and optimization APIs within this process."""
    attempts = collections.Counter()
    targets = [(c, m) for c in [RandomForestClassifier, RandomForestRegressor, LogisticRegression,
                                StandardScaler, Pipeline] for m in ['fit', 'fit_transform', 'partial_fit'] if hasattr(c, m)]
    targets += [(old, m) for m in ['fit_forest', 'calibrate', 'select_threshold', 'choose_grid', 'assign_folds']]
    targets += [(fusion, m) for m in ['fit_residual', 'spectral', 'choose_lambda', 'minimize', 'eigh']]
    targets += [(scipy.optimize, 'minimize'), (scipy.linalg, 'eigh'), (np.linalg, 'eigh')]
    with contextlib.ExitStack() as stack:
        for owner, method in targets:
            name = getattr(owner, '__name__', str(owner)) + '.' + method
            def denied(*args, _name=name, **kwargs):
                attempts[_name] += 1
                if record:
                    event('prohibited_learning_attempt', scope=scope, operation=_name)
                raise PermissionError('Phase2A fixed inference prohibits ' + _name)
            stack.enter_context(patch.object(owner, method, denied))
        if record:
            event('no_learning_guard_started', scope=scope, intercepted_apis=len(targets))
        try:
            yield attempts
        finally:
            if record:
                event('no_learning_guard_finished', scope=scope, attempted_calls=dict(attempts),
                      successful_research_fits=0, intercepted_apis=len(targets))


def positive_probability(model, x):
    classes = np.asarray(model.classes_)
    if set(classes.tolist()) != {0, 1}:
        raise ValueError('Expected both binary classes 0 and 1')
    return model.predict_proba(x)[:, int(np.flatnonzero(classes == 1)[0])]


def forest_equal(a, b):
    if a.get_params() != b.get_params() or not np.array_equal(a.classes_, b.classes_):
        return False
    if len(a.estimators_) != len(b.estimators_):
        return False
    for ta, tb in zip(a.estimators_, b.estimators_):
        aa, bb = ta.tree_.__getstate__(), tb.tree_.__getstate__()
        if aa.keys() != bb.keys() or any(not np.array_equal(aa[k], bb[k]) for k in aa):
            return False
    return True


def bits(array):
    return [DataStructs.CreateFromBitString(''.join(map(str, row))) for row in array]


class FixedInference:
    """Only structures/IDs enter query(); models are loaded exclusively from final contexts."""
    def __init__(self, root=ROOT):
        self.root = pathlib.Path(root)
        self.cfg1 = json.loads((self.root / 'config/phase1_execution.json').read_text(encoding='utf8'))
        self.cfg15 = json.loads((self.root / 'config/phase1p5_execution.json').read_text(encoding='utf8'))
        self.o, self.ob, self.od = old.load_dataset('O')
        self.ids = self.o.member_id.astype(str).tolist()
        self.oq = pd.read_csv(self.root / 'data/phase1p5/O_coverage.csv').set_index('member_id').loc[self.ids, ['q_' + t for t in TASKS]].to_numpy()
        self.oz = pd.read_csv(self.root / 'data/phase1/features/O_proxy_features.csv').set_index('member_id').loc[self.ids, ['z_' + t for t in TASKS]].to_numpy()
        self.models, self.aux, self.a, self.leaves, self.paths = {}, {}, {}, {}, []
        self.manifest = []
        for task in TASKS:
            ai, ab, ad = old.load_dataset('A_' + task)
            self.a[task] = (ai, ab, bits(ab))
            self.paths += [self.root / f'data/phase1/inputs/A_{task}_{suffix}' for suffix in ['index.csv', 'base.npz']]
            for branch in ['real', 'perm', 'chem']:
                path = self.root / f'models/phase1/auxiliary/{task}_{branch}.joblib'
                model = joblib.load(path)
                meta = json.loads(path.with_suffix('.json').read_text(encoding='utf8'))
                assert meta['train_ids'] == ai.member_id.tolist()
                assert model.n_features_in_ == 2048 and len(model.estimators_) == 500
                assert 1 <= model.n_jobs <= 4
                if branch != 'chem':
                    assert set(model.classes_) == {0, 1}
                self.aux[task, branch] = model
                lp = self.root / f'data/phase1p5/O_leaf_indices/{task}_{branch}.npz'
                with np.load(lp) as f:
                    assert f['ids'].tolist() == self.ids
                    self.leaves[task, branch] = f['leaves'].copy()
                self.paths += [path, path.with_suffix('.json'), lp]
        for name in MODELS:
            phase = 'phase1/dili' if name in old.MODELS else 'phase1p5'
            path = self.root / f'models/{phase}/{name}_final.joblib'
            obj = joblib.load(path)
            assert obj['task_order'] == TASKS
            assert obj['calibrator'] is not None and set(obj['calibrator'].classes_) == {0, 1}
            if name in old.MODELS:
                meta = obj['metadata']; model = obj['rf']
                assert meta['context'] == 'final' and not meta['validation_ids'] and meta['train_ids'] == self.ids
                assert obj['feature_schema'] == self.cfg1['features']
                assert model.n_features_in_ == (2056 if name == 'M0' else 2059)
                assert set(model.classes_) == {0, 1} and len(model.estimators_) == 500
            else:
                assert obj['model_id'] == name + '_final' and obj['train_ids'] == self.ids and not obj['query_ids']
                assert obj['config'] == self.cfg15 and obj['basis_reference_ids'] == self.ids
                if name.startswith('R_'):
                    assert obj['config']['residual']['epsilon'] == fusion.EPS
                    assert obj['basis']['vectors'].shape == (len(self.ids), len(obj['weights']))
                    assert np.all(obj['basis']['values'] > 0)
                else:
                    assert obj['model'].n_features_in_ == 2056 and set(obj['model'].classes_) == {0, 1}
            self.models[name] = obj
            self.paths.append(path)
            self.manifest.append({'model': name, 'source_phase': phase, 'path': str(path.relative_to(self.root)),
                'context': 'final_full_O', 'training_n': len(self.ids), 'base_dimension': 2056,
                'input_dimension': 2059 if name in ['M1', 'Mperm', 'Mchem'] else 2056,
                'descriptors': '|'.join(self.cfg1['features']['descriptors']), 'task_order': '|'.join(TASKS),
                'raw_definition': 'offset_logistic_residual' if name.startswith('R_') else 'positive_class_probability',
                'raw_threshold': obj['threshold'], 'threshold_operator': '>=', 'calibrator': type(obj['calibrator']).__name__,
                'calibration_coef_exact': repr(obj['calibrator'].coef_.tolist()),
                'calibration_intercept_exact': repr(obj['calibrator'].intercept_.tolist()),
                'lambda': obj.get('selected_lambda'), 'rank': len(obj.get('weights', [])),
                'epsilon': obj.get('config', {}).get('residual', {}).get('epsilon'),
                'reference_identity_file': 'data/phase1/inputs/O_index.csv',
                'reference_leaf_dir': obj.get('O_leaf_reference_dir', ''),
                'basis_and_weights': str(path.relative_to(self.root)) if name.startswith('R_') else '',
                'auxiliary_branch': {'M1':'real','Mperm':'perm','Mchem':'chem','R_gate_perm':'perm','R_gate_chem':'chem'}.get(name, 'real' if name in ['R_score','R_leaf','R_gate'] else 'none')})
        self.backbone_path = self.root / self.models['R_gate']['backbone_path']
        bb = joblib.load(self.backbone_path)
        assert bb['expected']['train_ids'] == self.ids and not bb['expected']['query_ids']
        self.backbone = bb['model']
        assert forest_equal(self.models['B_star']['model'], self.backbone)
        assert all(self.models[r]['backbone_path'] == self.models['R_gate']['backbone_path'] for r in fusion.ROUTES)
        self.alias_raw = forest_equal(self.models['M0']['rf'], self.models['B_star']['model'])
        m0, bs = self.models['M0'], self.models['B_star']
        self.alias_full = self.alias_raw and m0['threshold'] == bs['threshold'] and np.array_equal(m0['calibrator'].coef_, bs['calibrator'].coef_) and np.array_equal(m0['calibrator'].intercept_, bs['calibrator'].intercept_)
        for row in self.manifest:
            row['alias_of'] = 'M0' if row['model'] == 'B_star' and self.alias_full else ''
        self.obits = bits(self.ob)
        self.paths += [self.backbone_path, self.root / 'data/phase1/inputs/O_index.csv', self.root / 'data/phase1/inputs/O_base.npz',
                       self.root / 'data/phase1p5/O_coverage.csv', self.root / 'data/phase1/features/O_proxy_features.csv']

    def query(self, structures, *, role):
        if list(structures.columns) != QUERY_COLUMNS:
            raise ValueError('Query requires only ID, fixed structure and group; labels/extra columns forbidden')
        if role not in ['N', 'O_interface_check']:
            raise PermissionError('Only released N and small O interface checks are authorized')
        ids = structures.member_id.astype(str).tolist()
        if len(ids) != len(set(ids)) or not len(ids):
            raise ValueError('Unique nonempty query IDs required')
        if role == 'N':
            release = json.loads((self.root / 'config/phase2a_N_execution.json').read_text(encoding='utf8'))
            assert release['allow_fixed_N_query'] and not release['allow_research_fit']
            allowed = pd.read_csv(self.root / 'data/phase2a_N/N_structure_queries.csv').set_index('member_id')
        else:
            allowed = self.o.set_index('member_id')
        if not set(ids) <= set(allowed.index) or not np.array_equal(structures.model_input_smiles, allowed.loc[ids, 'model_input_smiles']):
            raise PermissionError('Query outside fixed released structures')
        event('query_started', role=role, ids=ids, n=len(ids))
        encoded = [old.structural_features(s) for s in structures.model_input_smiles]
        b = np.stack([v[0] for v in encoded]); d = np.stack([v[1] for v in encoded]); x = np.column_stack([b, d])
        query_bits = bits(b); q = np.zeros((len(ids), 3)); nn = []; leaf = {}; parts = {}; proxy = {'member_id': ids}
        for k, task in enumerate(TASKS):
            ai, ab, abits = self.a[task]; aids = ai.member_id.to_numpy(str)
            for i, bit in enumerate(query_bits):
                sim = np.asarray(DataStructs.BulkTanimotoSimilarity(bit, abits)); order = np.lexsort((aids, -sim))[:5]
                q[i, k] = sim[order].mean()
                nn.append({'member_id': ids[i], 'reference_role': 'A_' + task, 'reference_n': len(ai),
                    'q_top5': q[i,k], 'max_Tanimoto': sim[order[0]], 'neighbor_ids': '|'.join(aids[order]),
                    'neighbor_Tanimoto': '|'.join(map(str, sim[order])),
                    'same_Murcko_group_seen': str(structures.fold_group.iloc[i]) in set(ai.fold_group)})
            for branch, prefix in [('real', 'z'), ('perm', 'z_perm'), ('chem', 'z_chem')]:
                model = self.aux[task, branch]
                proxy[prefix + '_' + task] = model.predict(b) if branch == 'chem' else positive_probability(model, b)
                leaf[task, branch] = model.apply(b).astype(np.int32)
                parts[task, branch] = fusion.leaf_kernel(leaf[task, branch], self.leaves[task, branch])
                event('fixed_auxiliary_predict_apply', role=role, task=task, branch=branch, ids=ids, reference_role='fixed_A_' + task)
        for i, bit in enumerate(query_bits):
            sim = np.asarray(DataStructs.BulkTanimotoSimilarity(bit, self.obits)); order = np.lexsort((np.asarray(self.ids), -sim))[:5]
            nn.append({'member_id':ids[i], 'reference_role':'O', 'reference_n':len(self.ids), 'q_top5':None,
                'max_Tanimoto':sim[order[0]], 'neighbor_ids':'|'.join(np.asarray(self.ids)[order]),
                'neighbor_Tanimoto':'|'.join(map(str,sim[order])),
                'same_Murcko_group_seen':str(structures.fold_group.iloc[i]) in set(self.o.fold_group)})
        z = np.column_stack([proxy['z_' + t] for t in TASKS])
        kernels = {'R_score': z @ self.oz.T / 3, 'R_leaf': sum(parts[t,'real'] for t in TASKS) / 3,
                   'R_gate_struct': (q @ self.oq.T / 3) * fusion.tanimoto(b, self.ob)}
        for branch, route in [('real','R_gate'),('perm','R_gate_perm'),('chem','R_gate_chem')]:
            kernels[route] = sum(q[:,k,None] * self.oq[None,:,k] * parts[t,branch] for k,t in enumerate(TASKS)) / 3
        base = positive_probability(self.backbone, x)
        projections = {}; rows = []
        for name in MODELS:
            obj = self.models[name]; h = np.zeros(len(ids))
            if name in old.MODELS:
                prefix = {'M1':'z','Mperm':'z_perm','Mchem':'z_chem'}.get(name)
                inp = x if prefix is None else np.column_stack([x, *[proxy[prefix + '_' + t] for t in obj['task_order']]])
                raw = positive_probability(obj['rf'], inp)
            elif name in ['B_star','L2LR']:
                raw = positive_probability(obj['model'], x)
            else:
                assert kernels[name].shape == (len(ids), len(obj['basis_reference_ids']))
                projections[name] = fusion.project(kernels[name], obj['basis'])
                raw, h = fusion.residual_predict(base, projections[name], obj['weights'], obj['selected_lambda'])
            cal = positive_probability(obj['calibrator'], raw.reshape(-1,1))
            assert np.isfinite(raw).all() and np.isfinite(cal).all() and ((raw>=0)&(raw<=1)).all()
            for i, identity in enumerate(ids):
                rows.append({'member_id':identity, 'model':name, 'raw_score':raw[i], 'calibrated_score':cal[i],
                    'raw_threshold':obj['threshold'], 'predicted_class':int(raw[i]>=obj['threshold']),
                    'backbone_raw':base[i], 'logit_correction':h[i] if name.startswith('R_') else None,
                    'probability_correction':raw[i]-np.clip(base[i],fusion.EPS,1-fusion.EPS) if name.startswith('R_') else None,
                    'q_mean':q[i].mean(), 'failure_status':'ok', 'failure_reason':''})
            event('fixed_model_prediction', role=role, model=name, context='final_full_O', ids=ids,
                  spectral_query_shape=list(kernels[name].shape) if name.startswith('R_') else None)
        return {'predictions':pd.DataFrame(rows), 'b':b, 'd':d, 'proxy':pd.DataFrame(proxy),
                'coverage':pd.DataFrame({'member_id':ids, **{'q_'+t:q[:,k] for k,t in enumerate(TASKS)}}),
                'neighbors':pd.DataFrame(nn), 'leaves':leaf, 'kernels':kernels, 'projections':projections}


def merge_outcomes(predictions, labels, model_names=MODELS):
    if 'y_v2' in predictions or labels.member_id.duplicated().any() or predictions.duplicated(['model','member_id']).any():
        raise ValueError('Outcomes must be joined once with unique IDs')
    if set(predictions.model) != set(model_names):
        raise ValueError('Incomplete model coverage')
    for name, group in predictions.groupby('model'):
        if set(group.member_id) != set(labels.member_id) or not group.failure_status.eq('ok').all():
            raise ValueError('Cannot evaluate an incomplete or differing N subset')
    return predictions.merge(labels, on='member_id', how='left', validate='many_to_one')


def bootstrap_indices(ids, y, replicates=2000, seed=4109):
    ids=np.asarray(ids,str); y=np.asarray(y,int)
    if ids.tolist()!=sorted(ids.tolist()) or len(set(ids))!=len(ids):
        raise ValueError('Unique stable ascending IDs required')
    pos=np.flatnonzero(y==1);neg=np.flatnonzero(y==0)
    if not len(pos) or not len(neg):raise ValueError('Both classes required for stratified bootstrap')
    rng=np.random.Generator(np.random.PCG64(seed));out=np.empty((replicates,len(ids)),np.int32)
    for i in range(replicates):out[i]=np.r_[rng.choice(pos,len(pos),replace=True),rng.choice(neg,len(neg),replace=True)]
    return out


def metrics4(y, raw, calibrated):
    y=np.asarray(y); raw=np.asarray(raw);calibrated=np.asarray(calibrated);both=len(np.unique(y))==2
    return {'AP':float(average_precision_score(y,raw)) if both else np.nan,
            'ROC_AUC':float(roc_auc_score(y,raw)) if both else np.nan,
            'Brier_raw':float(np.mean((y-raw)**2)), 'Brier_calibrated':float(np.mean((y-calibrated)**2))}


def model_metrics(y, raw, cal, pred, ids):
    out=old.classification_metrics(y,raw,ids,cal,pred)
    out['MCC_library_return']=float(matthews_corrcoef(y,pred))
    return out


def percentile_interval(values, coverage):
    a=np.asarray(values,float);good=a[np.isfinite(a)];tail=(1-coverage)/2
    return {'low':float(np.quantile(good,tail,method='linear')) if len(good) else None,
            'high':float(np.quantile(good,1-tail,method='linear')) if len(good) else None,
            'valid_replicates':len(good),'invalid_replicates':len(a)-len(good),
            'NA_reason':'' if len(good)==len(a) else 'undefined statistic in resample'}
