"""Small local execution records; no uploads or model fitting."""
from pathlib import Path
from datetime import datetime, timezone
import json
import os

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs'

def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')

def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))

def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    tmp.replace(path)

def event(kind, **fields):
    path = OUT / 'logs/execution.jsonl'
    path.parent.mkdir(parents=True, exist_ok=True)
    item = {"utc":now(), "event":kind, **fields}
    with path.open('a', encoding='utf-8') as f:
        f.write(json.dumps(item, ensure_ascii=False, allow_nan=False)+'\n')
    return item

def require_remote_record():
    rec = read_json(ROOT/'protocol/PROTOCOL_RECORD.json')
    waiver=rec.get('author_authorized_local_execution')
    if waiver and waiver.get('instruction') and waiver.get('recorded_utc'):
        # The recorded author decision permits local execution without a remote protocol record.
        # This is not a remote record and must never populate its URL/time fields.
        return rec
    keys = ['repository_url', 'protocol_commit_url', 'repository_visibility',
            'remote_recorded_utc', 'remote_record_evidence']
    if any(not rec.get(k) for k in keys):
        raise RuntimeError('AUTHOR PREREQUISITE: verified remote protocol record is missing; no study-data fit allowed.')
    if not rec['protocol_commit_url'].startswith(rec['repository_url'].rstrip('/')+'/commit/'):
        raise RuntimeError('Protocol commit must belong to the designated repository.')
    stamp = datetime.fromisoformat(rec['remote_recorded_utc'])
    if stamp.tzinfo is None or stamp > datetime.now(timezone.utc):
        raise RuntimeError('Remote evidence time must be an actual past timezone-aware time.')
    return rec

def mark_first_fit():
    rec = require_remote_record()
    if rec.get('first_study_fit_utc') is None:
        rec['first_study_fit_utc'] = now()
        rec['status'] = ('author_authorized_local_plan_study_fitting_started' if rec.get('author_authorized_local_execution')
                         else 'remote_protocol_record_verified_study_fitting_started')
        write_json(ROOT/'protocol/PROTOCOL_RECORD.json', rec)
        event('first_study_fit_start', protocol_commit_url=rec['protocol_commit_url'])

def limit_threads():
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS']:
        os.environ[key] = '4'

