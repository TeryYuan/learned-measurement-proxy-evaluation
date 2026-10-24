#!/usr/bin/env python3
"""Export an explicit, source-only allowlist to a fresh directory.

This prepares reviewable repository contents. It does not run git, publish,
fit models, or copy study data, predictions, caches or manuscript context.
"""
from pathlib import Path
import argparse
import json
import shutil

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = [
    '.gitignore', 'requirements.txt', 'requirements.lock.txt',
    'protocol/PREMODEL_ANALYSIS_PLAN.md', 'protocol/PROTOCOL_RECORD.json',
    'protocol/PROTOCOL_zh.md', 'protocol/analysis_config.json',
    'protocol/feature_dictionary.json', 'protocol/STATISTICS.md',
    'protocol/CBC_INSTRUMENT_AUDIT.md', 'protocol/LMF_WEIGHT_AUDIT.md',
    'protocol/INTERPRETATION_zh.md', 'protocol/CHANGES_V2.md',
    'protocol/DEVIATIONS.md',
    'protocol/EXECUTION_CADENCE_AMENDMENT.json',
    'scripts/prepare_cohort.py', 'scripts/evaluate_fixed_predictions.py',
    'scripts/check_statistics.py', 'scripts/fetch_inputs.py',
    'scripts/pre_model_descriptives.py', 'scripts/export_code_only.py',
    'scripts/README.md', 'sources/SOURCES.md', 'sources/source_manifest.json',
]
# The official R reader is intentionally not redistributed; Python preparation
# uses the official URLs documented in sources/SOURCES.md and needs no R file.
# Add actual reproduction modules here when implemented; never add input/, data/,
# predictions, trained model objects, or unpublished manuscript files.
OPTIONAL_FUTURE_SOURCE = [
    'scripts/modeling.py', 'scripts/run_benchmark.py',
    'scripts/make_outputs.py', 'scripts/run_module.py',
    'scripts/modeling.py', 'scripts/runtime.py', 'scripts/preflight.py',
    'scripts/check_modeling.py', 'scripts/initialize_execution.py',
    'scripts/record_local_authorization.py',
    'scripts/saved_reports.py', 'scripts/render_saved_figures.py',
]


def export(destination):
    destination = Path(destination).expanduser().resolve()
    if destination.exists():
        raise FileExistsError('Use a fresh export destination; existing work is not overwritten.')
    if destination == ROOT or ROOT in destination.parents:
        raise ValueError('Use a sibling destination outside the local working kit.')
    missing = [name for name in REQUIRED if not (ROOT / name).is_file()]
    readme = ROOT / 'README_PUBLIC.md'
    if not readme.is_file():
        # The exported source tree already has its public README at this path.
        readme = ROOT / 'README.md'
    if missing or not readme.is_file():
        raise FileNotFoundError(f'Missing required source files: {missing}; public README present: {readme.is_file()}')
    names = list(dict.fromkeys(REQUIRED + [n for n in OPTIONAL_FUTURE_SOURCE if (ROOT / n).is_file()]))
    sources = [(ROOT / name, name) for name in names] + [(readme, 'README.md')]
    for source, _ in sources:
        if source.is_symlink() or ROOT not in source.resolve().parents:
            raise ValueError(f'Unexpected symlink or external source: {source.name}')
    destination.mkdir(parents=True)
    for source, relative in sources:
        dest = destination / relative
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
    report = {
        'status': 'source_export_prepared_not_published',
        'files': sorted(relative for _, relative in sources),
        'excluded': ['raw inputs', 'participant-level prepared data or predictions',
                     'trained models and caches', 'unpublished manuscript context',
                     'journal-specific manuscript integration drafts'],
        'git_operations': 0, 'study_model_fits': 0,
    }
    (destination / 'EXPORT_MANIFEST.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    return {'directory': str(destination), 'files': len(sources) + 1, 'published': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path,
                        default=ROOT.parent / 'nhanes_q3_code_only_v2_20261019')
    args = parser.parse_args()
    print(json.dumps(export(args.destination)))


if __name__ == '__main__':
    main()
