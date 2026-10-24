# NHANES clinical reproduction module

Current instructions for this collection. The [script guide](scripts/README.md) describes the current interfaces, inputs and outputs. [Historical source-export provenance](README_original_source_export.md) explains the smaller earlier export.

From this directory:
```sh
python -m pip install -r requirements.lock.txt
python scripts/render_saved_figures.py --report-only --report-dir saved_figures
python scripts/run_module.py --report-only --report-dir saved_report
python scripts/run_benchmark.py --outputs-only --report-dir saved_benchmark_report
```
The drawing entry reads `reference/outputs/figures/Figure6_source_data.json` and `reference/diagnostics/cbc_distributions_by_cycle.svg`. It rebuilds the historical 167 mm Figure 6/S12 layout and saves drawing source/layout records. The current 174 mm publication presentation is in [figure_patch/](../figure_patch/README.md), as described in [reproduction scope](../REPRODUCTION_SCOPE.md). S12 reuses displayed curve geometry, not person-level ECDF recalculation. The broader report entries use saved aggregate tables and Figure 6. They never enter fitting or statistical prediction evaluation. When an execution log exists it is guarded against modification; omission in this public candidate is recorded as absence, with no fabricated log.

For full reconstruction, use a separate writable copy, obtain the 40 official files from the recorded manifest, then run:
```sh
python scripts/fetch_inputs.py
python scripts/prepare_cohort.py
python scripts/run_module.py
```
Acquisition requires network access; fitting can run offline once official inputs are present. The original local protocol record and its absence of pre-model remote registration are unchanged. New output is written to `outputs/`; reference aggregates under `reference/` remain unchanged. Every successful benchmark invocation is followed by evaluation of its current predictions and fresh output generation, regardless of previous summary or completion markers. Benchmark cache compatibility checks remain in effect.

See [validation scope](../VALIDATION.md) for the earlier isolated full reconstruction, including all 24 candidates, and its limits. Current publication preparation does not rerun that training. Current release state is recorded in [RELEASE_STATUS.json](../RELEASE_STATUS.json). This is NHANES engineering reproduction under the existing protocol, not a new independent scientific study. The original observed dependency record is `requirements_observed_historical.txt`; the tested compatible lock uses NumPy 2.5.3, retaining the other selected scientific versions. CPU numerical thread limits are capped at four without changing candidates, seeds, expected values or tolerances.

`export_code_only.py --destination NEW_SIBLING_DIRECTORY` exports source and protocol files only; it does not include aggregate reference inputs, figures or the full distribution dependencies. Saved-output commands in that smaller export require supplying the documented aggregate files. The official R reader is omitted, is never required for re-export, and is not used by Python preparation. Obtain it only from the official URL in [sources/SOURCES.md](sources/SOURCES.md) if needed independently.
