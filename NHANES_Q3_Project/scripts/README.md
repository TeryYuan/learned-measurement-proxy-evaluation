# NHANES clinical module scripts

## Scope

This module evaluates learned RDW and hemoglobin proxies and substitution in one fitted mortality model under the existing NHANES analysis. The baseline has 16 specified variables; the endpoint is public-file-defined five-year all-cause mortality in adults aged 40 or older. The current collection includes saved aggregates under [reference/](../reference/) and the numerical source code. It omits participant-level inputs, prepared records, predictions, fitted objects and caches. See [the module README](../README.md) and [validation scope](../../VALIDATION.md).

Commands below run from `NHANES_Q3_Project/`, using the compatible [dependency lock](../requirements.lock.txt) in a separate writable environment. Historical observed versions remain in [requirements_observed_historical.txt](../requirements_observed_historical.txt).

## Script index

| Script | Interface and responsibility | Inputs | Outputs |
|---|---|---|---|
| [fetch_inputs.py](fetch_inputs.py) | Acquisition entry; `--list` lists sources without network requests | Source manifest and official files | Validated files at manifest paths |
| [prepare_cohort.py](prepare_cohort.py) | Deterministic cohort-preparation entry | 36 SAS XPT files, four public mortality files, field dictionary | `data/analysis_cohort.csv`, `data/design_frame.csv`, role/source checks |
| [run_module.py](run_module.py) | Orchestrates checks, benchmark, evaluation and outputs; separate `--report-only` branch | Prepared data/protocol for reconstruction; saved aggregates for reporting | Current `outputs/`, or requested report directory |
| [run_benchmark.py](run_benchmark.py) | Numerical fitting/selection/refit entry; `--outputs-only` branches before fitting imports | Prepared data, protocol and compatible stage caches; saved aggregates in reporting mode | Models, predictions and logs, or saved-result reports |
| [render_saved_figures.py](render_saved_figures.py) | Saved Figure 6/S12 rendering entry | Aggregate Figure 6 JSON and saved S12 SVG curves | Historical-layout PDF/SVG/PNG and drawing records |
| [saved_reports.py](saved_reports.py) | Shared reporting module; also accepts `--report-only` and `--report-dir` | Saved summary, Q1, figure JSON and aggregate CSV/JSON | `results_tables.md`, `source_data/`, Figure 6 and reuse record |
| [evaluate_fixed_predictions.py](evaluate_fixed_predictions.py) | Fixed-prediction evaluation entry; accepts `--predictions`, `--cohort`, `--design-frame`, `--output` | Five-route CSV and full survey-frame inputs | `outputs/tables/fixed_substitution_summary.json` by default |
| [make_outputs.py](make_outputs.py) | Output entry after current prediction evaluation | Current evaluated summaries and benchmark output | Tables and Figure 6 under `outputs/` |
| [pre_model_descriptives.py](pre_model_descriptives.py) | Descriptive entry; `--report-only` delegates to saved rendering | Prepared records in descriptive mode; saved aggregates/curves in report mode | Cycle/measurement diagnostics, or saved figures |
| [modeling.py](modeling.py) | Imported modeling module, not a CLI entry | Training subsets, explicit variables and candidates | Fitted preprocessors, proxies/outcome models, grouped OOF results |
| [runtime.py](runtime.py) | Imported runtime and record helpers | Paths, protocol and records | Thread limits, cache/record support and event logs when called |
| [preflight.py](preflight.py) | Preparation/schema check entry | Prepared cohort, survey frame and protocol | Preflight check record |
| [check_statistics.py](check_statistics.py) | Mathematical checks on artificial examples | Artificial inputs | Helper checks, not a study-result analysis |
| [check_modeling.py](check_modeling.py) | Modeling interface checks on artificial examples | Artificial inputs | Interface check record, not new study evidence |
| [initialize_execution.py](initialize_execution.py) | Historical one-time record utility | Existing protocol/cadence record | May write records if the amendment is absent; an existing amendment is retained |
| [record_local_authorization.py](record_local_authorization.py) | Historical one-time authorization-record utility | Existing author-authorization record | May write records if absent; not a reconstruction prerequisite for this completed collection |
| [export_code_only.py](export_code_only.py) | Source-subset export to a fresh external `--destination` | Allowlisted source, protocol, dependencies and documentation | Smaller source tree and `EXPORT_MANIFEST.json` |

## Rebuild reports from saved results

```sh
python scripts/render_saved_figures.py --report-only --report-dir saved_figures
python scripts/run_module.py --report-only --report-dir saved_report
python scripts/run_benchmark.py --outputs-only --report-dir saved_benchmark_report
```

The first reads [Figure6_source_data.json](../reference/outputs/figures/Figure6_source_data.json) and [cbc_distributions_by_cycle.svg](../reference/diagnostics/cbc_distributions_by_cycle.svg). Figure S12 reuses saved curve geometry, without recomputing person-level ECDFs. The other two call `saved_reports.generate`, which uses the saved tables and figure inputs under `reference/outputs/`. These branches do not fit models, query fitted objects or evaluate participant predictions. Choose destinations outside the saved reference tree. A missing execution log is treated as absent, not manufactured.

These commands reproduce the historical 167 mm module layout. The current 174 mm six-figure presentation, including the Figure 2/5 display edits, is in [figure_patch/](../../figure_patch/README.md). Its renderer uses saved drawing inputs. This distinction is recorded in [REPRODUCTION_SCOPE.md](../../REPRODUCTION_SCOPE.md).

## Reconstruct the numerical analysis

Use a separate writable copy and the [official source manifest](../sources/source_manifest.json):

```sh
python -m pip install -r requirements.lock.txt
python scripts/fetch_inputs.py
python scripts/prepare_cohort.py
python scripts/run_module.py
```

Acquisition requires network access; fitting can run offline once official inputs are available. Cohort preparation joins within cycle, applies the specified eligibility/outcome rules and retains the full MEC frame. The numerical module runs the benchmark, evaluates its current predictions and generates current outputs, with compatible-stage reuse under the existing cache contract. New results go to `outputs/`; bundled aggregate references remain under `reference/`.

This is the numerical reconstruction workflow, not a requirement for editing or inspecting supplied figures. `--prepare-only` is a distinct mode and cannot be combined with `--report-only`. The two historical record utilities are not part of the normal command sequence. Existing timestamps and null remote-registration fields describe the original chronology.

## Prediction schema and statistical conventions

The evaluation entry expects one row for each of the 7,477 eligible test IDs, with finite probabilities in [0, 1]:

```text
person_id,p_B,p_M,p_S,p_MEAN,p_P
```

`person_id` uses `cycle_SEQN`. The evaluator reads the outcome, weights and survey-design fields from the cohort/frame. Only the seven numeric and nine categorical X fields in [the feature dictionary](../protocol/feature_dictionary.json) are baseline predictors. Outcome, follow-up, mortality metadata, eligibility flags and survey weights are not additional prediction features. RDW and hemoglobin form the specified Z block.

Preprocessing and proxy selection stay within each training subset; grouped outer OOF construction repeats proxy selection inside each training partition. B, M and P are separately fitted routes selected from the specified ordered candidate family using the validation cycle. S and MEAN are alternative queries of the same fitted M, with preprocessing and coefficients fixed. MEAN uses development-weighted measurement means. Test data do not select models or provide calibration refits.

Single-cycle fits use WTMEC2YR; four-year development/test pools use WTMEC2YR/2. Each fit normalizes its current-fitting weights to mean one. Descriptive pooling across all four cycles uses WTMEC2YR/4. The convenience four-year `weight` column is not an eight-year national population weight.

Evaluation retains the full later-cycle MEC frame with zero contributions outside the analytic domain, paired Taylor covariance, pointwise 95% intervals and unpenalized Fieller sets. All returned set types, signs and finite/infinite endpoints remain meaningful. Primary and age 40–79 sensitivity analyses use their specified domains and common full test frame. These are fixed-prediction design approximations. Public follow-up boundaries and the MAXM-to-HMX transition constrain interpretation; they do not identify altered records or an instrument effect. Definitions are in [STATISTICS.md](../protocol/STATISTICS.md).

## Protocol and chronology

The analysis followed the [local V2 plan](../protocol/PREMODEL_ANALYSIS_PLAN.md). The [event record](../protocol/DEVIATIONS.md) records cancellation of the seven-day deadline at `2026-10-19T09:39:53+00:00`, followed by the author's decision to run the local analysis before repository publication at `2026-10-19T09:47:25+00:00`. The [protocol record](../protocol/PROTOCOL_RECORD.json) gives the first study-model fit as `2026-10-19T09:47:28+00:00` and `remote_recorded_utc` as null. No pre-model remote registration or formal preregistration was completed. Later publication does not change that ordering. Existing machine-readable records and scientific event bodies are retained.

## Source-only export limits

The current full collection includes `reference/` aggregates needed by saved-output commands. `python scripts/export_code_only.py --destination NEW_SIBLING_DIRECTORY` creates a smaller explicit source/protocol subset; it omits those references, manuscript integration materials, figures, raw inputs, participant records, predictions, models and caches. The destination must be new and outside the module. Reporting in the smaller tree requires the documented saved references supplied separately. Its manifest lists the actual exported paths.

The full collection's repository-level license and validation inventories are not automatically included in the smaller subset. Retain applicable notices when preparing a separate redistribution. Official source conditions are in [SOURCES.md](../sources/SOURCES.md); historical export provenance is in [README_original_source_export.md](../README_original_source_export.md).
