# Predicting a measurement is not replacing it: a matched evaluation of learned measurement proxies

This repository accompanies the study titled “Predicting a measurement is not replacing it: a matched evaluation of learned measurement proxies”. It evaluates measurement prediction, downstream increment, and fixed-model retention through matched comparisons in liver-safety and clinical-risk tasks. It provides NHANES saved-summary checks, saved Figure 6/S12 reconstruction and the NHANES numerical pipeline. The minimal liver-safety input collection now contains all 22 inputs required by the existing fixed-prediction replay, including the original scores, labels, sampling arrays and expected results. It does not provide a complete liver-safety training environment. Full historical liver-safety model training remains NOT_EVALUATED. Earlier isolated NHANES reconstruction completed the full 24-candidate workflow; current publication preparation reuses that evidence without retraining. See [validation scope](VALIDATION.md) and [current release status](RELEASE_STATUS.json).

Author: Yuan Cao, Independent Researcher, Beijing, China.

**Publicly available on GitHub.** The repository is available at https://github.com/TeryYuan/learned-measurement-proxy-evaluation. Original software is licensed under MIT; original documentation, clinical aggregates and artwork use CC BY 4.0 where listed. Source-related files retain their specific conditions, including the listed NC and SA author grants. See [LICENSE_SCOPE.md](LICENSE_SCOPE.md), [FILE_LICENSES.json](FILE_LICENSES.json) and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). The repository is public; no GitHub Release or Zenodo archive is asserted. Their separate current states are recorded in [RELEASE_STATUS.json](RELEASE_STATUS.json).

## Quick start

From the repository root, use CPython 3.14.3 and the tested compatible lock:

```sh
python -m pip install -r requirements.lock.txt
python tools/verify_nhanes_summary.py --input NHANES_Q3_Project/reference/outputs/tables/fixed_substitution_summary.json --output summary_check.json
python NHANES_Q3_Project/scripts/render_saved_figures.py --report-only --report-dir saved_figures
python NHANES_Q3_Project/scripts/run_module.py --report-only --report-dir saved_report
python tools/run_fixed_replay.py --bundle release_candidate/legacy_fixed_replay --output-dir liver_replay_results --timeout-seconds 1200
```

Choose new writable output locations. The summary entry checks 385 identities. Drawing reconstructs Figure 6 and S12 in PDF/SVG/PNG; report-only generates saved tables and Figure 6 without fitting. S12 uses the preserved curve geometry, not person-level ECDF recomputation. Liver replay checks 238 existing numerical comparisons at absolute tolerance `2e-12`, with zero new random draws or model queries; it may take several minutes. The wrapper excludes old output and runs from a fresh working directory. It is not itself a filesystem sandbox.

The bundled minimal liver replay needs no external download or private model cache. Five tables are column projections; all retained values and row order, all six original index files and all nine expected tables are unchanged. Original C/O grouping relations are preserved in an additional membership table. Omitted chemical strings and model-path annotations are documented in [DATA_ACCESS_AND_REPLAY.md](DATA_ACCESS_AND_REPLAY.md). Removing those fields does not remove the conditions on retained source labels or measurements.

## NHANES raw-data reconstruction

Use a separate writable checkout. Obtain the 40 official files specified in the module manifest and run the original numerical workflow:

```sh
cd NHANES_Q3_Project
python scripts/fetch_inputs.py
python scripts/prepare_cohort.py
python scripts/run_module.py
```

Acquisition needs network access. Raw records are not bundled. The full pipeline writes new predictions and evaluates them under `outputs/`; it does not reuse static statistical references merely because a summary exists. Preserved aggregates remain under `reference/`. See [module instructions](NHANES_Q3_Project/README.md) and [official sources](NHANES_Q3_Project/sources/SOURCES.md).

## Other entry points and scope

```sh
python -m unittest discover -s tests
python tools/fetch_liver_sources.py --sources SI002 SI004 --output-dir liver_sources
```

The optional source downloader validates exact historical workbooks; it does not recover historical predictions or replace the saved draws. FDA acquisition currently requires an exact legally obtained local copy when its historical URL is unavailable. TG-GATEs and Hallmark access for historical training are described separately; no complete automated training reconstruction is claimed.

Read [reproduction scope](REPRODUCTION_SCOPE.md), [current versus prior validation](VALIDATION.md), [source decisions](sources/REPLAY_INPUT_DECISIONS.csv), and [data access](DATA_ACCESS_AND_REPLAY.md). Cite this work using [CITATION.cff](CITATION.cff) and acknowledge the applicable upstream sources. The public repository URL is recorded in CITATION.cff. Release version, DOI and archival publication date remain unset until the corresponding records exist. The full manuscript, Supplementary Information and private author-review records are distributed separately from this repository.

## Author information

Author: Yuan Cao (曹源), Independent Researcher, Beijing, China. Correspondence: 18281658226@163.com. This research received no funding.

## Declaration of interests

The author is currently enrolled at the Chinese Academy of Medical Sciences & Peking Union Medical College, Beijing, China. This work was conducted independently of the author’s degree program, using personal resources and public data, without institutional support. The author declares no competing interests.

## Scientific result attachments

The [reader index](source_index.tsv) maps every distributed file to its evidence role, manuscript location, delivery and rights status. [Scientific attachments](scientific_SI/README.md) add 140 saved table files, seven chemical-coverage and three scaffold-stratified files; [naming_map.csv](naming_map.csv) retains original model labels. These additions were checked against original saved files, not generated by new training or bootstrap. Fourteen table files retain unresolved grouping cells; extended group-dependent replay remains limited. The historical 438-row author index is not the current public inventory.

See [release-preparation checks](validation/release_preparation.json), [prior scientific validation](VALIDATION.md) and [license scope](LICENSE_SCOPE.md). The exact delivered file set is listed in [PUBLIC_EXPORT_INVENTORY.json](PUBLIC_EXPORT_INVENTORY.json); [SHA256SUMS](SHA256SUMS) covers all files except itself. Private author material and development Git history are excluded. The GitHub repository is public; archive identifiers remain unset until real records exist.

## Current publication figures

The current six-figure presentation is in [figure_patch/main_figures](figure_patch/main_figures/), with its saved-result renderer, fonts instructions and per-file rights in [the patch README](figure_patch/README.md). The Figure 6 files under figures/ match that presentation exactly; Figure S12 is retained unchanged. Earlier module plotting commands still reconstruct the historical saved-result layout. This distribution includes the current 174 mm figure files; the module commands reconstruct the historical 167 mm Figure 6/S12 layout. Source SVG icon credits and retained license conditions accompany the patch.
