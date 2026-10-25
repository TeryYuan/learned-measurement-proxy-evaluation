# Data access and fixed-prediction replay

## Public distribution with complete fixed-prediction replay inputs

The public repository at https://github.com/TeryYuan/learned-measurement-proxy-evaluation supplies all 22 files used by the unchanged `release_candidate/legacy_fixed_replay/replay.py`. The [22-file decisions](sources/REPLAY_INPUT_DECISIONS.csv) retain the original file-level inclusion decisions. Current author grants and retained source conditions are documented in [LICENSE_STATUS.md](LICENSE_STATUS.md). Read [VALIDATION.md](VALIDATION.md) for actual current execution evidence and its limits.

From the repository root:

```sh
python -m pip install -r requirements.lock.txt
python tools/run_fixed_replay.py --bundle release_candidate/legacy_fixed_replay --output-dir liver_replay_results --timeout-seconds 1200
```

No external data, private cache or model is needed for this fixed-score command. Its output must report 238/238, absolute tolerance `2e-12`, zero new draws and zero new model queries. Existing `output` directories are excluded by the wrapper. Saved inputs are not produced by rerunning a seed, downloading a new database or inferring labels from target metrics.

## Exactly what changed in the input packaging

All six original index files (five NPZ and one CSV), both JSON contexts and nine expected-result tables remain byte-for-byte unchanged. NPZ arrays were inspected with `allow_pickle=False`; they contain original member IDs, labels and sampling positions, not positions alone. SI004 labels in three NPZ files retain their noncommercial source conditions. FDA labels in the other two remain attributed to FDA. NAM positions refer to the unchanged original parent order.

Five CSVs retain the original row order and every retained string token, including full-precision numbers. The [column record](sources/replay_projection_columns.json) lists exact included and excluded columns. This is not a new label or cohort definition. No original research file was overwritten.

| File | Included information | Excluded annotation scope |
|---|---|---|
| N_predictions.csv | All saved scores, thresholds, classes, correction values, failure fields, V2 labels/categories and member/model IDs | Chemical identity string and source-name annotations |
| BSEP_predictions.csv | Scores, labels, original SI004 row links, overlap/membership and eligibility flags, 201-member evaluated set | Chemical representations, chemical grouping strings, name/CAS and source identity strings |
| C_OOF.csv | Original 624 member IDs, LTKBIDs, label versions/categories, original folds and cohort flags, all model/training-version scores | Chemical strings, scaffold strings and source-name annotations |
| NAM_predictions.csv | All 72 conditions, 24 parents, 11 routes, original doses/source compound IDs/names, targets and predictions | Unused historical model-path column only; it cannot be resolved by public readers and is not needed by replay |
| bridge_membership.csv | Original 136 IDs, C/SI row mapping, FDA label versions, C outer fold, overlap flags, observed labels and scores | Chemical representations/group strings, name/CAS/source identity strings and parent SMILES |

The [group membership table](release_candidate/legacy_fixed_replay/group_membership.csv) retains equality of every original identity, fold and scaffold grouping across the four relevant tables. Each group is represented by the lexicographically first original `table:member_id` already belonging to that group. Original member IDs are never changed. This does not reassign folds or run a new grouping algorithm; numerical C fold assignments remain in C_OOF.csv. The original chemical-string crosswalk is retained privately. [Field-level excluded-material decisions](sources/REPLAY_EXCLUDED_FIELDS.csv) distinguish unresolved chemical redistribution from metadata-only model paths. Deleting annotations does not erase the upstream conditions on retained labels/measurements.

The original E172 array is a nested subset of the original E201 membership. ALT parent order and all paired positions are unchanged. Six ALT route names remain exactly equal; D_perm retains full precision. CSV consumers must use `float_precision="round_trip"` where applicable. No scientific formula, expected value or tolerance was adjusted.

## Optional acquisition of original source workbooks

The script below is a working source-acquisition/audit utility, not an assembler of missing historical predictions:

```sh
python tools/fetch_liver_sources.py --sources SI002 SI004 --output-dir liver_sources
python tools/fetch_liver_sources.py --input-dir legally_obtained_sources --output-dir validated_sources
```

The second command requires these exact files in the selected input directory: `tx4c00453_si_002.xlsx`, `tx4c00453_si_004.xlsx`, `dilirank_fda_current.xlsx`. It uses no network. [The pinned manifest](sources/liver_source_downloads.json) gives official URLs, sizes, exact historical version digests and filenames. It refuses HTML error pages, a changed current version, invalid XLSX structure, missing FDA V1/V2 sheets, missing files and an existing destination. It never reads the author's original project by default.

SI002 and SI004 version 1 were freshly downloaded from the official file endpoints in this stage and matched the preserved historical bytes. Their CC BY-NC 4.0 terms apply. The FDA historical download URL was unavailable in this stage; the legally held exact historical workbook was validated through the local-input route. A future downloader must obtain that exact version or stop; a current workbook with changed bytes is not substituted. Crucially, historical FDA V1 is its own worksheet, not V2 with new compounds removed. These original workbooks are not inside the public ZIP and are not required by the bundled minimal replay.

## Historical training and other sources

For NHANES, use the module's official 40-file manifest and the three commands in [README.md](README.md). The Python workflow does not require the omitted official R reader.

Open TG-GATEs source acquisition uses the [official catalog](https://dbarchive.biosciencedbc.jp/en/open-tggates/download.html); the study archived observations on 2026-10-08. A mutable LATEST URL is not a frozen experiment. Historical cohort construction also needs original condition selection, matched controls, expression processing, source crosswalks, training membership and model artifacts. Hallmark 2025.1.Hs must be obtained through an authorized MSigDB/declared distribution route with its terms; do not silently upgrade to current sets. Details and links are in [DATA_SOURCES.md](DATA_SOURCES.md).

A general raw-data-to-all-models liver training reconstruction is **not provided or validated** here. No placeholder `assemble_liver_replay_inputs.py` is supplied: the 22 actual replay inputs are already present, and downloading raw data cannot recreate the historical scores or saved arrays. Omitted enriched annotations affect extended provenance/structure inspection and historical retraining, not the existing 238-check entry. Approval of any broader annotation/raw/model distribution needs a new explicit manifest and the relevant source facts.
