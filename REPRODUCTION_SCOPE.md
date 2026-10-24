# Reproduction scope

| Layer | Candidate capability | Evidence |
|---|---|---|
| NHANES saved-summary identities | Included aggregate reference, 385-check entry | Passed in preceding QA; the same 385-check entry is included for local saved-summary verification |
| NHANES saved figures and reports | Figure 6/S12 rebuild and report-only without training | Passed in preceding QA; saved reference estimates are unchanged. These commands rebuild the historical module layout. |
| NHANES raw/cohort/training pipeline | Network acquisition of 40 official files, original pipeline in a writable workspace | Preceding isolated reconstruction/training passed; not newly trained here |
| Minimal liver fixed-prediction replay | All 22 required inputs supplied in the collection; independent command in README | Prior execution result is recorded in VALIDATION.md and validation/publication_preparation.json; no raw-data/model access is needed |
| Historical enriched liver annotations | Chemical strings and model paths omitted from five projections | Field-level remaining issues documented; group equivalence, original member IDs, folds and row order retained |
| Liver raw-data/full training | Historical source code is inspectable; raw inputs, complete caches/environment/model artifacts not supplied | NOT_EVALUATED; fixed-score replay does not establish full training |
| Filesystem isolation | Prior isolated NHANES evidence retained; prior minimal liver replay tested in a separate namespace root | Only the actually executed layer may be called isolated; the ordinary Python wrapper alone is not a sandbox |

The public file scope and original grants are specified in FILE_LICENSES.json. Publication state and real endpoints are recorded in RELEASE_STATUS.json. Scientific content, ethical facts, model candidates, label definitions, seeds, expected results and statistical tolerances remain unchanged; current authorship and availability prose are updated. No new model was trained or queried.

The expanded scientific attachment collection is a saved-result distribution with original-token verification. It is not covered as a whole by the 238-check replay. See [attachment limits](scientific_SI/README.md) and [unresolved grouping](scientific_SI/LOCAL_RECOVERY_REQUIREMENTS.json). The prior isolation result is inherited only for unchanged executable paths and original replay inputs; no fresh isolated reconstruction of every attachment is claimed.

## Current figure presentation

The six current 174 mm figures and saved-result drawing source are in [figure_patch/](figure_patch/README.md). Figure 2 now labels exact prevalence positions 37/201 and 34/172 as 0.184 and 0.198. Figure 5 moves the inset gamma label inside its axes and adds 12 pt of bottom margin. Their numerical points, interval endpoints and MCSE values are unchanged. Figures 1/3/4/6 retain their previous independent files. The module saved-figure commands continue to produce the historical 167 mm Figure 6/S12 layout; they do not recreate the newer six-figure presentation. Current presentation checks do not replace or extend the earlier numerical-validation scope.
