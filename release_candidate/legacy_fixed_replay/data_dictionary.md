# Fixed-results data dictionary

This dictionary describes existing saved records; it does not redefine labels, members or comparisons. `field_inventory.csv` enumerates every CSV column, stored dtype and missing-value count. `provenance.csv` maps each local file to its original project file. Full-precision numerical strings are read with pandas `float_precision="round_trip"`.

## Units, keys and score meanings

| File family | Grain and key | Outcome | Score |
|---|---|---|---|
| N_predictions | member_id × model; 125 members, 12 names | y_v2: V2 Most/Less=1, No=0 | raw_score: saved fixed full-O probability; calibrated_score: saved O-trained calibration; no N calibration |
| BSEP_predictions | member_id × model; 201 identities, two forests | y_BSEPi: observed author SI004 BSEP binary class | raw_score: fixed corresponding-assay forest probability; no DILI score |
| C_OOF | member_id × model × training_version; 624 × 4 × 2 | y_v1/y_v2 are separately retained; training_version chooses fitted label, scoring_version chooses evaluated label | raw_score: held-out C-fold probability; not all-O training prediction |
| NAM_predictions | condition_id × route_model; 72 conditions × 11 names, nested in 24 parent_id | y: mean(log2 treated-animal ALT) minus mean(log2 matched-control ALT) | prediction: fixed condition-level ALT log2-response; unconstrained continuous output |
| bridge_membership | identity_id; E_BSEP intersect C, 136 identities | y_BSEPi, DILI_y_v1, DILI_y_v2 are distinct endpoints/versions | BSEP_score_real/perm: saved assay scores; C_M0/C_M1 fields: old-label-trained C OOF scores |

The original enriched tables used chemical `identity_group`, `fold_group` and scaffold strings. Those representations and other excluded annotation columns are not distributed in the minimal projections. `group_membership.csv` preserves the exact original equivalence relations with an anchor from existing table/member IDs; it does not create new memberships. Original member IDs, LTKBIDs, C numerical folds and source row references remain where applicable. See the root DATA_ACCESS_AND_REPLAY.md and sources/replay_projection_columns.json for exact omissions. Some historical membership flags predate subsequent authorizations; the saved evaluated cohort is given by the corresponding file/ID arrays, not by treating an old flag as a new exclusion rule.

`raw_threshold` is the exact saved development-selected threshold; `predicted_class` applies it to raw_score. BSEP's `fixed_descriptive_threshold` is 0.5 and has a different role. `backbone_raw`, `logit_correction`, `probability_correction` and `q_mean` describe the saved kernel path; missing correction fields for other model families mean not applicable, not zero. E_BSEP is the main training-identity-cleared assay set; E_no_anyA is its nested 172-member sensitivity. O/C/N/A overlap columns are annotations, not additional independent validation populations. Empty failure/NA-reason fields indicate no reason was recorded; missing numerical measurements are never negative labels.

NAM dose_mg_kg, log10_rat_mg_kg and log10_PHH_uM retain source-system units; no molar-to-mg/kg conversion is implied. `error=prediction-y`, `absolute_error=abs(error)`, and `squared_error=error²`. The original model_path column was historical provenance only; it is omitted from the minimal public projection and is not needed to run this replay. Six names (D_cond, D_struct, D_meas_only, D_meas, D_swap, D_meanZ) are exactly the same final prediction vector. D_perm is distinct at full precision. M0/B_star are aliases in N.

## Original index matrices

All arrays are zero-based and must be applied in their explicitly saved identity order. NPZ files are loaded with `allow_pickle=False`.

- N_bootstrap_indices.npz contains ids, labels and indices; rows of indices are the original 2,000 shared class-stratified draws.
- E_BSEP_bootstrap_indices.npz and E_no_anyA_bootstrap_indices.npz have the same roles within their own fixed cohorts. Do not use one cohort's array as an independent replicate of the other.
- bridge_indices_BSEP/DILI.npz contain identity_id, stratification_label and indices. BSEP- and DILI-stratified arrays serve the different corresponding outcomes; the saved BSEP array also supports the original proportion-difference calculation.
- NAM_indices.csv has bootstrap_draw plus draw_position_0…23. Positions refer to parent_ID_order in NAM_context.json. Every selected parent's conditions stay together; repeated parent selections retain multiplicity.

No new RNG is used. Indices are never reconstructed from seed as part of this replay.

## Estimators and expected-table fields

AP is average_precision_score, not trapezoidal PR area. ROC_AUC is area under the ROC curve. Raw/calibrated Brier is mean squared probability error against the relevant binary outcome. N primary differences are M1−M0 and R_gate−B_star. BSEP real−perm and bridge paired_difference fields are differences of two metrics on shared identities, not single-model intervals. C training/scoring versions form four cells; only `scope=pooled_OOF` is replayed here.

NAM MAE/RMSE/R2 use equal parent weights after averaging conditions within each parent. The replay uses saved parent sampling indices and recomputes weighted target variance in each existing draw, with the recorded near-zero variance guard. `model_a−model_b` fixes contrast direction; a negative MAE difference favors the first. `confidence`/`coverage`, low/high (or ci_low/ci_high), and estimate/point retain the historical interval type. No P value is derived from bootstrap signs or from overlap of single-model intervals.

Bridge association rows a,b,c,d correspond to [[BSEP1,DILI1; BSEP1,DILI0],[BSEP0,DILI1; BSEP0,DILI0]]: [[25,1],[100,10]]. Conditional OR and its interval use scipy.stats.contingency.odds_ratio(kind="conditional"). The proportion-difference interval uses the saved BSEP-stratified draws, not the conditional-OR procedure. These remain post-hoc nominal descriptions.

Expected CSVs are comparison targets, not inputs to new fitting. `output/replay_checks.csv` records equality to the historical points/endpoints at absolute tolerance 2e-12. Existing supporting tables not covered by the 238 listed checks remain copied saved evidence, not newly recomputed claims. Full historical training is not reproduced by this lightweight layer.
