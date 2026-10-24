# Current release preparation

The original rightsholder is confirmed as Yuan Cao. Current grants and publication state are defined in FILE_LICENSES.json and RELEASE_STATUS.json. This release was assembled from an explicit file whitelist. Scientific tables and minimal replay inputs retain their reviewed bytes; privacy edits are limited to specified documentation and figure metadata. Figure metadata changes preserve rendered PDF pixels and all non-date SVG content. No scientific model, new bootstrap, new statistical comparison or full isolated training run is counted in this stage.

See validation/release_preparation.json for these bounded checks. The following dated records retain the facts of their own stages, including former pending-license states. They are not the current publication status. Irrelevant author-handoff and document-export debugging details were removed from selected public summaries; original private records are retained by the author.

---

# Handoff review, 2026-10-23

The later handoff review verified the 315-file incoming candidate against its standalone archive, embedded archive and reader index. The 140 saved scientific tables, seven coverage files and three scaffold-stratified files agree byte-for-byte with the supplied recovery inputs. Independent formulas reproduced 24 saved S41 AP/AUC point estimates, with maximum absolute difference 1.11e-16; S40 summary statistics agree with the saved 24/201/172 nearest-neighbor values. These checks did not reconstruct molecular fingerprints, refit models or resample confidence intervals. Fourteen table files still contain unresolved grouping relationships.

The column-treatment pointer now names the delivered RECOVERY_RECORD.json. The final local candidate has its own relative-path SHA256SUMS, excluding itself. At that handoff stage, this manifest verified file integrity and license/publication decisions remained pending. See [handoff review](validation/handoff_review.json). The earlier records below retain the scope and results of their own stages. The separately supplied author handoff now includes updated manuscript PDFs; the preceding stage's PDF-export failure remains a historical fact.

---

# Earlier independent-researcher preparation

This stage updates authorship, declarations, reader navigation and saved scientific attachments. [Current checks](validation/ir_preparation.json) concern original-byte/token preservation, metadata, references and privacy. No model fitting, model queries, new bootstrap or full numerical audit was performed. The 238 minimal replay checks, 385 NHANES summary checks, earlier code tests and earlier isolated NHANES training retain their own scopes. No new filesystem isolation run was performed for the expanded attachment collection.

The following dated account is preserved as prior-stage evidence; “new” in that account means its own earlier preparation stage.

---

# Validation: publication preparation, 2026-10-23

The original accepted 116-file candidate is preserved as the comparison baseline. It was independently copied before this work. This stage prepares licenses, data scope and publication material; it does not change model candidates, seeds, expected values, formulas, tolerances, paper conclusions or ethical statements.

## New execution in this stage

| Check | Actual result and scope |
|---|---|
| 22 liver inputs | PASS: 17 files byte-unchanged; five exact-token column projections with original row order and all retained round-trip values equal. Group equivalence preserved; original member IDs/folds, E172 within E201, ALT parent order and six-route exact identity retained |
| Minimal proposed public liver replay | PASS: 238/238; original absolute tolerance 2e-12; no new draws, no new model queries; unchanged replay.py and expected tables; fresh output from only candidate inputs |
| New acquisition utility | 7 targeted boundary tests passed; SI002/SI004 v1 freshly downloaded and exact historical bytes verified; FDA historical endpoint unavailable, so its exact legally held workbook was validated through the explicit local-input route |
| Filesystem isolation | PASS for the new minimal liver replay and new tests: WSL unprivileged Linux user/mount/network/PID namespaces plus chroot, capabilities dropped, read-only candidate/runtime/libraries, fresh temporary/results directories and unrelated working directory. Original project, host home/mounts and author-review inputs absent; no external data mount; network unreachable |
| Full historical liver training | NOT_EVALUATED; no training performed. Exact replay is saved-score evidence, not raw-data/model reconstruction |

Detailed safe records: [replay](validation/minimal_liver_replay.json), [238 checks](validation/minimal_liver_checks.csv), [isolation](validation/minimal_liver_isolation.json), [input contracts](validation/input_organization.json), and [publication preparation](validation/publication_preparation.json). The Python replay wrapper does not itself supply filesystem isolation. Technical completeness of the local proposed collection is not authorization to publish it.

## Earlier evidence carried forward, not rerun here

The preceding code QA passed 34 tests on Linux, 385 NHANES saved-summary checks, Figure 6/S12 reconstruction and report-only entries. Windows had one platform symlink fixture skipped; that fixture passed on Linux. The preceding NHANES rebuild acquired 40 official files, reproduced the 13,236-row cohort and 39,608-row design frame, and ran the original 24-candidate training protocol in isolation. Selected configurations/folds matched; maximum final prediction difference was 8.89e-16 and maximum summary difference 8.11e-15. No new NHANES training or figure rendering is counted in this publication-preparation stage.

The prior private full-input replay also passed 238/238; its scope remains private. The new minimal-candidate replay above is distinct evidence from the actual projected collection. Historical records are retained in [prior QA status](validation/prior_qa_run_status.json); their earlier public-replay BLOCKED entry describes the then-withheld inputs, not today's proposed collection. Other older records in validation retain their original scopes and dates.

The original submission ZIP and its original 96-file candidate remain preserved. The accepted 116-file baseline was the starting point for this stage. Fourteen NHANES reference files and six existing figure assets remain byte-for-byte unchanged in the updated candidate; the updated package as a whole is not identical to either original baseline. Original manuscript/SI/legend and ethics material was not edited. Existing PDF/PNG format and visual checks are carried forward only for identical bytes; newly included NPZ files are parsed, not skipped as opaque binaries.

Final syntax, official CFF schema, relative-reference, file-inventory, privacy and ZIP results are recorded in validation/publication_preparation.json and validation/release_review.json. Sensitive match text, private paths and tracebacks are not included in these public records. Automatic hits must be resolved by meaning, and an unsupported binary blocks clearance. These bounded checks are not a guarantee against arbitrary hidden information.

No GitHub push, Release, deposit, Zenodo upload, email or manuscript submission was performed. MIT, original-content licenses, exact release identity and the final approved distribution scope remain author decisions.
