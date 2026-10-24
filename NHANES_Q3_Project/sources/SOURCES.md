# Primary sources and input provenance

Checked / downloaded for this package on 2026-10-19. All study inputs are NCHS public-use files. Source links identify the official provider; they do not imply a new license issued by this package.

## 1. Mortality linkage and outcome definition

**NCHS, National Death Index / linked mortality files.**
https://www.cdc.gov/nchs/linked-data/mortality-files/index.html

The official page distinguishes the 2019 public-use files from the 2022 restricted-use linkage. This protocol uses the public 2019 release. For selected public records, follow-up time or underlying cause may have been replaced by synthetic values; vital status was not changed. This limits interpretation of an exact 60-month endpoint without making the baseline CBC measurements constructed outcomes.

**NCHS, Codebook for the 2019 Public-Use Linked Mortality Files**, created 2022-04-28.
https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-files-data-dictionary.pdf

The codebook identifies SEQN for NHANES linkage, ELIGSTAT for linkage eligibility, MORTSTAT for final vital status and PERMTH_EXM for months from the MEC examination. It identifies DIABETES and HYPERTEN as multiple-cause-of-death flags; they must not enter the baseline history features.

**Official R fixed-width read-in example**, May 2022.
https://ftp.cdc.gov/pub/Health_statistics/NCHS/datalinkage/linked_mortality/R_ReadInProgramAllSurveys.R

The official example is intentionally omitted from this candidate; obtain it from the official URL above if needed independently. The preparation script uses the official NHANES positions: SEQN 1–6, ELIGSTAT 15, MORTSTAT 16, PERMTH_EXM 46–48, expressed as zero-based half-open positions for pandas. The URL is retained for provenance; Python preparation does not require R or that file.

## 2. CBC measurements and cycle choice

Official CBC documentation for the selected cycles:

- 2003–2004: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/L25_C.htm
- 2005–2006: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/CBC_D.htm
- 2007–2008: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/CBC_E.htm
- 2009–2010: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/CBC_F.htm

The two prespecified measurements are LBXRDW (red cell distribution width, %) and LBXHGB (hemoglobin, g/dL). Their actual fields exist in each of the four downloaded CBC files.

**2013–2014 CBC analytic notes**, used to make the cycle boundary decision:
https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/CBC_H.htm

The earlier draft considered the 2012 transition, but that was incomplete. The official 2007–2008 change note explicitly states MAXM in 2005–2006 and HMX in 2007–2008. The 2009–2010 laboratory manual also specifies HMX. Thus the retained 2003–2006 versus 2007–2010 split crosses MAXM→HMX. Generic MAXM boilerplate on the 2007/2009 webpages must not override the explicit change note and manuals. V2 retains the split and makes this a transport limitation; observed cycle distributions cannot identify an instrument effect or justify test-based harmonization.

Official method manuals:
- 2003: https://wwwn.cdc.gov/nchs/data/nhanes/public/2003/labmethods/l25_c_met_complete_blood_count.pdf
- 2005: https://wwwn.cdc.gov/nchs/data/nhanes/public/2005/labmethods/cbc_d_met.pdf
- 2007: https://wwwn.cdc.gov/nchs/data/nhanes/public/2007/labmethods/cbc_e_met.pdf
- 2009: https://wwwn.cdc.gov/nchs/data/nhanes/public/2009/labmethods/CBC_F_met_HE.pdf

See protocol/CBC_INSTRUMENT_AUDIT.md for the per-cycle evidence and the limits of the calibration claim.

## 3. Prior scientific motivation, independent of this performance analysis

Perlstein TS, Weuve J, Pfeffer MA, Beckman JA. **Red blood cell distribution width and mortality risk in a community-based prospective cohort.** Archives of Internal Medicine. 2009;169(6):588–594.
https://doi.org/10.1001/archinternmed.2009.55
https://pubmed.ncbi.nlm.nih.gov/19307522/
https://pmc.ncbi.nlm.nih.gov/articles/PMC3387573/

This NHANES III study provides prior evidence of an RDW–mortality association. It justifies considering RDW before inspecting the proposed benchmark performance. It does not establish incremental Brier improvement over the specified 16-feature baseline, or show that RDW can be replaced by a learned proxy. Hemoglobin is included as a second prespecified, interpretable measurement from the same routine CBC panel; the pair has not been selected by this task's observed prediction gain.

## 4. Survey uncertainty

**CDC, NHANES Variance Estimation tutorial.**
https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx

Supports using weights, public strata/PSU variables and Taylor linearization, retaining the full examination design frame for subdomain variance, and calculating domain degrees of freedom from represented PSUs and strata. Public design fields are masked variance units; the package does not claim access to actual identifiable geographic clusters.

**R survey package, summary statistics documentation.**
https://r-survey.r-forge.r-project.org/pkgdown/docs/reference/surveysummary.html

Reference for survey means, covariance and design degrees of freedom. The Python implementation in this package computes the specified fixed-prediction mean contrasts directly and does not require installing a second statistical runtime.

**Fieller EC. Some Problems in Interval Estimation.** Journal of the Royal Statistical Society, Series B. 1954;16(2):175–185.
https://doi.org/10.1111/j.2517-6161.1954.tb00159.x
https://academic.oup.com/jrsssb/article/16/2/175/7026708

Source of the ratio confidence-set approach. The protocol uses unpenalized Fieller sets with the joint design-based covariance and a t critical value; this application is an approximation conditional on fitted predictions. Unbounded sets are retained, not “repaired” with denominator penalties.

## 5. Software APIs

Use the installed supported API versions, recording them in the executed report:

- HistGradientBoostingClassifier: https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.HistGradientBoostingClassifier.html
- HistGradientBoostingRegressor: https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.HistGradientBoostingRegressor.html
- Ridge: https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.Ridge.html
- SplineTransformer: https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.SplineTransformer.html
- LogisticRegression: https://scikit-learn.org/stable/modules/generated/sklearn.linear_model.LogisticRegression.html

Explicit early_stopping=False and ordered finite candidate lists prevent software defaults from introducing additional unreported validation choices. These links are implementation references, not evidence of model performance in this study.

## 6. Data use and completeness

**NCHS Data User Agreement.**
https://www.cdc.gov/nchs/policy/data-user-agreement.html

Retain official public-use provenance and statistical-use conditions. Do not relabel the source files as newly generated CC0/MIT data. No restricted or identified records are included or needed.

`source_manifest.json` contains all 40 exact download addresses, the corresponding codebook links, expected downloaded byte lengths and package-relative locations. It is the executable input list used by `fetch_inputs.py` if a file needs recovery.

Actual preparation results are local analysis provenance, not claims taken from a published article:

- `data_checks/cohort_flow.json`: raw join validation and successive exclusions by cycle;
- `data_checks/feasibility_profile.json`: analytic sample sizes, events, X missingness, and survey design counts;
- `scripts/prepare_cohort.py`: exact code that produces these counts and the two delivered data tables.

The current candidate omits the 40 participant-level source files; acquire the exact specified releases using `scripts/fetch_inputs.py`. Codebook and provider links are retained. An earlier isolated full NHANES reconstruction and 24-candidate numerical run are recorded in [the comparison](../../validation/nhanes_rebuild_comparison.json) and [validation history](../../VALIDATION.md). They are completed prior evidence, not new training performed during this editorial preparation.


## 7. V2 weight, perturbation and protocol-record sources

- NCHS weighting tutorial, including WTMEC2YR divided by the number of combined post-2001 two-year cycles: https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx
- NCHS 2019 public-use LMF description, including synthetic follow-up and the unchanged vital status: https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-file-description.pdf
- Continuous NHANES public-use linked mortality data citation specified by NCHS: https://doi.org/10.15620/cdc:117142
- Public-use LMF data dictionary: https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-files-data-dictionary.pdf — it has no record-level perturbation flag or unperturbed follow-up time. Boundary counts do not identify affected records.
- Git commit documentation (author/committer dates can be set): https://git-scm.com/docs/git-commit
- OSF registration overview and secondary-data registration option: https://help.osf.io/article/330-welcome-to-registrations — an optional alternative record, not one created by this package.

V2 considered but did not add albumin/creatinine. The extra processing context includes the 2005–2006 creatinine analytic notes and the 2007–2008 LX20/DxC800 chemistry instrument change:
- https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/BIOPRO_D.htm
- https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/BIOPRO_E.htm

V2 local descriptive evidence is reproducible from scripts/pre_model_descriptives.py. It produced the cycle distribution table/ECDF and public follow-up boundary counts without fitting study models. These analyses were inspected before recording the V2 plan, and are disclosed as such.
