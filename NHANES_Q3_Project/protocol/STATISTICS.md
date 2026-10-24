# Fixed-prediction survey comparison and Fieller confidence sets

This document specifies the new clinical module only. Existing ALT MAE, simulation, and DILI AP intervals retain their original definitions.

## Estimand

For route r, the primary loss is the MEC-weighted mean Brier loss in the complete, linked analytic domain:

\[
\widehat L_r=\frac{\sum_i w_i I_i(Y_i-p_{ri})^2}{\sum_i w_i I_i}.
\]

I denotes analytic eligibility. All route comparisons use the same I, Y, and w. Keep the full final-test MEC sampling frame with valid positive weights and design fields; participants outside the analytic domain contribute zero to linearized totals. They need no model predictions.

The final four-year test weight is WTMEC2YR/2, and the four-year development weight is separately WTMEC2YR/2. A single two-year training/validation cycle or single-cycle description uses WTMEC2YR; an all-four-cycle description uses WTMEC2YR/4. Derive weights from WTMEC2YR using the number of prespecified cycles in the relevant pool. The delivered CSV column `weight` is a four-year-period convenience column, not a universal eight-year population weight.

Every proxy and downstream fit, including nested CV and OOF, uses current-pool weights normalized to the mean of its current fitting subset. Proxy selection and downstream validation are weighted as specified in config. Ordinary training medians, scaler moments and spline knots remain the fixed preprocessing conventions; they do not turn fitting losses into unweighted losses.

This is an analytic-cohort benchmark. MEC weighting and full-frame domain variance do not automatically correct selection through CBC missingness or linkage ineligibility, or establish national representativeness. The evaluation crosses the documented MAXM-to-HMX transition, and conditions on the public-use follow-up times.

## Joint Taylor covariance

Let v_i be the vector of the five Brier losses, W=sum(w_i I_i), and theta the vector of their weighted means. For stratum h and PSU j:

\[
u_{hj}=\frac{1}{W}\sum_{i\in(h,j)}w_i I_i(v_i-\widehat\theta).
\]

With m_h public PSUs in the full design frame, and no supplied finite-population correction:

\[
\widehat V(\widehat\theta)=\sum_h\frac{m_h}{m_h-1}
\sum_j(u_{hj}-\bar u_h)(u_{hj}-\bar u_h)^\top.
\]

Use the actual PSU count; at least one selected NHANES cycle includes a stratum with three PSUs. IDs nest cycle, stratum and PSU. Do not treat SDMVPSU=1/2 as two clusters for the entire dataset.

The t degrees of freedom are the number of PSUs containing domain observations minus the number of strata containing domain observations. Record both the full-frame and active-domain counts. With no domain-empty PSU these coincide with the usual full-frame difference. Do not automatically invent a singleton correction if a contributing full-frame stratum has only one PSU; diagnose the input/definition.

For any coefficient vector k, contrast estimate=k' theta and variance=k' V k. Use estimate +/- t(df,0.975)*SE for a pointwise 95% interval. The covariance and df apply to the fixed selected predictions. They do not incorporate re-fitting or public-file linkage/measurement error.

The reference implementation is `scripts/evaluate_fixed_predictions.py`. Its finite hand-calculated tests are `scripts/check_statistics.py`. It does not query or refit models. A single optional comparison with R `survey::svymean` may be used if R is already available; it is not a reason to require a second runtime or perform repeated resampling.

## Contrasts and signs

Route order is B, M, S, MEAN, P:

| Estimand | Coefficients | Interpretation of positive values |
| --- | --- | --- |
| G = L_B-L_M | 1,-1,0,0,0 | measured inputs improve Brier |
| D = L_S-L_M | 0,-1,1,0,0 | fixed substitution worsens Brier |
| H = L_B-L_S | 1,0,-1,0,0 | substituted route retains gain over B |
| A = L_MEAN-L_S | 0,0,-1,1,0 | proxy beats a training-mean substitution |
| G_P = L_B-L_P | 1,0,0,0,-1 | retrained proxy route improves over B |

The R_fixed numerator is H, not G_P. The identity H=G-D must hold numerically. Replacing a measurement with its training mean need not reproduce B, because B is a separately fitted function; therefore MEAN does not automatically define 0% retention.

## Fieller confidence set

For G>epsilon (epsilon=1e-8 Brier units), R_fixed=H/G. This is a ratio of two survey-weighted mean differences, not a mean of person-level ratios. It differs in aggregation from the old simulation's mean across repeat-level ratios; both definitions must remain explicit.

From the joint covariance obtain V_G, V_H and C_GH. Set q=t(df,0.975)^2. The unpenalized Fieller confidence set consists of all r satisfying:

\[
(G^2-qV_G)r^2-2(GH-qC_{GH})r+(H^2-qV_H)\leq0.
\]

Report bounded, one-sided, disjoint, or all-real sets as returned. This is a design-based normal/t approximation, not exact finite-sample coverage. Do not penalize the denominator to force boundedness, discard nonpositive bootstrap denominators, or clamp the interval to [0,1].

If G<=epsilon, the point retention is NA (nonpositive and small-positive reasons separate). A diagnostic Fieller set may still be retained, but it is not an interpretable useful-measurement retention interval. If G>epsilon while its interval crosses zero, preserve a descriptive point and the full confidence set, label them unstable, and do not headline a retention percentage. A positive G interval and a bounded Fieller set permit a quantitative retention statement, with its interval.

Pointwise intervals do not claim simultaneous family-wise coverage. No clinical noninferiority margin is supplied. A substitution-loss interval spanning zero, or a retention interval containing one, does not establish equivalence or authorize dropping a laboratory test.

## Fixed sensitivities

1. Unweighted Brier and the same signed contrasts on exactly the same predictions: descriptive, no new fits.
2. The 40-79-year final-test subdomain: same predictions, same full test frame, domain weights and Taylor covariance recalculated. Report domain n and events.

## Primary sources

- CDC NHANES Variance Estimation tutorial: https://wwwn.cdc.gov/nchs/nhanes/tutorials/varianceestimation.aspx
- R survey package, survey summary statistics and design df: https://r-survey.r-forge.r-project.org/pkgdown/docs/reference/surveysummary.html
- Fieller EC. Some Problems in Interval Estimation. JRSS B. 1954;16(2):175-185. https://doi.org/10.1111/j.2517-6161.1954.tb00159.x



Additional official weight and mortality-file sources are listed in LMF_WEIGHT_AUDIT.md; the deterministic-proxy Brier derivation is in INTERPRETATION_zh.md.
