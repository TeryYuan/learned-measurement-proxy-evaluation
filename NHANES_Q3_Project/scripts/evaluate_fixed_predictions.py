"""Evaluate the five predeclared routes on identical later-cycle NHANES records.

This helper never trains, selects, calibrates, or queries a model. Its input is
the completed prediction table with person_id,p_B,p_M,p_S,p_MEAN,p_P. The
covariance is a stratified-PSU Taylor linearization of domain-weighted means,
conditional on the supplied predictions. No resampling or hashes are used.
See protocol/STATISTICS.md for the estimand and Fieller confidence-set rules.
"""
from pathlib import Path
import argparse
import json
import math

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

ROOT = Path(__file__).resolve().parents[1]
ROUTES = ["B", "M", "S", "MEAN", "P"]
CONTRASTS = {
    "G_measured_gain": [1, -1, 0, 0, 0],
    "D_substitution_loss": [0, -1, 1, 0, 0],
    "H_remaining_gain": [1, 0, -1, 0, 0],
    "A_proxy_over_mean": [0, 0, -1, 1, 0],
    "G_P_retrained_proxy_gain": [1, 0, 0, 0, -1],
}


def survey_mean_cov(frame, domain_values):
    """Return joint domain means, first-stage design covariance, and domain df.

    frame: one row per person, with person_id,weight,stratum_id,psu_id;
           includes the whole valid MEC frame, including domain exclusions.
    domain_values: DataFrame indexed by a subset of frame.person_id; all its
                   numeric columns are observed in every domain record.
    Domain exclusions contribute exactly zero to the linearized PSU totals.
    FPC is not supplied. The actual number of PSUs per stratum is used.
    """
    required = ["person_id", "weight", "stratum_id", "psu_id"]
    if frame["person_id"].duplicated().any() or domain_values.index.duplicated().any():
        raise ValueError("Duplicate participant IDs are not an independent sample")
    f = frame[required].copy()
    f["person_id"] = f["person_id"].astype(str)
    f = f.set_index("person_id")
    values = domain_values.copy()
    values.index = values.index.astype(str)
    if not values.index.isin(f.index).all():
        raise ValueError("Domain contains IDs absent from the survey design frame")
    if not np.isfinite(values.to_numpy(dtype=float)).all():
        raise ValueError("Domain outcomes/losses must be finite; do not silently drop a route")
    if not np.isfinite(f["weight"]).all() or (f["weight"] <= 0).any():
        raise ValueError("Positive finite MEC weights are required")
    w = f.loc[values.index, "weight"].to_numpy(dtype=float)
    total = w.sum()
    if len(w) == 0 or total <= 0:
        raise ValueError("Empty analytic domain")
    v = values.to_numpy(dtype=float)
    estimate = (w[:, None] * v).sum(axis=0) / total
    lin = pd.DataFrame(0.0, index=f.index, columns=values.columns)
    lin.loc[values.index, :] = w[:, None] * (v - estimate[None, :]) / total
    lin["stratum_id"] = f["stratum_id"]
    lin["psu_id"] = f["psu_id"]
    psu = lin.groupby(["stratum_id", "psu_id"], sort=True)[list(values.columns)].sum()
    covariance = np.zeros((v.shape[1], v.shape[1]), dtype=float)
    for _, part in psu.groupby(level=0, sort=True):
        a = part.to_numpy(dtype=float)
        m = len(a)
        if m < 2:
            if np.all(a == 0):
                continue
            raise ValueError("A contributing stratum has one public PSU; no automatic variance fix is allowed")
        centered = a - a.mean(axis=0, keepdims=True)
        covariance += (m / (m - 1)) * centered.T @ centered
    active = f.loc[values.index, ["stratum_id", "psu_id"]].drop_duplicates()
    df = int(len(active) - active["stratum_id"].nunique())
    if df <= 0:
        raise ValueError("No positive survey design degrees of freedom in this domain")
    return estimate, (covariance + covariance.T) / 2, df


def fieller_set(G, H, var_G, var_H, cov_GH, df, level=0.95, guard=1e-8):
    """Solve the unpenalized Fieller quadratic; keep unbounded/disjoint sets."""
    q = float(student_t.ppf((1 + level) / 2, df) ** 2)
    a = G * G - q * var_G
    b = -2 * (G * H - q * cov_GH)
    c = H * H - q * var_H
    scale = max(abs(G * G), abs(q * var_G), abs(G * H), abs(q * cov_GH), abs(H * H), abs(q * var_H), 1e-300)
    a, b, c = a / scale, b / scale, c / scale
    eps = 1e-12
    result = {
        "point": float(H / G) if G > guard else None,
        "point_reason": "positive_measured_gain" if G > guard else ("nonpositive_measured_gain" if G <= 0 else "small_positive_measured_gain"),
        "guard": guard,
        "confidence_level": level,
        "method": "unpenalized_Fieller_using_joint_survey_Taylor_covariance_and_t",
    }
    if abs(a) <= eps:
        if abs(b) <= eps:
            kind, intervals = ("all_real", [["-inf", "inf"]]) if c <= eps else ("empty", [])
        elif b > 0:
            kind, intervals = "lower_halfline", [["-inf", float(-c / b)]]
        else:
            kind, intervals = "upper_halfline", [[float(-c / b), "inf"]]
    else:
        disc = b * b - 4 * a * c
        disc_tol = eps * max(abs(b * b), abs(4 * a * c), 1e-300)
        if disc < -disc_tol:
            kind, intervals = ("all_real", [["-inf", "inf"]]) if a < 0 else ("empty", [])
        elif a < 0 and abs(disc) <= disc_tol:
            # A downward-opening quadratic touching zero is <=0 everywhere.
            kind, intervals = "all_real", [["-inf", "inf"]]
        else:
            root = math.sqrt(max(disc, 0.0))
            roots = sorted([(-b - root) / (2 * a), (-b + root) / (2 * a)])
            if a > 0:
                kind, intervals = "bounded", [[float(roots[0]), float(roots[1])]]
            else:
                kind, intervals = "disjoint", [["-inf", float(roots[0])], [float(roots[1]), "inf"]]
    result.update({"set_type": kind, "intervals": intervals})
    return result


def ci_summary(value, variance, df, level=0.95):
    if variance < -1e-12:
        raise ValueError("Materially negative variance")
    se = math.sqrt(max(float(variance), 0.0))
    half = float(student_t.ppf((1 + level) / 2, df)) * se
    return {"estimate": float(value), "SE": se, "lower": float(value - half), "upper": float(value + half), "confidence_level": level, "df": df}


def evaluate_domain(frame, cohort, predictions, guard=1e-8):
    c = cohort.set_index("person_id")
    p = predictions.set_index("person_id").loc[c.index]
    y = c["y_5y"].to_numpy(dtype=float)
    loss = pd.DataFrame({r: (p["p_" + r].to_numpy(dtype=float) - y) ** 2 for r in ROUTES}, index=c.index)
    means, cov, df = survey_mean_cov(frame, loss)
    route_results = {r: ci_summary(means[i], cov[i, i], df) for i, r in enumerate(ROUTES)}
    contrasts = {}
    for name, coefficients in CONTRASTS.items():
        k = np.asarray(coefficients, dtype=float)
        contrasts[name] = ci_summary(k @ means, k @ cov @ k, df)
    g = np.asarray(CONTRASTS["G_measured_gain"], dtype=float)
    h = np.asarray(CONTRASTS["H_remaining_gain"], dtype=float)
    G, H = float(g @ means), float(h @ means)
    retention = fieller_set(G, H, float(g @ cov @ g), float(h @ cov @ h), float(g @ cov @ h), df, guard=guard)
    retention["denominator_supported_positive"] = contrasts["G_measured_gain"]["lower"] > guard
    retention["headline_percentage_allowed"] = bool(retention["denominator_supported_positive"] and retention["set_type"] == "bounded")
    unweighted = loss.mean(axis=0).to_numpy(dtype=float)
    descriptive_unweighted = {
        "route_Brier": dict(zip(ROUTES, map(float, unweighted))),
        "contrasts": {name: float(np.asarray(k) @ unweighted) for name, k in CONTRASTS.items()},
        "scope": "descriptive_same_predictions_no_refit_no_inferential_interval",
    }
    weight = c["weight"].to_numpy(dtype=float)
    observed = float(np.average(y, weights=weight))
    return {
        "n": len(c), "events_5y": int(y.sum()), "weighted_event_fraction": observed,
        "df": df, "route_order": ROUTES,
        "Brier": route_results, "contrasts": contrasts, "retention": retention,
        "joint_covariance_of_route_losses": cov.tolist(),
        "unweighted_sensitivity": descriptive_unweighted,
        "scope": "pointwise_design_approximate_intervals_conditional_on_fitted_models_and_public_use_labels",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, default=ROOT / "outputs/predictions/test_predictions.csv")
    parser.add_argument("--cohort", type=Path, default=ROOT / "data/analysis_cohort.csv")
    parser.add_argument("--design-frame", type=Path, default=ROOT / "data/design_frame.csv")
    parser.add_argument("--output", type=Path, default=ROOT / "outputs/tables/fixed_substitution_summary.json")
    args = parser.parse_args()
    cohort = pd.read_csv(args.cohort, dtype={"person_id": str, "stratum_id": str, "psu_id": str})
    frame = pd.read_csv(args.design_frame, dtype={"person_id": str, "stratum_id": str, "psu_id": str})
    predictions = pd.read_csv(args.predictions, dtype={"person_id": str})
    cohort = cohort.loc[cohort["role"].eq("test")].copy()
    frame = frame.loc[frame["role"].eq("test")].copy()
    required = ["person_id"] + ["p_" + r for r in ROUTES]
    if not set(required).issubset(predictions.columns):
        raise ValueError("Required columns: " + ",".join(required))
    if predictions["person_id"].duplicated().any():
        raise ValueError("Duplicate prediction IDs")
    if set(cohort["person_id"]) != set(predictions["person_id"]):
        raise ValueError("Every route must cover exactly the predeclared final-test cohort")
    a = predictions[["p_" + r for r in ROUTES]].to_numpy(dtype=float)
    if not np.isfinite(a).all() or ((a < 0) | (a > 1)).any():
        raise ValueError("All supplied probabilities must be finite and in [0,1]")
    cfg = json.loads((ROOT / "protocol/analysis_config.json").read_text())
    guard = cfg["statistics"]["ratio_guard_brier_units"]
    result = {
        "protocol_id": cfg["protocol_id"],
        "primary_40plus": evaluate_domain(frame, cohort, predictions, guard),
        "sensitivity_40_79": evaluate_domain(frame, cohort.loc[cohort["age_cap80"] < 80], predictions, guard),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(args.output), "n": result["primary_40plus"]["n"], "status": "fixed_predictions_evaluated"}))


if __name__ == "__main__":
    main()
