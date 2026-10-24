#!/usr/bin/env python3
"""Prepare the prespecified NHANES cohort and report counts only.

Dependencies: pandas and numpy.  Run from any working directory:
    python /path/to/kit/scripts/prepare_cohort.py

This script does not fit a model, impute predictors, examine feature-outcome
associations, or select a task from its observed performance.  Public-use
mortality timing may be perturbed for selected records by NCHS.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
CYCLES = (
    (2003, "C", "L25_C", "train", "development"),
    (2005, "D", "CBC_D", "validation", "development"),
    (2007, "E", "CBC_E", "test", "test"),
    (2009, "F", "CBC_F", "test", "test"),
)
RAW_COLUMNS = {
    "DEMO": [
        "SEQN", "RIDSTATR", "RIDAGEYR", "RIAGENDR", "RIDRETH1",
        "DMDEDUC2", "INDFMPIR", "WTMEC2YR", "SDMVSTRA", "SDMVPSU",
    ],
    "BMX": ["SEQN", "BMXBMI", "BMXWAIST"],
    "BPX": [
        "SEQN", "BPXSY1", "BPXSY2", "BPXSY3", "BPXSY4",
        "BPXDI1", "BPXDI2", "BPXDI3", "BPXDI4", "BPXPLS",
    ],
    "SMQ": ["SEQN", "SMQ020", "SMQ040"],
    "DIQ": ["SEQN", "DIQ010"],
    "MCQ": [
        "SEQN", "MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E",
        "MCQ160F", "MCQ220",
    ],
    "BPQ": ["SEQN", "BPQ020"],
    "HUQ": ["SEQN", "HUQ010"],
    "CBC": ["SEQN", "LBXRDW", "LBXHGB"],
}
X_NUMERIC = [
    "age_cap80", "income_poverty_ratio", "bmi", "waist_cm",
    "sbp_mmhg", "dbp_mmhg", "pulse_bpm",
]
X_CATEGORICAL = [
    "sex", "race", "education", "smoking", "general_health",
    "diabetes", "cvd", "hypertension", "cancer",
]
Z_COLUMNS = ["rdw_pct", "hemoglobin_g_dl"]
YES_NO = {1: "yes", 2: "no"}


def finite_positive(series: pd.Series) -> pd.Series:
    return pd.Series(np.isfinite(series), index=series.index) & series.gt(0)


def valid_design_number(series: pd.Series) -> pd.Series:
    return finite_positive(series) & series.eq(np.floor(series))


def endpoint_5y(mortstat: pd.Series, months: pd.Series) -> pd.Series:
    """Death at month 0 or 60 is an event; alive before month 60 is unknown."""
    valid = mortstat.isin([0, 1]) & np.isfinite(months) & months.ge(0)
    case = valid & mortstat.eq(1) & months.le(60)
    control = valid & (
        months.gt(60) | (mortstat.eq(0) & months.ge(60))
    )
    y = pd.Series(np.nan, index=mortstat.index, dtype="float64")
    y.loc[case] = 1.0
    y.loc[control] = 0.0
    return y


def smoking_status(ever: pd.Series, current: pd.Series) -> pd.Series:
    """Honor the structural skip of SMQ040 in respondents with SMQ020 == 2."""
    result = pd.Series(pd.NA, index=ever.index, dtype="string")
    result.loc[ever.eq(2)] = "never"
    result.loc[ever.eq(1) & current.isin([1, 2])] = "current"
    result.loc[ever.eq(1) & current.eq(3)] = "former"
    return result


def categorical(series: pd.Series, mapping: dict[int, str]) -> pd.Series:
    # Field-specific maps leave refusals, unknowns and missing values missing.
    # Never recode 7/9 globally: they can be valid continuous measurements.
    return series.map(mapping).astype("string")


def positive_mean(frame: pd.DataFrame) -> pd.Series:
    valid = np.isfinite(frame) & frame.gt(0)
    return frame.where(valid).mean(axis=1)


def boundary_checks() -> None:
    mortality = pd.Series([1, 1, 0, 0, 1, 1, 9, 1], dtype=float)
    months = pd.Series([0, 60, 59, 60, 61, -1, 100, np.nan], dtype=float)
    expected = np.array([1, 1, np.nan, 0, 0, np.nan, np.nan, np.nan])
    assert np.allclose(endpoint_5y(mortality, months), expected, equal_nan=True)
    actual = smoking_status(
        pd.Series([2, 1, 1, 1, 1, 9], dtype=float),
        pd.Series([np.nan, 1, 2, 3, np.nan, np.nan], dtype=float),
    ).fillna("missing").tolist()
    assert actual == ["never", "current", "current", "former", "missing", "missing"]


def validate_seqn(frame: pd.DataFrame, filename: str) -> pd.DataFrame:
    if not valid_design_number(frame["SEQN"]).all():
        raise ValueError(f"Invalid SEQN in {filename}")
    frame["SEQN"] = frame["SEQN"].astype("int64")
    if frame["SEQN"].duplicated().any():
        raise ValueError(f"Duplicate SEQN in {filename}; one-to-one merge required")
    return frame


def read_xpt(path: Path, component: str) -> tuple[pd.DataFrame, dict]:
    frame = pd.read_sas(path, format="xport")
    required = RAW_COLUMNS[component]
    absent = [column for column in required if column not in frame]
    if absent:
        raise ValueError(f"{path.name} is missing required columns: {absent}")
    frame = frame[required].copy().apply(pd.to_numeric, errors="raise")
    frame = validate_seqn(frame, path.name)
    metadata = {
        "file": path.name, "format": "SAS_XPORT", "rows": int(len(frame)),
        "required_columns": required, "missing_required_columns": [],
    }
    return frame, metadata


def read_mortality(path: Path) -> tuple[pd.DataFrame, dict]:
    # NCHS R_ReadInProgramAllSurveys.R uses 1-based inclusive positions:
    # SEQN 1:6; ELIGSTAT 15; MORTSTAT 16; PERMTH_EXM 46:48.
    # Do not load death-certificate DIABETES or HYPERTEN as baseline predictors.
    columns = ["SEQN", "ELIGSTAT", "MORTSTAT", "PERMTH_EXM"]
    frame = pd.read_fwf(
        path, colspecs=[(0, 6), (14, 15), (15, 16), (45, 48)],
        names=columns, na_values=["", "."], keep_default_na=True,
    ).apply(pd.to_numeric, errors="raise")
    frame = validate_seqn(frame, path.name)
    return frame, {
        "file": path.name, "format": "NCHS_fixed_width", "rows": int(len(frame)),
        "required_columns": columns, "missing_required_columns": [],
    }


def masks_for(raw: pd.DataFrame) -> list[tuple[str, pd.Series]]:
    return [
        ("mec_examined", raw["RIDSTATR"].eq(2)),
        ("positive_mec_weight", finite_positive(raw["WTMEC2YR"])),
        ("valid_survey_design", valid_design_number(raw["SDMVSTRA"])
         & valid_design_number(raw["SDMVPSU"])),
        ("age_40_or_older", np.isfinite(raw["RIDAGEYR"]) & raw["RIDAGEYR"].ge(40)),
        ("mortality_linkage_eligible", raw["ELIGSTAT"].eq(1)),
        ("both_cbc_measurements_observed_positive",
         finite_positive(raw["LBXRDW"]) & finite_positive(raw["LBXHGB"])),
        ("valid_mortality_status", raw["MORTSTAT"].isin([0, 1])),
        ("valid_nonnegative_mec_followup",
         np.isfinite(raw["PERMTH_EXM"]) & raw["PERMTH_EXM"].ge(0)),
        ("five_year_outcome_determinable",
         endpoint_5y(raw["MORTSTAT"], raw["PERMTH_EXM"]).notna()),
    ]


def build_frame(raw: pd.DataFrame, cycle: str, role: str, period: str) -> tuple[pd.DataFrame, dict]:
    conditions = masks_for(raw)
    keep = pd.Series(True, index=raw.index)
    steps = [{"step": "demographic_records", "n_remaining": int(len(raw)), "n_excluded": 0}]
    base_mask = None
    for name, condition in conditions:
        previous_n = int(keep.sum())
        keep &= condition.fillna(False)
        steps.append({
            "step": name, "n_remaining": int(keep.sum()),
            "n_excluded": previous_n - int(keep.sum()),
        })
        if name == "valid_survey_design":
            base_mask = keep.copy()
    assert base_mask is not None

    # Keep every age in the design frame.  Age >=40 and subsequent criteria
    # define the analysis domain; dropping them from the frame can lose PSUs.
    source = raw.loc[base_mask].copy()
    frame = pd.DataFrame(index=source.index)
    frame["person_id"] = cycle + "_" + source["SEQN"].astype(str)
    frame["SEQN"] = source["SEQN"]
    frame["cycle"] = cycle
    frame["period"] = period
    frame["role"] = role
    frame["SDMVSTRA"] = source["SDMVSTRA"].astype("int64")
    frame["SDMVPSU"] = source["SDMVPSU"].astype("int64")
    frame["stratum_id"] = cycle + "_" + frame["SDMVSTRA"].astype(str)
    frame["psu_id"] = frame["stratum_id"] + "_" + frame["SDMVPSU"].astype(str)
    frame["WTMEC2YR"] = source["WTMEC2YR"]
    # Development and test each pool two 2-year cycles.  Training algorithms
    # can subsequently normalize these weights inside their own fit splits.
    frame["weight"] = source["WTMEC2YR"] / 2.0
    frame["is_eligible"] = keep.loc[source.index].astype(bool)
    frame["exclusion_reason"] = "included"
    for name, condition in conditions[3:]:
        first_failure = frame["exclusion_reason"].eq("included") & ~condition.loc[source.index].fillna(False)
        frame.loc[first_failure, "exclusion_reason"] = name

    frame["RIDAGEYR"] = source["RIDAGEYR"]
    frame["ELIGSTAT"] = source["ELIGSTAT"]
    frame["MORTSTAT"] = source["MORTSTAT"]
    frame["PERMTH_EXM"] = source["PERMTH_EXM"]
    frame["y_5y"] = endpoint_5y(source["MORTSTAT"], source["PERMTH_EXM"]).astype("Int64")

    frame["age_cap80"] = source["RIDAGEYR"].clip(upper=80)
    income = source["INDFMPIR"]
    frame["income_poverty_ratio"] = income.where(np.isfinite(income) & income.between(0, 5))
    for output, original in [("bmi", "BMXBMI"), ("waist_cm", "BMXWAIST"), ("pulse_bpm", "BPXPLS")]:
        frame[output] = source[original].where(finite_positive(source[original]))
    frame["sbp_mmhg"] = positive_mean(source[[f"BPXSY{i}" for i in range(1, 5)]])
    frame["dbp_mmhg"] = positive_mean(source[[f"BPXDI{i}" for i in range(1, 5)]])
    frame["sex"] = categorical(source["RIAGENDR"], {1: "male", 2: "female"})
    frame["race"] = categorical(source["RIDRETH1"], {
        1: "mexican_american", 2: "other_hispanic", 3: "non_hispanic_white",
        4: "non_hispanic_black", 5: "other_including_multiracial",
    })
    frame["education"] = categorical(source["DMDEDUC2"], {
        1: "less_than_9th_grade", 2: "9_to_11th_grade_no_diploma",
        3: "high_school_or_GED", 4: "some_college_or_AA", 5: "college_graduate_or_above",
    })
    frame["smoking"] = smoking_status(source["SMQ020"], source["SMQ040"])
    frame["general_health"] = categorical(source["HUQ010"], {
        1: "excellent", 2: "very_good", 3: "good", 4: "fair", 5: "poor",
    })
    frame["diabetes"] = categorical(source["DIQ010"], {1: "yes", 2: "no", 3: "borderline"})
    cvd = source[["MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E", "MCQ160F"]]
    frame["cvd"] = pd.Series(pd.NA, index=source.index, dtype="string")
    frame.loc[cvd.eq(1).any(axis=1), "cvd"] = "yes"
    frame.loc[cvd.eq(2).all(axis=1), "cvd"] = "no"
    frame["hypertension"] = categorical(source["BPQ020"], YES_NO)
    frame["cancer"] = categorical(source["MCQ220"], YES_NO)
    frame["rdw_pct"] = source["LBXRDW"].where(finite_positive(source["LBXRDW"]))
    frame["hemoglobin_g_dl"] = source["LBXHGB"].where(finite_positive(source["LBXHGB"]))

    metadata_columns = [
        "person_id", "SEQN", "cycle", "period", "role", "SDMVSTRA", "SDMVPSU",
        "stratum_id", "psu_id", "WTMEC2YR", "weight", "is_eligible", "exclusion_reason",
        "RIDAGEYR", "ELIGSTAT", "MORTSTAT", "PERMTH_EXM", "y_5y",
    ]
    frame = frame[metadata_columns + X_NUMERIC + X_CATEGORICAL + Z_COLUMNS]
    return frame.reset_index(drop=True), {
        "cycle": cycle, "role": role, "period": period,
        "steps": steps, "design_frame_n": int(base_mask.sum()),
        "analysis_domain_n": int(keep.sum()),
    }


def group_counts(frame: pd.DataFrame, column: str) -> dict:
    result = {}
    for value, group in frame.groupby(column, sort=False):
        eligible = group.loc[group["is_eligible"]]
        eligible_psus = int(eligible["psu_id"].nunique())
        eligible_strata = int(eligible["stratum_id"].nunique())
        result[str(value)] = {
            "design_frame_n": int(len(group)),
            "age_40_or_older_in_frame_n": int(group["RIDAGEYR"].ge(40).sum()),
            "eligible_n": int(len(eligible)),
            "events_within_5_years_n": int(eligible["y_5y"].eq(1).sum()),
            "controls_at_5_years_n": int(eligible["y_5y"].eq(0).sum()),
            "design_frame_strata_n": int(group["stratum_id"].nunique()),
            "design_frame_psus_n": int(group["psu_id"].nunique()),
            "eligible_strata_n": eligible_strata,
            "eligible_psus_n": eligible_psus,
            "eligible_psus_minus_strata_df": eligible_psus - eligible_strata,
        }
    return result


def feature_missingness(frame: pd.DataFrame) -> dict:
    return {
        column: {"missing_n": int(frame[column].isna().sum()), "denominator_n": int(len(frame))}
        for column in X_NUMERIC + X_CATEGORICAL + Z_COLUMNS
    }


def write_json(path: Path, content: dict) -> None:
    path.write_text(json.dumps(content, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")


def main() -> None:
    boundary_checks()
    manifest = json.loads((ROOT / "sources/source_manifest.json").read_text(encoding="utf-8"))
    paths = {item["file"]: ROOT / item["local_path"] for item in manifest["data_files"]}
    dictionary = json.loads((ROOT / "protocol/feature_dictionary.json").read_text(encoding="utf-8"))
    assert dictionary["predictor_columns"]["numeric"] == X_NUMERIC
    assert dictionary["predictor_columns"]["categorical"] == X_CATEGORICAL
    assert dictionary["measurement_columns"] == Z_COLUMNS

    frames, flows, file_checks = [], [], []
    for year, suffix, cbc_stem, role, period in CYCLES:
        cycle = f"{year}-{year + 1}"
        raw, metadata = read_xpt(paths[f"DEMO_{suffix}.xpt"], "DEMO")
        file_checks.append(metadata)
        identifiers = set(raw["SEQN"])
        for component in ["BMX", "BPX", "SMQ", "DIQ", "MCQ", "BPQ", "HUQ", "CBC"]:
            filename = f"{cbc_stem}.xpt" if component == "CBC" else f"{component}_{suffix}.xpt"
            data, metadata = read_xpt(paths[filename], component)
            extra = set(data["SEQN"]) - identifiers
            if extra:
                raise ValueError(f"{filename} has {len(extra)} identifiers absent from its DEMO file")
            metadata["unmatched_to_demographics_n"] = 0
            file_checks.append(metadata)
            raw = raw.merge(data, on="SEQN", how="left", validate="one_to_one", sort=False)
        filename = f"NHANES_{year}_{year + 1}_MORT_2019_PUBLIC.dat"
        mortality, metadata = read_mortality(paths[filename])
        if set(mortality["SEQN"]) != identifiers:
            raise ValueError(f"{filename} does not contain exactly the DEMO identifiers")
        metadata["unmatched_to_demographics_n"] = 0
        file_checks.append(metadata)
        raw = raw.merge(mortality, on="SEQN", how="left", validate="one_to_one", sort=False)
        frame, flow = build_frame(raw, cycle, role, period)
        frames.append(frame)
        flows.append(flow)

    design = pd.concat(frames, ignore_index=True)
    analysis = design.loc[design["is_eligible"]].copy().reset_index(drop=True)
    if design["person_id"].duplicated().any():
        raise ValueError("Duplicate person_id after combining cycles")
    assert analysis["RIDAGEYR"].ge(40).all()
    assert analysis["y_5y"].isin([0, 1]).all()
    assert analysis[Z_COLUMNS].notna().all().all()

    (ROOT / "data").mkdir(exist_ok=True)
    (ROOT / "data_checks").mkdir(exist_ok=True)
    design.to_csv(ROOT / "data/design_frame.csv", index=False, na_rep="")
    analysis.to_csv(ROOT / "data/analysis_cohort.csv", index=False, na_rep="")

    generated = datetime.now(timezone.utc).isoformat(timespec="seconds")
    cohort_flow = {
        "generated_utc": generated,
        "stage": "prespecified_cohort_preparation_counts_only",
        "model_fitting_performed": False,
        "feature_outcome_associations_examined": False,
        "flow_order": "MEC -> positive weight -> valid design (full frame) -> age >=40 -> linkage -> measured Z -> mortality status -> follow-up -> 5-year outcome",
        "cycles": flows,
        "total_design_frame_n": int(len(design)),
        "total_analysis_domain_n": int(len(analysis)),
        "excluded_domain_n": int((~design["is_eligible"]).sum()),
        "exclusion_reasons_in_design_frame": {
            str(key): int(value) for key, value in design["exclusion_reason"].value_counts().items()
        },
    }
    profile = {
        "generated_utc": generated,
        "stage": "feasibility_counts_and_missingness_only",
        "model_fitting_performed": False,
        "feature_outcome_associations_examined": False,
        "effect_estimates_computed": False,
        "input_file_count": int(len(file_checks)),
        "input_file_checks": file_checks,
        "predictor_count": 16,
        "predictor_columns": {"numeric": X_NUMERIC, "categorical": X_CATEGORICAL},
        "measurement_columns": Z_COLUMNS,
        "outcome_column": "y_5y",
        "prediction_columns_present": False,
        "counts_by_cycle": group_counts(design, "cycle"),
        "counts_by_role": group_counts(design, "role"),
        "counts_by_period": group_counts(design, "period"),
        "feature_missingness_before_imputation_in_analysis_cohort": {
            "all": feature_missingness(analysis),
            "by_role": {
                str(role): feature_missingness(group)
                for role, group in analysis.groupby("role", sort=False)
            },
        },
        "measurements_missing_in_full_design_frame_n": {
            column: int(design[column].isna().sum()) for column in Z_COLUMNS
        },
        "survey_domain_note": "The frame contains all ages with MEC examination, positive weight and valid design IDs. Rows outside is_eligible receive zero linearized contribution in later domain variance estimation; do not fit or invent predictions for them.",
        "weight_note": "weight = WTMEC2YR / 2 for each four-year development/test pool. Derive single-cycle (/1) or all-four-cycle (/4) descriptive weights from retained WTMEC2YR; all fits normalize current-pool weights to current-fit mean 1.",
        "temporal_note": "Roles separate recruitment cycles in a retrospective study; five-year training outcomes were not all available before the first test recruitment date.",
        "public_mortality_note": "NCHS perturbs follow-up time for selected public-use records; vital status is not perturbed. No individual perturbation corrections were attempted.",
    }
    write_json(ROOT / "data_checks/cohort_flow.json", cohort_flow)
    write_json(ROOT / "data_checks/feasibility_profile.json", profile)
    print(json.dumps({
        "stage": profile["stage"],
        "counts_by_role": profile["counts_by_role"],
        "counts_by_period": profile["counts_by_period"],
        "model_fitting_performed": False,
        "effect_estimates_computed": False,
    }, indent=2))


if __name__ == "__main__":
    main()
