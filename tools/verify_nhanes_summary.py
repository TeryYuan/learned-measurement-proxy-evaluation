#!/usr/bin/env python3
"""Check algebra in a saved NHANES summary, without participant-level data.

This verifies stored estimates against stored route losses and covariance.
It does NOT reconstruct a survey covariance, refit a model, or certify a clean
project run. No input paths, source values, or exception messages enter reports.
Scalar comparisons use atol=1e-12 and rtol=1e-10. Matrix roundoff and Fieller
degeneracies use 64 * machine epsilon relative to their own numerical scale.
"""

import argparse
import json
import math
import os
from pathlib import Path
import sys
import tempfile

try:
    import numpy as np
    from scipy.stats import t
except ImportError:
    print('{"status":"ERROR","error":"dependencies_unavailable"}')
    raise SystemExit(2)


ATOL = 1e-12
RTOL = 1e-10
ROUND = 64 * np.finfo(float).eps
ROUTES = ("B", "M", "S", "MEAN", "P")
DOMAINS = ("primary_40plus", "sensitivity_40_79")
CONTRASTS = {
    "G_measured_gain": ("B", "M"),
    "D_substitution_loss": ("S", "M"),
    "H_remaining_gain": ("B", "S"),
    "A_proxy_over_mean": ("MEAN", "S"),
    "G_P_retrained_proxy_gain": ("B", "P"),
}


def _psd(matrix):
    scale = float(np.max(np.abs(matrix)))
    return scale == 0 or bool(np.linalg.eigvalsh(matrix / scale).min() >= -ROUND)


def _zero_roundoff(value, scale):
    if not math.isfinite(value) or not math.isfinite(scale):
        return value
    return 0.0 if abs(value) <= ROUND * scale else value


def fieller_set(numerator, denominator, var_numerator, var_denominator,
                covariance, critical_value):
    """Solve the unpenalized Fieller quadratic inequality for H/G.

    Returns canonical set_type and closed intervals with +/-inf endpoints.
    bounded includes a singleton; half_line is the linear degeneracy;
    all_real and empty also cover constant and zero-discriminant cases.
    Invalid covariance or numerical overflow raises ValueError, never an
    invented finite interval. This function does not apply a point-ratio guard.
    """
    values = np.asarray([numerator, denominator, var_numerator,
                         var_denominator, covariance, critical_value], dtype=float)
    if not np.isfinite(values).all() or critical_value <= 0:
        raise ValueError("invalid_numeric_input")
    if var_numerator < 0 or var_denominator < 0:
        raise ValueError("negative_variance")
    joint = np.array([[var_numerator, covariance], [covariance, var_denominator]])
    if not _psd(joint):
        raise ValueError("invalid_covariance")
    h, g, vh, vg, chg, q = map(float, values)
    with np.errstate(over="ignore", invalid="ignore"):
        q2 = q * q
        a = _zero_roundoff(g * g - q2 * vg, abs(g * g) + abs(q2 * vg))
        b = _zero_roundoff(-2 * (h * g - q2 * chg),
                          2 * (abs(h * g) + abs(q2 * chg)))
        c = _zero_roundoff(h * h - q2 * vh, abs(h * h) + abs(q2 * vh))
    if not all(math.isfinite(x) for x in (a, b, c)):
        raise ValueError("numeric_range_exceeded")
    scale = max(abs(a), abs(b), abs(c))
    if scale:
        a, b, c = a / scale, b / scale, c / scale
    if a == 0:
        if b == 0:
            return {"set_type": "all_real" if c <= 0 else "empty",
                    "intervals": [[-math.inf, math.inf]] if c <= 0 else []}
        endpoint = -c / b
        intervals = [[-math.inf, endpoint]] if b > 0 else [[endpoint, math.inf]]
        return {"set_type": "half_line", "intervals": intervals}
    discriminant = _zero_roundoff(b * b - 4 * a * c,
                                  abs(b * b) + abs(4 * a * c))
    if discriminant < 0:
        return {"set_type": "empty" if a > 0 else "all_real",
                "intervals": [] if a > 0 else [[-math.inf, math.inf]]}
    if discriminant == 0:
        root = -b / (2 * a)
        return {"set_type": "bounded" if a > 0 else "all_real",
                "intervals": [[root, root]] if a > 0 else [[-math.inf, math.inf]]}
    # The product-of-roots form avoids cancellation in one of the roots.
    qroot = -0.5 * (b + math.copysign(math.sqrt(discriminant), b))
    roots = sorted((qroot / a, c / qroot))
    if not all(math.isfinite(x) for x in roots):
        raise ValueError("numeric_range_exceeded")
    if a > 0:
        return {"set_type": "bounded", "intervals": [roots]}
    return {"set_type": "disjoint",
            "intervals": [[-math.inf, roots[0]], [roots[1], math.inf]]}


class Check:
    def __init__(self):
        self.count = 0
        self.errors = []

    def test(self, condition, field):
        self.count += 1
        if not bool(condition) and field not in self.errors:
            self.errors.append(field)
        return bool(condition)

    def number(self, obj, key, field):
        value = obj.get(key) if isinstance(obj, dict) else None
        valid = isinstance(value, (int, float)) and not isinstance(value, bool)
        if self.test(valid and math.isfinite(value), field):
            return float(value)
        return None

    def close(self, actual, expected, field):
        return self.test(actual is not None and
                         np.isclose(actual, expected, atol=ATOL, rtol=RTOL), field)

    def interval_record(self, record, expected, variance, df, field):
        if not self.test(isinstance(record, dict), field):
            return
        point = self.number(record, "estimate", field + ".estimate")
        self.close(point, expected, field + ".estimate")
        se = self.number(record, "SE", field + ".SE")
        level = self.number(record, "confidence_level", field + ".confidence_level")
        rdf = self.number(record, "df", field + ".df")
        self.close(level, 0.95, field + ".confidence_level")
        self.close(rdf, df, field + ".df")
        lower = self.number(record, "lower", field + ".lower")
        upper = self.number(record, "upper", field + ".upper")
        if variance is not None:
            expected_se = math.sqrt(variance)
            self.close(se, expected_se, field + ".SE")
            delta = float(t.ppf(0.975, df)) * expected_se
            self.close(lower, expected - delta, field + ".lower")
            self.close(upper, expected + delta, field + ".upper")

    def report(self):
        return {
            "schema_version": 1,
            "status": "FAIL" if self.errors else "PASS",
            "scope": "saved_summary_algebra_only",
            "participant_level_reconstruction": False,
            "survey_covariance_reconstruction": False,
            "clean_environment_project_acceptance": False,
            "checks_performed": self.count,
            "mismatch_count": len(self.errors),
            "mismatches": sorted(self.errors),
            "tolerance": {"absolute": ATOL, "relative": RTOL,
                          "roundoff_scale_factor": 64},
        }


def _covariance(check, domain, field):
    try:
        matrix = np.asarray(domain.get("joint_covariance_of_route_losses"), dtype=float)
    except (TypeError, ValueError, OverflowError):
        check.test(False, field + ".numeric")
        return None
    if not check.test(matrix.shape == (5, 5), field + ".shape"):
        return None
    if not check.test(np.isfinite(matrix).all(), field + ".finite"):
        return None
    scale = float(np.max(np.abs(matrix)))
    if not check.test(np.max(np.abs(matrix - matrix.T)) <= ROUND * scale,
                      field + ".symmetric"):
        return None
    # Symmetrization is limited to matrix roundoff, after the symmetry check.
    matrix = (matrix + matrix.T) / 2
    if not check.test(_psd(matrix), field + ".positive_semidefinite"):
        return None
    return matrix


def _variance(matrix, vector):
    # A PSD matrix can produce a tiny negative result through cancellation.
    return max(0.0, float(vector @ matrix @ vector))


def _check_set(check, stored, expected, field):
    aliases = {"unbounded": "all_real", "linear": "half_line",
               "half-line": "half_line", "singleton": "bounded"}
    actual_type = stored.get("set_type")
    if not isinstance(actual_type, str):
        actual_type = None
    check.test(aliases.get(actual_type, actual_type) == expected["set_type"],
               field + ".set_type")
    intervals = stored.get("intervals")
    if not check.test(isinstance(intervals, list) and
                      len(intervals) == len(expected["intervals"]), field + ".intervals"):
        return
    for index, (actual, wanted) in enumerate(zip(intervals, expected["intervals"])):
        label = field + ".intervals." + str(index)
        if not check.test(isinstance(actual, list) and len(actual) == 2, label):
            continue
        for side, target in enumerate(wanted):
            value = actual[side]
            # JSON null may represent an infinite endpoint, determined by side.
            if value is None:
                value = -math.inf if side == 0 else math.inf
            numeric = isinstance(value, (int, float)) and not isinstance(value, bool)
            same = numeric and ((value == target) if math.isinf(target) else
                                 (math.isfinite(value) and
                                  np.isclose(value, target, atol=ATOL, rtol=RTOL)))
            check.test(same, label + (".lower" if side == 0 else ".upper"))


def _check_domain(check, domain, label):
    if not check.test(isinstance(domain, dict), label):
        return
    order = domain.get("route_order")
    valid_order = (isinstance(order, list) and len(order) == 5 and
                   all(isinstance(x, str) for x in order) and set(order) == set(ROUTES))
    if not check.test(valid_order, label + ".route_order"):
        return
    df = check.number(domain, "df", label + ".df")
    if not check.test(df is not None and df > 0, label + ".df"):
        return
    n = check.number(domain, "n", label + ".n")
    events = check.number(domain, "events_5y", label + ".events_5y")
    check.test(n is not None and n > 0 and n.is_integer(), label + ".n")
    check.test(events is not None and events.is_integer() and n is not None and
               0 <= events <= n, label + ".events_5y")
    fraction = check.number(domain, "weighted_event_fraction", label + ".weighted_event_fraction")
    check.test(fraction is not None and 0 <= fraction <= 1,
               label + ".weighted_event_fraction")
    matrix = _covariance(check, domain, label + ".joint_covariance_of_route_losses")
    brier = domain.get("Brier")
    if not check.test(isinstance(brier, dict), label + ".Brier"):
        return
    means = [check.number(brier.get(route), "estimate", label + ".Brier." + route + ".estimate")
             for route in order]
    if any(value is None for value in means):
        return
    means = np.array(means)
    for i, route in enumerate(order):
        check.test(0 <= means[i] <= 1, label + ".Brier." + route + ".estimate_range")
        variance = _variance(matrix, np.eye(5)[i]) if matrix is not None else None
        check.interval_record(brier[route], means[i], variance, df, label + ".Brier." + route)
    contrasts = domain.get("contrasts")
    if not check.test(isinstance(contrasts, dict), label + ".contrasts"):
        return
    vectors, expected = {}, {}
    for name, (positive, negative) in CONTRASTS.items():
        vector = np.zeros(5)
        vector[order.index(positive)] = 1
        vector[order.index(negative)] = -1
        vectors[name] = vector
        expected[name] = float(vector @ means)
        variance = _variance(matrix, vector) if matrix is not None else None
        check.interval_record(contrasts.get(name), expected[name], variance, df,
                              label + ".contrasts." + name)
    check.close(expected["H_remaining_gain"],
                expected["G_measured_gain"] - expected["D_substitution_loss"],
                label + ".identity.H_equals_G_minus_D")
    retention = domain.get("retention")
    if check.test(isinstance(retention, dict), label + ".retention"):
        field = label + ".retention"
        guard = check.number(retention, "guard", field + ".guard")
        check.close(guard, 1e-8, field + ".guard")
        level = check.number(retention, "confidence_level", field + ".confidence_level")
        check.close(level, 0.95, field + ".confidence_level")
        g, h = expected["G_measured_gain"], expected["H_remaining_gain"]
        if guard is not None and guard > 0:
            if g > guard:
                point = check.number(retention, "point", field + ".point")
                check.close(point, h / g, field + ".point")
            else:
                check.test(retention.get("point") is None, field + ".point")
        if matrix is not None:
            cg, ch = vectors["G_measured_gain"], vectors["H_remaining_gain"]
            vg, vh = _variance(matrix, cg), _variance(matrix, ch)
            q = float(t.ppf(0.975, df))
            supported = bool(g - q * math.sqrt(vg) > 0)
            check.test(retention.get("denominator_supported_positive") is supported,
                       field + ".denominator_supported_positive")
            allowed = supported and guard is not None and g > guard
            check.test(retention.get("headline_percentage_allowed") is bool(allowed),
                       field + ".headline_percentage_allowed")
            try:
                result = fieller_set(h, g, vh, vg, float(ch @ matrix @ cg), q)
                _check_set(check, retention, result, field)
            except (ValueError, OverflowError, FloatingPointError):
                check.test(False, field + ".fieller_numeric_domain")
    unweighted = domain.get("unweighted_sensitivity")
    if not check.test(isinstance(unweighted, dict), label + ".unweighted_sensitivity"):
        return
    unweighted_routes = unweighted.get("route_Brier")
    uvalues = {route: check.number(unweighted_routes, route,
               label + ".unweighted_sensitivity.route_Brier." + route) for route in ROUTES}
    if any(value is None for value in uvalues.values()):
        return
    for route, value in uvalues.items():
        check.test(0 <= value <= 1, label + ".unweighted_sensitivity.route_Brier." + route)
    for name, (positive, negative) in CONTRASTS.items():
        field = label + ".unweighted_sensitivity.contrasts." + name
        observed = check.number(unweighted.get("contrasts"), name, field)
        check.close(observed, uvalues[positive] - uvalues[negative], field)


def verify_summary(payload):
    """Return a sanitized report; all field labels come from this code."""
    check = Check()
    if check.test(isinstance(payload, dict), "summary.object"):
        for label in DOMAINS:
            _check_domain(check, payload.get(label), label)
    return check.report()


class SafeParser(argparse.ArgumentParser):
    def error(self, message):
        print('{"status":"ERROR","error":"arguments_invalid"}')
        raise SystemExit(2)


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate_json_key")
        result[key] = value
    return result


def _reject_constant(value):
    raise ValueError("nonfinite_json_constant")


def main(argv=None):
    parser = SafeParser(prog="verify_nhanes_summary.py", description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    report = None
    try:
        if args.input.resolve() == args.output.resolve():
            print('{"status":"ERROR","error":"input_output_collision"}')
            return 2
        with args.input.open(encoding="utf-8-sig") as handle:
            payload = json.load(handle, object_pairs_hook=_unique_object,
                                parse_constant=_reject_constant)
    except (OSError, ValueError, TypeError, RecursionError):
        report = Check().report()
        report.update(status="ERROR", error="input_unreadable_or_invalid")
    if report is None:
        try:
            report = verify_summary(payload)
        except Exception:
            # Fail closed without copying an input-derived exception into a log.
            report = Check().report()
            report.update(status="ERROR", error="verification_internal_error")
    temporary = None
    try:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8",
                                         dir=args.output.parent, prefix="summary-check-",
                                         suffix=".tmp", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(report, handle, ensure_ascii=True, indent=2, allow_nan=False)
            handle.write("\n")
        os.replace(temporary, args.output)
    except (OSError, ValueError, TypeError):
        print('{"status":"ERROR","error":"report_write_failed"}')
        return 2
    finally:
        if temporary is not None:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass
    fields = ("status", "checks_performed", "mismatch_count", "mismatches", "error")
    print(json.dumps({key: report[key] for key in fields if key in report}))
    return {"PASS": 0, "FAIL": 1, "ERROR": 2}[report["status"]]


if __name__ == "__main__":
    raise SystemExit(main())
