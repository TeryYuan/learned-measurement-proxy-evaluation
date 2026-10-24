"""A finite check of scientific boundary cases, using artificial values only."""
import numpy as np
import pandas as pd
from evaluate_fixed_predictions import survey_mean_cov, fieller_set


def main():
    frame = pd.DataFrame({"person_id": ["a", "b", "c", "d"], "weight": [1., 1., 1., 1.], "stratum_id": ["h1", "h1", "h2", "h2"], "psu_id": ["p1", "p2", "p3", "p4"]})
    values = pd.DataFrame({"loss": [0., .2, .4, .6]}, index=["a", "b", "c", "d"])
    # Hand-computed: each stratum contributes .0025, total variance .005.
    m, cov, df = survey_mean_cov(frame, values)
    assert np.allclose(m, [.3]) and np.allclose(cov, [[.005]]) and df == 2
    # A domain-excluded PSU still supplies its zero linearized contribution.
    m, cov, df = survey_mean_cov(frame, values.loc[["b", "c", "d"]])
    assert np.allclose(m, [.4]) and np.allclose(cov, [[2 / 225]]) and df == 1
    # For identical measured and substituted losses, robust G implies R=1.
    f = fieller_set(.02, .02, 1e-6, 1e-6, 1e-6, 30)
    assert f["set_type"] == "bounded" and np.allclose(f["intervals"], [[1., 1.]])
    # A noisy denominator must not be converted into a finite percentage band.
    f = fieller_set(.001, .2, 1e-4, 1e-4, 0., 30)
    assert f["set_type"] == "disjoint"
    f = fieller_set(0., .01, 1e-4, 1e-4, 0., 30)
    assert f["point"] is None and f["set_type"] == "all_real"
    # If H=G but G is not separated from zero, touching half-lines form all R.
    f = fieller_set(.001, .001, 1e-4, 1e-4, 1e-4, 30)
    assert f["set_type"] == "all_real"
    print("PASS: hand-calculated design variance, domain handling, exact substitution identity, and unstable Fieller sets; no study outcomes used")


if __name__ == "__main__":
    main()
