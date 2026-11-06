"""Phase 5 (Agent 1): descriptive METABRIC survival tools. Hand-computed worked examples and data-coding checks."""
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy import stats

from dormancy.km_descriptive import (administrative_censor, holm, km_curve, km_point, landmark_incidence, logrank_test,
                                     pairwise_logrank, reliable_horizon_mask)

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

# worked example: times 1,2,3,4,5,6 ; events 1,0,1,1,0,1  ->  S(1)=5/6, S(3)=5/6*3/4... see below
T = np.array([1., 2., 3., 4., 5., 6.])
E = np.array([1, 0, 1, 1, 0, 1])


def test_km_hand_computed():
    # at t=1: n=6,d=1 -> 5/6 ; t=2 censored ; t=3: n=4,d=1 -> *3/4 ; t=4: n=3,d=1 -> *2/3 ; t=6: n=1,d=1 -> *0
    assert km_point(T, E, 1)["S"] == pytest.approx(5 / 6)
    assert km_point(T, E, 2.5)["S"] == pytest.approx(5 / 6)
    assert km_point(T, E, 3)["S"] == pytest.approx(5 / 6 * 3 / 4)
    assert km_point(T, E, 4)["S"] == pytest.approx(5 / 6 * 3 / 4 * 2 / 3)
    assert km_point(T, E, 6)["S"] == pytest.approx(0.0)
    r = km_point(T, E, 4)
    assert r["cum_inc"] == pytest.approx(1 - r["S"])
    assert r["ci_lo"] <= r["cum_inc"] <= r["ci_hi"]


def test_km_greenwood_se():
    # Greenwood: var = S^2 * sum d/(n(n-d))
    r = km_point(T, E, 4)
    S = 5 / 6 * 3 / 4 * 2 / 3
    var = S ** 2 * (1 / (6 * 5) + 1 / (4 * 3) + 1 / (3 * 2))
    assert r["se"] == pytest.approx(np.sqrt(var))


def test_censoring_not_counted_as_event():
    r = km_point(T, np.zeros(6, int), 6)
    assert r["S"] == 1.0 and r["n_events_by_horizon"] == 0


def test_at_risk_and_event_counts_at_horizon():
    r = km_point(T, E, 3)
    assert r["n_at_risk"] == 4          # times >= 3
    assert r["n_events_by_horizon"] == 2  # events at 1 and 3


def test_administrative_censoring():
    t2, e2 = administrative_censor([1, 5, 10], [1, 1, 1], 6)
    assert list(t2) == [1, 5, 6] and list(e2) == [1, 1, 0]
    # an event exactly at the horizon is kept
    t3, e3 = administrative_censor([6.0], [1], 6.0)
    assert e3[0] == 1


def test_horizon_months_equal_12_times_years():
    rng = np.random.default_rng(1)
    t = rng.exponential(100, 500); e = rng.integers(0, 2, 500)
    for y in (5, 10, 15, 20):
        a = km_point(t, e, 12 * y)
        assert a["n_at_risk"] == int(np.sum(t >= 12 * y))


def test_km_curve_monotone_and_ends_at_tmax():
    rng = np.random.default_rng(2)
    t = rng.exponential(50, 300); e = (rng.random(300) < 0.7).astype(int)
    tc, S, lo, hi = km_curve(t, e, t_max=100)
    assert tc[0] == 0 and S[0] == 1 and tc.max() <= 100
    assert np.all(np.diff(S) <= 1e-12)
    assert np.all(lo <= S + 1e-12) and np.all(S <= hi + 1e-12)


def test_reliability_rule():
    assert reliable_horizon_mask(50, 100, 5) and reliable_horizon_mask(50, 100, 10)
    assert not reliable_horizon_mask(49, 500, 5)
    assert not reliable_horizon_mask(100, 99, 10)
    assert reliable_horizon_mask(50, 30, 20) and not reliable_horizon_mask(50, 29, 20)
    with pytest.raises(ValueError):
        reliable_horizon_mask(100, 100, 15)


def test_landmark_conditional():
    # patients with time > 3: times 4,5,6 -> shifted 1,2,3 ; events 1,0,1
    r = landmark_incidence(T, E, 3.0, 3.0)
    assert r["n_landmark"] == 3
    # S(3) shifted: t=1 n=3 d=1 ->2/3 ; t=2 censored ; t=3 n=1 d=1 -> 0
    assert r["S"] == pytest.approx(0.0)
    r2 = landmark_incidence(T, E, 3.0, 1.0)
    assert r2["cum_inc"] == pytest.approx(1 / 3)


def test_logrank_identical_groups_p_near_one():
    rng = np.random.default_rng(3)
    t = rng.exponential(10, 400); e = (rng.random(400) < 0.8).astype(int)
    g = np.array(["a"] * 200 + ["b"] * 200)
    # make the two groups exact copies of one sample
    t[200:] = t[:200]; e[200:] = e[:200]
    chi2, df, p, O, Ex = logrank_test(t, e, g)
    assert chi2 == pytest.approx(0.0, abs=1e-9) and p == pytest.approx(1.0, abs=1e-6)


def test_logrank_detects_difference_and_matches_scipy_when_two_groups():
    rng = np.random.default_rng(4)
    t = np.concatenate([rng.exponential(10, 300), rng.exponential(25, 300)]); e = np.ones(600, int)
    g = np.array(["a"] * 300 + ["b"] * 300)
    chi2, df, p, O, Ex = logrank_test(t, e, g)
    assert p < 1e-6 and df == 1
    assert O.sum() == pytest.approx(Ex.sum())        # observed total = expected total
    # permutation sanity: shuffled labels give a non-significant result
    chi2s, _, ps, _, _ = logrank_test(t, e, rng.permutation(g))
    assert ps > 0.001


def test_holm_adjustment():
    adj = holm([0.01, 0.04, 0.03])
    assert adj == pytest.approx([0.03, 0.06, 0.06])
    assert np.all(adj >= np.array([0.01, 0.04, 0.03]))
    assert holm([0.5])[0] == 0.5 and holm([0.9, 0.9, 0.9]).max() <= 1.0


def test_pairwise_holm_in_range():
    rng = np.random.default_rng(5)
    t = rng.exponential(10, 300); e = np.ones(300, int); g = np.repeat(["a", "b", "c"], 100)
    rows = pairwise_logrank(t, e, g)
    assert len(rows) == 3 and all(0 <= r["p_holm"] <= 1 and r["p_holm"] >= r["p_raw"] - 1e-12 for r in rows)


# ---------------------------------------------------------------- data coding (uses the raw export, never modified)
RAW = ROOT / "data" / "raw" / "brca_metabric" / "data_clinical_patient.txt"


@pytest.mark.skipif(not RAW.exists(), reason="raw METABRIC export not present")
def test_metabric_event_coding_and_cohort():
    import phase5_metabric_km as p5
    d = p5.load()
    assert len(d) == 2380 and int(d.event.sum()) == 964          # Phase 1/4 cohort
    assert set(d.event.unique()) == {0, 1}
    assert (d.RFS_MONTHS >= 0).all()
    # coding: '1:Recurred' -> 1, '0:Not Recurred' -> 0
    assert (d.event[d.RFS_STATUS.str.startswith("1")] == 1).all() and (d.event[d.RFS_STATUS.str.startswith("0")] == 0).all()
    # 5-year estimate reproduces the Phase 1 table
    r = km_point(d.RFS_MONTHS.to_numpy(), d.event.to_numpy(), 60)
    assert r["cum_inc"] == pytest.approx(0.2628, abs=5e-4) and r["n_at_risk"] == 1581


@pytest.mark.skipif(not RAW.exists(), reason="raw METABRIC export not present")
def test_subgroup_selection_er():
    import phase5_metabric_km as p5
    d = p5.load(); m = p5.group_masks(d)
    pos, neg = m["ER (IHC) positive"], m["ER (IHC) negative"]
    assert not (pos & neg).any()
    assert int(pos.sum()) == 1735 and int(neg.sum()) == 573


def test_poisson_interval_exact():
    import phase5_metabric_km as p5
    lo, hi = p5.poisson_ci(10, 100.0)
    assert lo == pytest.approx(stats.chi2.ppf(0.025, 20) / 2 / 100)
    assert lo < 0.10 < hi
    assert p5.poisson_ci(0, 50.0)[0] == 0.0 and p5.poisson_ci(0, 50.0)[1] == pytest.approx(3.6889 / 50, rel=1e-3)
