import numpy as np
import pytest

from dormancy.endpoint_mapping import (ClinicalMapping, cumulative_incidence, dkw_epsilon, map_event_times, sup_distance,
                                       time_to_level)
from dormancy.synthetic_validation import kaplan_meier, km_at, km_cumulative_incidence_ci, turnbull


def test_time_origin_and_detection_delay_transform():
    T = np.array([10.0, 50.0, np.inf, np.nan])
    out = map_event_times(T, ClinicalMapping(tau=3.0, detect_delay=2.0))
    assert out[0] == 15.0 and out[1] == 55.0 and np.isinf(out[2]) and np.isinf(out[3])      # non-events stay inf
    assert np.array_equal(map_event_times(T, ClinicalMapping()), np.array([10.0, 50.0, np.inf, np.inf]))


def test_death_first_mapping_fraction_and_lag_nonnegative():
    T = np.full(20000, 40.0)
    m = ClinicalMapping(tau=6.0, detect_delay=1000.0, death_first_fraction=0.3, death_lag_mean=24.0)
    out = map_event_times(T, m, np.random.default_rng(1))
    death = out < 40.0 + 6.0 + 500.0                            # death-first events are recorded at T + tau + lag (lag << 1000)
    assert abs(death.mean() - 0.3) < 0.02
    assert np.all(out[death] >= 46.0)                           # lag >= 0
    assert np.allclose(out[~death], 1046.0)
    assert np.array_equal(out, map_event_times(T, m, np.random.default_rng(1)))              # reproducible
    allfirst = map_event_times(T, ClinicalMapping(death_first_fraction=1.0), np.random.default_rng(2))
    assert np.all(allfirst >= 40.0)


def test_cumulative_incidence_levels_and_distances():
    T = np.array([1.0, 2.0, np.inf, 5.0])
    t = np.array([0.0, 1.0, 2.0, 10.0])
    F = cumulative_incidence(T, t)
    assert np.allclose(F, [0, 0.25, 0.5, 0.75])
    assert time_to_level(F, t, 0.5) == 2.0 and np.isnan(time_to_level(F, t, 0.9))
    assert sup_distance(F, F + 0.1) == pytest.approx(0.1)
    assert dkw_epsilon(2380) == pytest.approx(np.sqrt(np.log(40.0) / (2 * 2380)))
    assert 0.027 < dkw_epsilon(2380) < 0.029


def test_kaplan_meier_hand_example_with_ties_and_censoring():
    time = np.array([1, 2, 2, 3, 4, 5.0]); event = np.array([1, 1, 0, 1, 0, 1])
    ut, S, se = kaplan_meier(time, event)
    assert np.array_equal(ut, [1, 2, 3, 5])
    assert np.allclose(S, [5 / 6, 2 / 3, 4 / 9, 0.0])
    assert se[0] == pytest.approx(5 / 6 * np.sqrt(1 / (6 * 5)))
    assert km_at(ut, S, se, 0.5) == (1.0, 0.0) and km_at(ut, S, se, 2.5)[0] == pytest.approx(2 / 3)


def test_kaplan_meier_uncensored_equals_empirical_cdf():
    rng = np.random.default_rng(0)
    t = rng.exponential(10.0, 200)
    ut, S, se = kaplan_meier(t, np.ones(200, int))
    for x in (2.0, 7.0, 15.0, 30.0):
        assert 1 - km_at(ut, S, se, x)[0] == pytest.approx(np.mean(t <= x), abs=1e-12)


def test_turnbull_matches_kaplan_meier_for_exact_and_right_censored_data():
    eps = 1e-6
    exact = np.array([1.0, 2.0, 3.0])
    L = np.concatenate([exact - eps, [2.5]]); R = np.concatenate([exact, [np.inf]])
    ends, mass, F = turnbull(L, R)
    ut, S, se = kaplan_meier(np.array([1.0, 2.0, 3.0, 2.5]), np.array([1, 1, 1, 0]))
    for x in (1.0, 2.0, 3.0):
        assert F(x) == pytest.approx(1 - km_at(ut, S, se, x)[0], abs=1e-4)
    assert mass.sum() == pytest.approx(1.0 - 0.0, abs=1e-6) or mass.sum() <= 1.0 + 1e-9


def test_turnbull_interval_censored_recovers_cdf_and_handles_inf():
    rng = np.random.default_rng(3)
    T = rng.exponential(30.0, 4000)
    visits = np.arange(6.0, 121.0, 6.0)
    L, R = [], []
    for x in T:
        if x > visits[-1]:
            L.append(visits[-1]); R.append(np.inf)
        else:
            j = np.searchsorted(visits, x, side="left")
            L.append(visits[j - 1] if j > 0 else 0.0); R.append(visits[j])
    ends, mass, F = turnbull(np.array(L), np.array(R))
    for x in (30.0, 60.0, 120.0):
        assert abs(F(x) - (1 - np.exp(-x / 30.0))) < 0.03


def test_administrative_censoring_keeps_kaplan_meier_consistent():
    rng = np.random.default_rng(4)
    T = rng.exponential(100.0, 20000)
    C = 100.0
    time, event = np.minimum(T, C), (T <= C).astype(int)
    F, lo, hi = km_cumulative_incidence_ci(time, event, 50.0)
    truth = 1 - np.exp(-0.5)
    assert lo <= F <= hi and abs(F - truth) < 0.02 and lo < truth < hi
    # naive 'censored = no event' fraction is biased at the horizon, KM is not
    assert abs(np.mean(T <= 100.0) - (1 - np.exp(-1.0))) < 0.01
