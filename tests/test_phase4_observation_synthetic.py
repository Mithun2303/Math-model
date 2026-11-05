import numpy as np
import pytest

from dormancy.observation import (ObservationAssumptions, exact_censor, generate_dataset, interval_censor, observation_times,
                                  simulate_hidden)
from dormancy.synthetic_validation import (FitProblem, FitResult, fit_from, make_bounds, multistart, prediction_errors,
                                           residuals_jacobian)

G1 = dict(r=0.12, q=0.04, a=0.005, delta_D=0.001, beta=0.04, rho=0.04, eta=0.08, N0=1.01e-4, fD=0.990099, e0=1.0)


def _a(**kw):
    base = dict(observed="N", noise_sigma_ln=0.1, lod=1e-9, schedule_months=12.0, follow_up_months=240.0)
    base.update(kw)
    return ObservationAssumptions(**base)


def test_dataset_reproducible_and_seed_dependent():
    a = _a()
    d1, d2, d3 = generate_dataset(G1, a, 7), generate_dataset(G1, a, 7), generate_dataset(G1, a, 8)
    assert np.array_equal(d1.y["N"], d2.y["N"]) and np.array_equal(d1.times, d2.times)
    assert not np.array_equal(d1.y["N"], d3.y["N"])


def test_observation_truncated_at_detection_and_baseline_visit():
    a = _a()
    d = generate_dataset(G1, a, 1)
    assert 200 < d.T_detect < 240 and d.times[0] == 0.0 and d.times[-1] <= d.T_detect
    full = observation_times(a.with_(truncate_at_detection=False), d.T_detect)
    assert full[-1] > d.T_detect                               # untruncated schedule would continue to the follow-up
    assert observation_times(a.with_(include_baseline=False), d.T_detect)[0] == 12.0
    # no crossing within follow-up -> observed to the follow-up limit
    th = dict(G1, beta=0.15)
    d2 = generate_dataset(th, a, 2)
    assert np.isinf(d2.T_detect) and d2.times[-1] == 240.0


def test_noise_statistics_match_sigma_and_lod_flags():
    sig = 0.2
    a = _a(noise_sigma_ln=sig, schedule_months=1.0, follow_up_months=120.0, truncate_at_detection=False, lod=1e-12)
    d = generate_dataset(G1, a, 3)
    true = simulate_hidden(G1, d.times).N
    z = np.log(d.y["N"]) - np.log(true)
    assert len(z) > 100 and abs(np.std(z) - sig) < 0.05 and abs(np.mean(z)) < 0.05
    assert not d.censored["N"].any()
    hi = generate_dataset(G1, a.with_(lod=3e-4), 3)             # N0 = 1e-4 < LOD: early points left-censored (NaN, flagged)
    c = hi.censored["N"]
    assert c[0] and (~c).any() and np.all(np.isnan(hi.y["N"][c])) and np.all(np.isfinite(hi.y["N"][~c]))


def test_interval_and_exact_censoring():
    v = np.array([6.0, 12.0, 18.0, 24.0])
    assert interval_censor(10.0, v, 24.0) == (6.0, 12.0, 1)
    assert interval_censor(6.0, v, 24.0) == (0.0, 6.0, 1)       # event exactly at a visit: in (0, 6]
    assert interval_censor(24.0, v, 24.0) == (18.0, 24.0, 1)
    L, R, ev = interval_censor(30.0, v, 24.0)
    assert (L, ev) == (24.0, 0) and np.isinf(R)
    L, R, ev = interval_censor(np.inf, v, 24.0)
    assert (L, ev) == (24.0, 0) and np.isinf(R)
    assert exact_censor(50.0, 240.0) == (50.0, 1)
    assert exact_censor(300.0, 240.0) == (240.0, 0) and exact_censor(np.inf, 240.0) == (240.0, 0)


def _problem(free, sigma=1e-4, lod=1e-12, bounds_factor=30.0, **akw):
    a = _a(noise_sigma_ln=sigma, lod=lod, schedule_months=6.0, **akw)
    d = generate_dataset(G1, a, 5)
    truth = dict(G1, c=1.0)
    fixed = {k: truth[k] for k in list(truth) if k not in free}
    lo, hi = make_bounds(truth, free, bounds_factor)
    return FitProblem(d, free, fixed, lo, hi), truth


def test_residual_jacobian_matches_finite_differences():
    free = ["r", "q", "beta", "N0", "fD"]
    prob, truth = _problem(free, sigma=0.1)
    x = prob.to_x(truth) + 0.05
    r, J = residuals_jacobian(prob, x, True)
    Jfd = np.zeros_like(J)
    for j in range(len(x)):
        xp, xm = x.copy(), x.copy()
        xp[j] += 1e-6; xm[j] -= 1e-6
        Jfd[:, j] = (residuals_jacobian(prob, xp, False)[0] - residuals_jacobian(prob, xm, False)[0]) / 2e-6
    assert np.max(np.abs(J - Jfd)) / np.max(np.abs(Jfd)) < 1e-6


def test_tobit_censored_residuals_are_finite():
    free = ["r", "q", "beta"]
    prob, truth = _problem(free, sigma=0.1, lod=3e-4)
    assert prob.data.censored["N"].any()
    r, J = residuals_jacobian(prob, prob.to_x(truth), True)
    assert np.all(np.isfinite(r)) and np.all(np.isfinite(J)) and np.all(r[prob.data.censored["N"]] >= 0)


def test_noise_free_identifiable_subset_is_recovered():
    free = ["beta", "q"]
    prob, truth = _problem(free)
    x0 = prob.to_x(truth) + np.log(1.3)
    fit = fit_from(prob, x0, max_nfev=100)
    est = np.exp(fit.x)
    for nm, v in zip(free, est):
        assert abs(v / truth[nm] - 1) < 0.01, (nm, v)
    assert not fit.bound_hit.any()


def test_bounds_respected_and_bound_hit_flagged():
    free = ["beta"]
    prob, truth = _problem(free)
    prob.lo, prob.hi = np.log([0.01]), np.log([0.02])           # truth (0.04) outside the box
    fit = fit_from(prob, np.log([0.015]), max_nfev=50)
    assert prob.lo[0] - 1e-9 <= fit.x[0] <= prob.hi[0] + 1e-9
    assert fit.bound_hit[0] and abs(fit.x[0] - prob.hi[0]) < 1e-3


def test_multistart_sorted_prediction_keys_and_solver_convergence():
    free = ["beta", "q"]
    prob, truth = _problem(free)
    ms = multistart(prob, 2, seed=1, max_nfev=40)
    assert [m.cost for m in ms] == sorted(m.cost for m in ms)
    pe = prediction_errors(prob, ms[0], truth)
    assert {"rmse_lnN_in_window", "abs_lnN_error_at_120", "abs_lnN_error_at_240", "T_true", "T_pred", "T_rel_error"} <= set(pe)
    x0 = prob.to_x(truth) + np.log(1.3)
    f1 = fit_from(prob, x0, max_nfev=100)
    prob.fit_solver = {"method": "DOP853", "rtol": 1e-11, "atol": 1e-14}
    f2 = fit_from(prob, x0, max_nfev=100)
    assert np.max(np.abs(f1.x - f2.x)) < 1e-3
