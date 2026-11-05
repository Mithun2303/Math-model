import math

import numpy as np
import pytest
from scipy import stats

from dormancy import NormalizedParams, simulate
from dormancy.uncertainty import (Prior, build_global_prior, build_prior, clopper_pearson, default_options, evaluate_row,
                                  event_summary, load_phase3_config, parallel_map, regime_of, row_to_model, run_trajectory,
                                  upper_bound_rare, wilson_interval)

CFG = load_phase3_config()
BASE = CFG["baseline"]


def base_row(**kw):
    r = dict(BASE)
    r.update(kw)
    return r


# ----------------------------------------------------------------------------------------- sampling
def test_sampling_is_reproducible_and_seed_dependent():
    pr = build_prior(CFG, "P1_base")
    for method in ("mc", "lhs", "sobol"):
        a, b, c = pr.sample(64, 11, method), pr.sample(64, 11, method), pr.sample(64, 12, method)
        assert a == b
        assert a != c


def test_mc_prefix_property_nested_runs():
    pr = build_prior(CFG, "P1_base")
    assert pr.sample(1000, 5, "mc") == pr.sample(5000, 5, "mc")[:1000]        # the 1000-run is a prefix of the 5000-run


def test_lhs_is_stratified():
    """Each marginal's CDF value u, binned into n equal strata, must contain exactly one point per stratum."""
    pr = build_global_prior(CFG)
    n = 200
    rows = pr.sample(n, 3, "lhs")
    for name in pr.names:
        spec = CFG["global_ranges"][name]
        x = np.array([r[name] for r in rows])
        lo, hi = spec["lo"], spec["hi"]
        if spec["scale"] == "log":
            u = (np.log(x) - np.log(lo)) / (np.log(hi) - np.log(lo))
        elif spec["scale"] == "logit":
            L = lambda v: np.log(v / (1 - v))
            u = (L(x) - L(lo)) / (L(hi) - L(lo))
        else:
            u = (x - lo) / (hi - lo)
        assert sorted(np.floor(np.clip(u, 0, 1 - 1e-12) * n).astype(int)) == list(range(n)), name
    # and a plain Monte Carlo sample is NOT stratified (the test has power)
    mc = pr.sample(n, 3, "mc")
    x = np.array([r["r"] for r in mc]); lo, hi = CFG["global_ranges"]["r"]["lo"], CFG["global_ranges"]["r"]["hi"]
    u = (np.log(x) - np.log(lo)) / (np.log(hi) - np.log(lo))
    assert sorted(np.floor(u * n).astype(int)) != list(range(n))


def test_samples_respect_ranges_and_positivity():
    g = CFG["global_ranges"]
    rows = build_global_prior(CFG).sample(500, 1, "lhs")
    for r in rows:
        for n, spec in g.items():
            assert spec["lo"] * (1 - 1e-9) <= r[n] <= spec["hi"] * (1 + 1e-9), n
        assert 0 < r["fD"] < 1 and r["N0"] < min(CFG["thresholds"])          # initial burden strictly below every threshold
    assert max(g["N0"]["hi"], 0) < g["N_detect"]["lo"]


def test_parameter_range_validation():
    ok = {"q": {"family": "loguniform", "lo": 0.01, "hi": 0.1}}
    Prior(ok)
    with pytest.raises(ValueError):
        Prior({"q": {"family": "loguniform", "lo": 0.1, "hi": 0.01}})
    with pytest.raises(ValueError):
        Prior({"q": {"family": "loguniform", "lo": 0.0, "hi": 1.0}})
    with pytest.raises(ValueError):
        Prior({"q": {"family": "lognormal", "median": -1.0, "factor95": 2.0}})
    with pytest.raises(ValueError):
        Prior({"q": {"family": "lognormal", "median": 1.0, "factor95": 1.0}})
    with pytest.raises(ValueError):
        Prior({"fD": {"family": "logit", "lo": 0.5, "hi": 1.2}})
    with pytest.raises(ValueError):
        Prior({"fD": {"family": "beta", "alpha": 0.0, "beta": 1.0}})
    with pytest.raises(ValueError):
        Prior({"nonsense": {"family": "uniform", "lo": 0, "hi": 1}})
    with pytest.raises(ValueError):
        Prior({"q": ok["q"], "a": {"family": "loguniform", "lo": 1e-3, "hi": 1e-1}}, copula=[["q", "a", 1.2]])
    with pytest.raises(ValueError):
        Prior({"q": ok["q"]}, copula=[["q", "a", 0.3]])                       # copula variable not sampled
    # config sanity
    for n, spec in CFG["global_ranges"].items():
        assert spec["lo"] < spec["hi"]


def test_lognormal_factor95_and_truncation():
    pr = Prior({"q": {"family": "lognormal", "median": 0.04, "factor95": 3.0}})
    q = np.array([r["q"] for r in pr.sample(200000, 1, "mc")])
    assert np.median(q) == pytest.approx(0.04, rel=0.02)
    assert np.quantile(q, 0.975) == pytest.approx(0.04 * 3.0, rel=0.03)
    assert np.quantile(q, 0.025) == pytest.approx(0.04 / 3.0, rel=0.03)
    prt = Prior({"N0": {"family": "lognormal", "median": 1e-4, "factor95": 30.0, "lo": 1e-5, "hi": 1e-3}})
    n0 = np.array([r["N0"] for r in prt.sample(20000, 2, "mc")])
    assert n0.min() >= 1e-5 * (1 - 1e-9) and n0.max() <= 1e-3 * (1 + 1e-9)


def test_gaussian_copula_induces_requested_rank_correlation():
    pr = build_prior(CFG, "P4_corr_pos")
    rows = pr.sample(20000, 4, "mc")
    rho_qa = stats.spearmanr([r["q"] for r in rows], [r["a"] for r in rows])[0]
    rho_br = stats.spearmanr([r["beta"] for r in rows], [r["rho"] for r in rows])[0]
    # Gaussian-copula Spearman = (6/pi) asin(r/2) = 0.4826 for r = 0.5
    assert rho_qa == pytest.approx(6 / math.pi * math.asin(0.25), abs=0.03)
    assert rho_br == pytest.approx(6 / math.pi * math.asin(0.25), abs=0.03)
    ind = build_prior(CFG, "P1_base").sample(20000, 4, "mc")
    assert abs(stats.spearmanr([r["q"] for r in ind], [r["a"] for r in ind])[0]) < 0.03


def test_all_configured_priors_build_and_sample_valid_models():
    for name in CFG["priors"]:
        pr = build_prior(CFG, name)
        for r in pr.sample(30, 1, "mc"):
            p, y0, thr = row_to_model(r, CFG["fixed"])
            assert np.all(y0 >= 0) and y0[0] + y0[1] < thr
            assert p.E_free == pytest.approx(1.0)                              # normalised gauge: tumour-free level e = 1


# ----------------------------------------------------------------------------------------- probabilities / censoring
@pytest.mark.parametrize("k,n", [(0, 10), (1, 10), (5, 10), (10, 10), (37, 1000), (0, 5000)])
def test_interval_bounds_are_valid_probabilities(k, n):
    for f in (wilson_interval, clopper_pearson):
        lo, hi = f(k, n)
        assert 0.0 <= lo <= k / n <= hi <= 1.0
    assert upper_bound_rare(k, n) >= k / n


def test_rare_event_upper_bound_rule_of_three():
    n = 5000
    assert upper_bound_rare(0, n) == pytest.approx(1 - 0.05 ** (1 / n), rel=1e-9)
    assert upper_bound_rare(0, n) == pytest.approx(3 / n, rel=0.01)


def test_event_summary_censoring_and_medians():
    T = np.array([10, 50, 100, 200, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan], float)
    ev = np.array([1, 1, 1, 1, 0, 0, 0, 0, 0, 0], float)
    s = event_summary(T, ev, (60, 120, 240), 240.0)
    assert s["P_by_60"] == 0.2 and s["P_by_120"] == 0.3 and s["P_by_240"] == 0.4
    assert s["RMST_240"] == pytest.approx((10 + 50 + 100 + 200 + 6 * 240) / 10)    # censored at the horizon contribute 240
    assert s["RMST_60"] == pytest.approx((10 + 50 + 2 * 60 + 6 * 60) / 10)   # events at 100, 200 and all censored count as 60
    assert s["conditional_median_T"] == 75.0
    assert math.isnan(s["unconditional_median_T"])                                 # P(event) < 0.5: no unconditional median
    ev2 = np.ones(4); T2 = np.array([1, 2, 3, 4.0])
    assert event_summary(T2, ev2, (60,), 240.0)["unconditional_median_T"] == 2.5
    # start-above samples (NaN event) are excluded, never counted as events
    s3 = event_summary(np.array([5.0, np.nan, np.nan]), np.array([1.0, np.nan, 0.0]), (60,), 240.0)
    assert s3["n"] == 2 and s3["P_by_60"] == 0.5


# ----------------------------------------------------------------------------------------- first crossing and exceedance
@pytest.mark.parametrize("beta", [0.02, 0.04, 0.0567, 0.1])
def test_dense_runner_matches_verified_solver(beta):
    p, y0, thr = row_to_model(base_row(beta=beta), CFG["fixed"])
    t = run_trajectory(p, y0, [0.1], 240.0, main_threshold=0.1)
    s = simulate(p, y0, 0.1, 240.0)
    assert t["event_0p1"] == s.event
    if s.event:
        assert t["T_0p1"] == pytest.approx(s.time, abs=1e-5)
    else:
        assert math.isnan(t["T_0p1"])                                          # no-event trajectory: no time, never 0


def test_no_event_trajectory_is_censored_not_zero():
    p, y0, _ = row_to_model(base_row(beta=0.2), CFG["fixed"])                  # beta > beta_c: controlled
    t = run_trajectory(p, y0, [0.1], 240.0, main_threshold=0.1)
    assert t["event_0p1"] == 0 and math.isnan(t["T_0p1"]) and t["N_240"] < 1e-3 and t["solver_ok"]


def test_start_above_threshold_is_not_an_event():
    p, y0, _ = row_to_model(base_row(N0=0.2), CFG["fixed"])
    assert simulate(p, y0, 0.1, 240.0).time == 0.0                               # the legacy convention
    t = run_trajectory(p, y0, [0.1, 0.5], 240.0, main_threshold=0.1)
    assert t["start_above_0p1"] and math.isnan(t["event_0p1"])                   # excluded in Phase 3
    assert t["event_0p5"] in (0, 1) and not t["start_above_0p5"]
    d = evaluate_row((base_row(N0=0.2), default_options(CFG, thresholds=[])))
    assert d["start_above_main"] and math.isnan(d["event_main"])


def test_crossing_times_increase_with_threshold():
    p, y0, _ = row_to_model(base_row(beta=0.02), CFG["fixed"])
    ths = [0.01, 0.02, 0.05, 0.1, 0.2]
    t = run_trajectory(p, y0, ths, 240.0)
    Ts = [t[f"T_{th:g}".replace(".", "p")] for th in ths]
    assert all(np.isfinite(Ts)) and Ts == sorted(Ts)


def test_persistent_exceedance_at_baseline():
    p, y0, _ = row_to_model(base_row(), CFG["fixed"])
    t = run_trajectory(p, y0, [0.1], 240.0, main_threshold=0.1)
    assert t["event_0p1"] == 1 and t["persistent_0p1"] == 1 and t["transient_0p1"] == 0 and math.isnan(t["T_down_0p1"])


def test_transient_exceedance_fixture_crosses_then_falls_back():
    # found numerically (Phase 3): strong immune stimulation, two stable nonzero states, N0 = 0.033 starts below 0.05,
    # rises above it, and relaxes to the LOW attractor (N* = 0.0028): first crossing is registered, exceedance is transient
    row = dict(r=0.12, q=0.0015841984, a=0.3062386609, delta_D=0.001, beta=0.06, rho=0.2404914173, eta=1.6373830589)
    p = NormalizedParams(**row).to_params()
    f, N0, e0 = 0.4930310852709215, 0.03280264166486928, 0.6059047928963627
    y0 = [N0 * (1 - f), N0 * f, e0 * p.E_free]
    t = run_trajectory(p, y0, [0.05], 240.0, main_threshold=0.05)
    assert t["event_0p05"] == 1 and t["transient_0p05"] == 1 and t["persistent_0p05"] == 0
    assert t["T_0p05"] == pytest.approx(8.7498, abs=1e-3) and t["T_down_0p05"] == pytest.approx(18.681, abs=1e-3)
    assert t["T_down_0p05"] > t["T_0p05"] and t["N_240"] < 0.05 and t["N_max"] > 0.05
    s = simulate(p, np.array(y0), 0.05, 240.0)
    assert s.event == 1 and s.time == pytest.approx(t["T_0p05"], abs=1e-5)        # the verified solver registers the same first crossing


def test_regime_classification_of_known_cases():
    base = dict(r=0.12, q=0.04, a=0.005, delta_D=0.001, beta=0.04, rho=0.04, eta=0.08)
    assert regime_of(NormalizedParams(**base).to_params())["regime"] == "NZ only"
    assert regime_of(NormalizedParams(**dict(base, beta=0.2)).to_params())["regime"] == "TF only"
    assert regime_of(NormalizedParams(**dict(base, eta=0.6, beta=0.4 * 0.572083602)).to_params())["regime"] == "TF + NZ"
    assert regime_of(NormalizedParams(**dict(base, eta=0.164, beta=0.4 * 0.27)).to_params())["regime"] == "NZ + NZ"


# ----------------------------------------------------------------------------------------- numerical robustness / convergence
def test_solver_tolerances_do_not_change_events():
    rows = build_prior(CFG, "P1_base").sample(40, 8, "mc")
    base = [evaluate_row((r, default_options(CFG, thresholds=[]))) for r in rows]
    for kw in (dict(rtol=1e-6, atol=1e-9), dict(rtol=1e-12, atol=1e-15), dict(method="Radau")):
        alt = [evaluate_row((r, default_options(CFG, thresholds=[], **kw))) for r in rows]
        assert [a["event_main"] for a in base] == [a["event_main"] for a in alt]
        for a, b in zip(base, alt):
            if a["event_main"] == 1:
                assert a["T_main"] == pytest.approx(b["T_main"], abs=0.05)


def test_trajectories_stay_nonnegative_and_finite_over_wide_prior():
    for r in build_global_prior(CFG).sample(60, 2, "lhs"):
        out = evaluate_row((r, default_options(CFG, thresholds=[])))
        assert out["solver_ok"] and out["min_state"] > -1e-10 and np.isfinite(out["N_240"])


def test_monte_carlo_estimate_converges_and_parallel_equals_serial():
    pr = build_prior(CFG, "P1_base")
    opt = default_options(CFG, thresholds=[], check_against_simulate=False)
    rows = pr.sample(1200, CFG["seeds"]["mc_base"], "mc")
    par = parallel_map(evaluate_row, [(r, opt) for r in rows], workers=4)
    ser = parallel_map(evaluate_row, [(r, opt) for r in rows[:60]], workers=1)
    assert [a["event_main"] for a in par[:60]] == [a["event_main"] for a in ser]    # scheduling does not matter
    ev = np.array([a["event_main"] for a in par]); T = np.array([a["T_main"] for a in par])
    big = event_summary(T, ev, (240,), 240.0)["P_by_240"]
    widths = []
    for n in (150, 600, 1200):
        s = event_summary(T[:n], ev[:n], (240,), 240.0)
        widths.append(s["P_by_240_wilson_hi"] - s["P_by_240_wilson_lo"])
        assert s["P_by_240_wilson_lo"] <= big <= s["P_by_240_wilson_hi"] + 0.02
    assert widths[0] > widths[1] > widths[2]
    assert widths[0] / widths[2] == pytest.approx(math.sqrt(1200 / 150), rel=0.25)    # 1/sqrt(n) scaling
