import math

import numpy as np
import pytest

from dormancy import NormalizedParams
from dormancy.basins import bisect_boundary, boundary_mask, classify_ic, tight_options
from dormancy.sensitivity import (bootstrap_ci, flip_margin, ishigami, log_derivative, morris_design, morris_effects,
                                  perturb_row, prcc, saltelli_unit_matrices, sobol_from_outputs, spearman_columns)
from dormancy.uncertainty import default_options, load_phase3_config, regime_of, run_trajectory, row_to_model

CFG = load_phase3_config()
BASE = dict(CFG["baseline"])


def test_sobol_estimators_reproduce_ishigami_analytic_indices():
    d, N = 3, 8192
    A, B, AB, BA = saltelli_unit_matrices(d, N, seed=7)
    tr = lambda U: -np.pi + 2 * np.pi * U
    f = lambda M: ishigami(tr(M))
    res = sobol_from_outputs(f(A), f(B), np.array([f(m) for m in AB]), np.array([f(m) for m in BA]))
    assert res["S1"] == pytest.approx([0.3139, 0.4424, 0.0], abs=0.04)
    assert res["ST"] == pytest.approx([0.5576, 0.4424, 0.2437], abs=0.04)
    assert res["S2"][0, 2] == pytest.approx(0.2437, abs=0.05)       # the x1-x3 interaction
    assert abs(res["S2"][0, 1]) < 0.05 and abs(res["S2"][1, 2]) < 0.05


def test_sobol_additive_model_has_no_interaction():
    d, N = 3, 4096
    A, B, AB, BA = saltelli_unit_matrices(d, N, seed=3)
    f = lambda M: 1.0 * M[:, 0] + 2.0 * M[:, 1] + 0.0 * M[:, 2]
    res = sobol_from_outputs(f(A), f(B), np.array([f(m) for m in AB]), np.array([f(m) for m in BA]))
    expect = np.array([1.0, 4.0, 0.0]) / 5.0
    assert res["S1"] == pytest.approx(expect, abs=0.03) and res["ST"] == pytest.approx(expect, abs=0.03)
    assert abs(res["S2"][0, 1]) < 0.03


def test_log_derivative_exact_for_power_law_and_flags_problems():
    th0 = 2.0
    mk = lambda f, h: {s: f(th0 * math.exp(s)) for s in (h, -h, h / 2, -h / 2)}
    d, st = log_derivative(mk(lambda t: t ** 3, 0.05), 0.05, log_output=True)
    assert st == "ok" and d == pytest.approx(3.0, rel=1e-3)
    d, st = log_derivative(mk(lambda t: 5.0, 0.05), 0.05, log_output=True)
    assert st == "zero" and d == 0.0
    vals = mk(lambda t: 1.0, 0.05); vals[-0.05] = float("nan")
    assert log_derivative(vals, 0.05, log_output=False)[1] == "undefined"             # missing stencil output (e.g. no event)
    step = {0.05: 10.0, -0.05: 1.0, 0.025: 10.0, -0.025: 1.0}                          # a discontinuity inside the stencil
    assert log_derivative(step, 0.05, log_output=False)[1] in ("unstable", "ok")
    jump = {0.05: 10.0, -0.05: 1.0, 0.025: 1.0, -0.025: 1.0}
    assert log_derivative(jump, 0.05, log_output=False)[1] == "unstable"               # the two step sizes disagree


def test_event_time_derivative_not_computed_across_event_boundary():
    # baseline reference: event at ~220 months; a 5% decrease of r removes the event within 240 months
    row = dict(BASE)
    T = {}
    for s in (0.05, -0.05, 0.025, -0.025):
        p, y0, thr = row_to_model(perturb_row(row, "r", s), CFG["fixed"])
        t = run_trajectory(p, y0, [0.1], 240.0, main_threshold=0.1)
        T[s] = t["T_0p1"]
    assert not all(np.isfinite(list(T.values())))                                      # at least one stencil point has no event
    assert log_derivative(T, 0.05, log_output=True)[1] == "undefined"


def test_flip_margin_synthetic_and_perturb_row():
    ind = lambda r: float(r["x"] > 1.7)
    fm = flip_margin(ind, {"x": 1.0}, "x", 0, max_log=math.log(10.0), n_scan=80)
    assert fm["up_factor"] == pytest.approx(1.7, rel=1e-6) and math.isinf(fm["down_factor"])
    r = perturb_row(BASE, "P0", math.log(2.0))
    P0, D0 = BASE["N0"] * (1 - BASE["fD"]), BASE["N0"] * BASE["fD"]
    assert r["N0"] * (1 - r["fD"]) == pytest.approx(2 * P0) and r["N0"] * r["fD"] == pytest.approx(D0)
    r = perturb_row(BASE, "D0", math.log(0.5))
    assert r["N0"] * r["fD"] == pytest.approx(0.5 * D0) and r["N0"] * (1 - r["fD"]) == pytest.approx(P0)


def test_spearman_prcc_and_bootstrap():
    rng = np.random.default_rng(0)
    n = 4000
    x1 = rng.normal(size=n); x2 = x1 + 0.3 * rng.normal(size=n); x3 = rng.normal(size=n)
    y = x1 + 0.2 * rng.normal(size=n)
    X = np.column_stack([x1, x2, x3])
    rho = spearman_columns(X, y)
    pr = prcc(X, y)
    assert rho[1] > 0.7 and abs(pr[1]) < 0.1 and pr[0] > 0.8 and abs(pr[2]) < 0.05    # PRCC removes the confounded input
    lo, hi = bootstrap_ci(spearman_columns, X, y, B=100, seed=1)
    assert np.all(lo <= rho + 1e-9) and np.all(rho <= hi + 1e-9)
    # monotone transform invariance of rank statistics
    assert spearman_columns(X, np.exp(y)) == pytest.approx(rho)


def test_morris_elementary_effects_on_linear_model():
    d, r = 4, 30
    X, orders, signs = morris_design(d, r, 4, seed=5)
    assert X.shape == (r, d + 1, d) and X.min() >= 0 and X.max() <= 1
    coef = np.array([4.0, 2.0, 1.0, 0.0])
    Y = X @ coef
    res = morris_effects(Y, orders, signs, 4 / 6)
    assert res["mu_star"] == pytest.approx(np.abs(coef), abs=1e-9) and np.all(res["sigma"] < 1e-9)


def test_local_model_sensitivity_signs():
    # more killing -> smaller burden; faster growth -> larger burden; N_detect cannot change the trajectory
    def N240(**kw):
        row = dict(BASE); row.update(kw)
        p, y0, thr = row_to_model(row, CFG["fixed"])
        return run_trajectory(p, y0, [0.1], 240.0)["N_240"]
    n0 = N240()
    assert N240(beta=BASE["beta"] * 1.05) < n0 < N240(beta=BASE["beta"] * 0.95)
    assert N240(r=BASE["r"] * 0.95) < n0 < N240(r=BASE["r"] * 1.05)
    assert N240(N_detect=0.3) == pytest.approx(n0, rel=1e-12)


# ----------------------------------------------------------------------------------------- basins
def test_boundary_mask_marks_label_changes_only():
    L = np.array([[0, 0, 1], [0, 0, 1], [0, 0, 0]])
    m = boundary_mask(L)
    # cells with a 4-neighbour of different label: (0,1) (0,2) (1,1) (1,2) (2,2); (2,1) has equal-label neighbours only
    assert m.tolist() == [[False, True, True], [False, True, True], [False, False, True]]
    assert not boundary_mask(np.zeros((4, 4), int)).any()


def test_basin_classification_reproduces_phase2_bistability():
    base = dict(r=0.12, q=0.04, a=0.005, delta_D=0.001, beta=0.4 * 0.5720836024531797, rho=0.04, eta=0.6)
    p = NormalizedParams(**base).to_params()
    reg = regime_of(p)
    assert reg["regime"] == "TF + NZ"
    opts = default_options(CFG, long_floor=1e5)
    small = classify_ic((p, [1e-6, 1e-4, p.E_free], 0.1, opts, reg))          # Phase 2: small residual disease -> tumour-free
    big = classify_ic((p, [0.3, 0.3, 0.1 * p.E_free], 0.1, opts, reg))        # Phase 2: large burden -> nonzero equilibrium
    assert small["attractor"] == "TF" and big["attractor"] == "NZ"
    assert small["event240"] == 0 and np.isnan(small["T240"])                  # no event
    assert np.isnan(big["event240"]) and big["start_above"]                    # starts above 0.1: not a relapse event
    assert big["N_inf"] == pytest.approx(reg["stable_N"][0], rel=1e-3)
    # tighter numerics reproduce the labels
    t = tight_options(opts)
    assert classify_ic((p, [1e-6, 1e-4, p.E_free], 0.1, t, reg))["attractor"] == "TF"
    assert classify_ic((p, [0.3, 0.3, 0.1 * p.E_free], 0.1, t, reg))["attractor"] == "NZ"
    # a locally stable tumour-free state does not guarantee control: the NZ basin is non-empty
    b = bisect_boundary(p, [1e-4 * 0.13, 1e-4 * 0.87, p.E_free], [0.5 * 0.13, 0.5 * 0.87, p.E_free], 0.1, opts, n_iter=25)
    assert b["found"] and 0 < b["t_star"] < 1 and {b["label_a"], b["label_b"]} == {"TF", "NZ"}


def test_two_nonzero_basins_and_unresolved_flag():
    base = dict(r=0.12, q=0.0015841984, a=0.3062386609, delta_D=0.001, beta=0.06, rho=0.2404914173, eta=1.6373830589)
    p = NormalizedParams(**base).to_params()
    reg = regime_of(p)
    assert reg["regime"] == "NZ + NZ" and not reg["tf_stable"]
    opts = default_options(CFG, long_floor=1e5)
    lo = classify_ic((p, [1e-4 * 0.5, 1e-4 * 0.5, p.E_free], 0.1, opts, reg))
    hi = classify_ic((p, [0.5, 0.4, p.E_free], 0.1, opts, reg))
    assert lo["attractor"] == "NZ-low" and hi["attractor"] == "NZ-high"
    assert lo["N_inf"] == pytest.approx(reg["stable_N"][0], rel=1e-3) and hi["N_inf"] == pytest.approx(reg["stable_N"][1], rel=1e-3)
