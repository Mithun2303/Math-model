from fractions import Fraction

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from dormancy.identifiability import (PARAMS, Dual, GF, PDEModel, forward_sensitivities, identifiability_matrix,
                                      linear_regime_symbolic, nullspace_q, random_point, rank_gf, structural_report, to_gf)

TH = dict(r=0.12, q=0.04, a=0.005, delta_D=0.001, beta=0.04, rho=0.04, eta=0.08, h=0.02, mu=0.05, K=1.0)
X0 = [1e-6, 1e-4, 1.0]


def _rhs(t, y, p):
    P, D, e = y
    N = P + D
    return [p["r"] * P * (1 - N / p["K"]) - p["q"] * P + p["a"] * D - p["beta"] * e * P,
            p["q"] * P - (p["a"] + p["delta_D"]) * D,
            p["mu"] * (1 - e) + p["rho"] * P / (p["h"] + P) * e - p["eta"] * N * e]


def _sim(p, x0, times):
    return solve_ivp(_rhs, (0, max(times)), x0, args=(p,), method="DOP853", rtol=1e-12, atol=1e-16, t_eval=times).y.T


def test_random_point_deterministic_and_positive():
    a, b, c = random_point(3), random_point(3), random_point(4)
    assert a == b and a != c
    assert all(v > 0 for v in a.values())


def test_taylor_coefficients_match_ode_solution():
    pt = random_point(5)
    order = 14
    par = {k: Dual(Fraction(pt[k]), None) for k in PARAMS}
    x0 = [Dual(Fraction(pt[k]), None) for k in ("P0", "D0", "e0")]
    P, D, e = PDEModel().coefficients(order, par, x0, Dual(Fraction(1), None), Dual(Fraction(0), None))
    pf = {k: float(v) for k, v in pt.items()}
    t = 0.3
    ref = solve_ivp(_rhs, (0, t), [pf["P0"], pf["D0"], pf["e0"]], args=(pf,), method="DOP853", rtol=1e-13, atol=1e-16).y[:, -1]
    approx = [sum(float(c.v) * t ** k for k, c in enumerate(S)) for S in (P, D, e)]
    assert np.max(np.abs(np.array(approx) - ref)) < 1e-12


def test_forward_sensitivities_match_finite_differences_all_parameters():
    times = np.linspace(0, 200, 6)
    X, S = forward_sensitivities(TH, X0, PARAMS, times)
    for j, nm in enumerate(PARAMS):
        h = 1e-6 * max(TH[nm], 1e-3)
        pp, pm = dict(TH), dict(TH)
        pp[nm] += h; pm[nm] -= h
        fd = (_sim(pp, X0, times) - _sim(pm, X0, times)) / (2 * h)
        rel = np.max(np.abs(S[:, :, j] - fd)) / (np.max(np.abs(fd)) + 1e-30)
        assert rel < 1e-6, (nm, rel)


def test_initial_state_sensitivity_matches_finite_differences():
    times = np.linspace(0, 150, 5)
    # sensitivity w.r.t. ln N0 with fD fixed: dx0/dlnN0 = (P0, D0, 0)
    X, S = forward_sensitivities(TH, X0, ["N0"], times, ic_sens={"N0": np.array([X0[0], X0[1], 0.0])})
    h = 1e-6
    up = _sim(TH, [X0[0] * np.exp(h), X0[1] * np.exp(h), 1.0], times)
    dn = _sim(TH, [X0[0] * np.exp(-h), X0[1] * np.exp(-h), 1.0], times)
    fd = (up - dn) / (2 * h)
    assert np.max(np.abs(S[:, :, 0] - fd)) / np.max(np.abs(fd)) < 1e-6


def test_full_model_generically_identifiable_from_total_burden():
    U = ["r", "q", "a", "delta_D", "beta", "rho", "eta", "P0", "D0"]
    rep = structural_report(U, "N", max_order=22, n_points=2)
    assert rep["rank"] == 9 and rep["flat_last_orders"] and rep["structurally_locally_identifiable"]


def test_linear_submodel_has_documented_deficiencies():
    lin = PDEModel(logistic=False, immune=False)
    assert structural_report(["r", "q", "a", "delta_D", "beta", "P0", "D0"], "N", model=lin, max_order=20, n_points=2)["rank"] == 4
    assert structural_report(["r", "beta"], "N", model=lin, max_order=20, n_points=2)["rank"] == 1   # only r - beta enters
    assert structural_report(["q", "a", "delta_D"], "N", model=lin, max_order=20, n_points=2)["rank"] == 3
    assert structural_report(["r", "q", "a", "delta_D", "beta", "P0", "D0"], "PD", model=lin, max_order=20, n_points=2)["rank"] == 6


def test_scale_gauge_when_burden_scale_and_signal_scale_unknown():
    base = ["r", "q", "a", "delta_D", "beta", "rho", "eta", "P0", "D0", "h", "mu", "e0", "c"]
    assert structural_report(base, "N", max_order=26, n_points=1)["deficiency"] == 0               # c alone: identifiable
    assert structural_report(base + ["K"], "N", max_order=26, n_points=1)["deficiency"] == 1       # c and K: scale gauge


def test_linear_regime_symbolic_rank():
    s = linear_regime_symbolic()
    assert s["rank"] == 4 and s["null_dim"] == 2


def test_rank_gf_on_explicit_matrices_and_toy_models():
    assert rank_gf([[GF(1), GF(2)], [GF(2), GF(4)]]) == 1
    assert rank_gf([[GF(1), GF(0)], [GF(0), GF(1)]]) == 2

    def series_rows(unknown_idx, vals, step_fn, order=8):
        """rows = gradients of Taylor coefficients of x(t) for a scalar toy ODE over GF(p); vals = (p1, p2, x0)."""
        n = 3
        mk = lambda i: Dual(to_gf(Fraction(vals[i])), [GF(1) if j == i else GF(0) for j in range(n)])
        pars = [mk(0), mk(1), mk(2)]
        x = step_fn(pars, order)
        return [list(c.g) for c in x]

    def decay(pars, order):                # x' = -a b x : only the product a*b is identifiable
        a, b, x0 = pars
        x = [x0]
        for k in range(order):
            x.append(x[k] * (a * b) * (-1) * to_gf(Fraction(1, k + 1)))
        return x

    def logistic(pars, order):             # x' = r x (1 - x/K): (r, K, x0) all identifiable
        r, K, x0 = pars
        x = [x0]
        for k in range(order):
            conv = x[0] * x[k]
            for i in range(1, k + 1):
                conv = conv + x[i] * x[k - i]
            x.append((r * x[k] - r * K.inv() * conv) * to_gf(Fraction(1, k + 1)))
        return x

    assert rank_gf(series_rows(None, (Fraction(3, 7), Fraction(5, 11), Fraction(2, 3)), decay)) == 2       # 3 unknowns, rank 2
    assert rank_gf(series_rows(None, (Fraction(3, 7), Fraction(5, 11), Fraction(2, 3)), logistic)) == 3


def test_gf_rank_agrees_with_exact_rational_rank():
    lin = PDEModel(logistic=False, immune=False)
    U = ["q", "a", "delta_D", "P0", "D0"]
    pt = random_point(11)
    rows_q = identifiability_matrix(lin, U, pt, 10, "N", "q")
    rows_gf = identifiability_matrix(lin, U, pt, 10, "N", "gf")
    vecs, rk = nullspace_q(rows_q)
    assert rk == rank_gf(rows_gf) == 4 and len(vecs) == 1
