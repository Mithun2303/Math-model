import numpy as np
import pytest
import sympy as sp

from dormancy import Params, critical_k, jacobian, tumour_free
from dormancy.bifurcation import (degenerate_eta, degenerate_rho, extended_root_near_zero, kprime0, scan_branch,
                                  transcritical_report)
from dormancy.continuation import continue_branch, trace_all_branches
from dormancy.equilibria import branch_k, enumerate_equilibria, fold_P
from dormancy.model import hessian_form, jacobian_k, rhs


def test_hessian_matches_sympy():
    P, D, E = sp.symbols("P D E")
    p = Params(k=0.2)
    F = sp.Matrix([
        p.r * P * (1 - (P + D) / p.K) - p.q * P + p.a * D - p.k * E * P,
        p.q * P - (p.a + p.delta_D) * D,
        p.s + p.rho * P / (p.h + P) * E - p.mu * E - p.eta * (P + D) * E])
    y = np.array([0.03, 0.2, 0.3]); w = np.array([0.7, -0.4, 1.3])
    t = sp.symbols("t")
    Fy = F.subs({P: y[0] + t * w[0], D: y[1] + t * w[1], E: y[2] + t * w[2]})
    d2 = np.array([float(sp.diff(f, t, 2).subs(t, 0)) for f in Fy])
    assert np.allclose(hessian_form(y, w, p), d2, atol=1e-12)


def test_jacobian_k_matches_finite_difference():
    p = Params(k=0.2); y = np.array([0.03, 0.2, 0.3]); h = 1e-6
    fd = (jacobian(y, p.with_(k=0.2 + h)) - jacobian(y, p.with_(k=0.2 - h))) / (2 * h)
    assert np.allclose(jacobian_k(y, p), fd, atol=1e-8)


def test_sotomayor_conditions_at_baseline():
    r = transcritical_report(Params())
    assert r.simple_zero_eigenvalue
    assert abs(r.F_k_projection) < 1e-14                 # F_k = 0 at the tumour-free state
    assert r.alpha < -1e-3                               # transversality: eigenvalue decreases through 0 as k grows
    assert r.beta < -1e-3                                # nondegeneracy (quadratic coefficient nonzero)
    assert r.direction.startswith("forward")
    assert r.is_transcritical
    # centre-manifold normal form predicts exactly the slope of the exact nonzero branch
    assert r.dP_dk_normal_form == pytest.approx(r.dP_dk_exact, rel=1e-10)


def test_eigenvalue_crosses_zero_with_predicted_speed():
    p = Params(); kc = critical_k(p); r = transcritical_report(p)

    def lam(k):
        pk = p.with_(k=k)
        return np.max(np.linalg.eigvals(jacobian(tumour_free(pk).state, pk)).real)

    assert lam(kc * 0.999) > 0 > lam(kc * 1.001)
    assert lam(kc) == pytest.approx(0.0, abs=1e-14)
    h = 1e-6
    assert (lam(kc + h) - lam(kc - h)) / (2 * h) == pytest.approx(r.alpha, rel=1e-6)


def test_branches_intersect_and_stability_is_exchanged():
    p = Params(); kc = critical_k(p)
    # intersection: the nonzero branch passes through the tumour-free state at k_c
    assert float(branch_k(0.0, p)) == pytest.approx(kc, rel=1e-14)
    for sgn in (-1, +1):
        pk = p.with_(k=kc * (1 + sgn * 1e-3))
        x, eigs, phys = extended_root_near_zero(pk)
        tf = np.linalg.eigvals(jacobian(tumour_free(pk).state, pk))
        assert (np.max(eigs.real) < 0) != (np.max(tf.real) < 0)       # exactly one of the two is stable
        assert phys == (sgn < 0)                                       # nonzero state admissible only below k_c
        assert (np.max(tf.real) < 0) == (sgn > 0)
        # eigenvalues are opposite to leading order (normal form: lambda_nz = -lambda_TF)
        assert np.max(eigs.real) == pytest.approx(-np.max(tf.real), rel=0.05)


def test_degenerate_eta_rho_zero_the_slope_and_flip_direction():
    p = Params()
    es = degenerate_eta(p)
    assert abs(kprime0(p.with_(eta=es))) < 1e-9
    assert kprime0(p.with_(eta=0.8 * es)) < 0 < kprime0(p.with_(eta=1.2 * es))
    rs = degenerate_rho(p)
    assert abs(kprime0(p.with_(rho=rs))) < 1e-9
    assert kprime0(p.with_(rho=1.5 * rs)) < 0 < kprime0(p.with_(rho=0.5 * rs))
    assert scan_branch(p.with_(eta=1.2 * es)).direction.startswith("backward")


def test_continuation_tracks_the_exact_branch_and_locates_kc():
    p = Params(); kc = critical_k(p)
    brs = trace_all_branches(p, (0.0, 3 * kc), n_seed=9)
    nz = [b for b in brs if np.any(np.abs(b.x[:, 0]) > 1e-12)]
    assert nz
    for b in nz:
        pos = np.abs(b.x[:, 0]) > 1e-12
        assert np.max(np.abs(b.lam[pos] - branch_k(b.x[pos, 0], p))) < 1e-12
        assert b.residual.max() < 1e-12
    ev = [e for b in brs for e in b.events if e.kind.startswith("branch point")]
    assert ev and all(abs(e.lam - kc) < 1e-6 for e in ev)
    # no fold and no Hopf in the baseline physical region
    assert not [e for b in brs for e in b.events if e.kind in ("fold", "Hopf candidate") and e.x[0] > 0]


def test_continuation_finds_fold_and_bistability_for_large_eta():
    p = Params(eta=0.6); kc = critical_k(p)
    brs = trace_all_branches(p, (0.0, 1.2), n_seed=21)
    folds = [e for b in brs for e in b.events if e.kind == "fold" and e.x[0] > 0]
    assert len(folds) == 1
    Pf = fold_P(p)[0]
    assert folds[0].lam == pytest.approx(float(branch_k(Pf, p)), abs=1e-8)
    assert folds[0].x[0] == pytest.approx(Pf, rel=1e-6)
    assert scan_branch(p).bistable
    # at k inside (k_c, k_fold) there are two stable states (tumour-free + one nonzero) and one saddle
    pk = p.with_(k=0.5 * (kc + folds[0].lam))
    eqs = enumerate_equilibria(pk)
    n_unst = [int(np.sum(np.linalg.eigvals(jacobian(e.state, pk)).real > 0)) for e in eqs]
    assert sorted(n_unst) == [0, 1]


def test_fold_stability_changes_across_fold():
    p = Params(eta=0.6)
    brs = trace_all_branches(p, (0.0, 1.2), n_seed=21)
    b = [b for b in brs if any(e.kind == "fold" for e in b.events)][0]
    i = np.argmax(b.lam)                                  # the fold is the rightmost point of this branch
    assert b.max_real[max(i - 5, 0)] * b.max_real[min(i + 5, len(b.lam) - 1)] < 0     # sign change through the fold


def test_no_hopf_in_baseline_scan_and_scan_matches_baseline_structure():
    sc = scan_branch(Params())
    assert sc.direction.startswith("forward") and not sc.folds and not sc.hopf_candidates and not sc.bistable


def test_continuation_numerical_convergence_in_step_size():
    p = Params(); kc = critical_k(p)
    seed = enumerate_equilibria(p.with_(k=0.5 * kc))[0]
    bps = []
    for ds in (0.1, 0.02):
        b = continue_branch(p.with_(k=0.5 * kc), seed.state, "k", +1, (0.0, 1.0), ds0=min(ds, 0.02), ds_max=ds, lam0=0.5 * kc)
        bps.append([e.lam for e in b.events if e.kind.startswith("branch point")][0])
        pos = np.abs(b.x[:, 0]) > 1e-12
        assert np.max(np.abs(b.lam[pos] - branch_k(b.x[pos, 0], p))) < 1e-12
    assert all(abs(v - kc) < 1e-7 for v in bps)
    assert abs(bps[0] - bps[1]) < 1e-7


def test_continuation_in_other_parameter_agrees_with_exact_branch():
    # continue in eta (not k) and compare with the exact branch relation k(P; eta) = k_fixed
    p = Params(k=0.1)
    seed = enumerate_equilibria(p)[0]
    b = continue_branch(p, seed.state, "eta", +1, (0.0, 0.5), ds0=0.01, ds_max=0.05, lam0=p.eta)
    for lam, x in zip(b.lam[::10], b.x[::10]):
        assert branch_k(x[0], p.with_(eta=lam)) == pytest.approx(0.1, abs=1e-10)


def test_nonzero_nonzero_bistability_below_kc_and_tumour_free_coexistence_above():
    from dormancy.bifurcation import coexisting_stable_nonzero, stable_overlap_windows
    # eta = 0.164: folds at 0.2561 (< k_c) and 0.2840 (> k_c)
    p = Params(eta=0.164); kc = critical_k(p)
    sc = scan_branch(p, n=1201)
    assert len(sc.folds) == 2 and sc.folds[0][0] < kc < sc.folds[1][0]
    assert coexisting_stable_nonzero(p, sc)
    lo, hi = stable_overlap_windows(sc)[0]
    assert hi == pytest.approx(kc, abs=1e-4)                       # window ends at k_c (stable branch loses stability there)
    pk = p.with_(k=0.27)                                             # inside the nonzero-nonzero window
    labels = sorted(int(np.sum(np.linalg.eigvals(jacobian(e.state, pk)).real > 0)) for e in enumerate_equilibria(pk))
    assert labels == [0, 0, 1]
    assert np.max(np.linalg.eigvals(jacobian(tumour_free(pk).state, pk)).real) > 0      # tumour-free is NOT stable here
    pk2 = p.with_(k=0.2837)                                          # above k_c: tumour-free stable + one stable nonzero
    assert np.max(np.linalg.eigvals(jacobian(tumour_free(pk2).state, pk2)).real) < 0
    assert sorted(int(np.sum(np.linalg.eigvals(jacobian(e.state, pk2)).real > 0)) for e in enumerate_equilibria(pk2)) == [0, 1]
    # eta = 0.12: both folds below k_c, so no tumour-free coexistence at all
    p12 = Params(eta=0.12)
    sc12 = scan_branch(p12, n=1201)
    assert all(f[0] < critical_k(p12) for f in sc12.folds) and not sc12.bistable and coexisting_stable_nonzero(p12, sc12)
    # baseline: neither
    assert not coexisting_stable_nonzero(Params())
