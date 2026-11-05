import numpy as np
import pytest

from dormancy import Params, critical_k, enumerate_equilibria, jacobian, rhs
from dormancy.equilibria import (branch_k, branch_k_prime, branch_state, cubic_coefficients, fold_P,
                                 multistart_equilibria, nonzero_equilibria, residual)

RNG = np.random.default_rng(7)


def random_params(n):
    out = []
    while len(out) < n:
        p = Params(a=10 ** RNG.uniform(-3.5, 0), q=10 ** RNG.uniform(-3, 0), rho=10 ** RNG.uniform(-3, 0.3),
                   eta=10 ** RNG.uniform(-3, 0.7), k=RNG.uniform(0, 0.8))
        out.append(p)
    return out


def test_equilibrium_residuals_are_machine_small():
    for p in random_params(200):
        for e in enumerate_equilibria(p):
            assert residual(e, p) < 1e-12
            assert min(e.P, e.D, e.E) > 0


def test_branch_inverse_function_matches_roots():
    p = Params(k=0.1)
    for e in enumerate_equilibria(p):
        assert branch_k(e.P, p) == pytest.approx(p.k, rel=1e-10)


def test_cubic_enumeration_agrees_with_independent_multistart_and_phase1_grid():
    for f in (0.05, 0.5, 0.95):
        p = Params().with_(k=f * critical_k(Params()))
        a = enumerate_equilibria(p)
        b = multistart_equilibria(p, n_starts=150, seed=3)
        c = nonzero_equilibria(p)
        assert len(a) == len(b) == len(c) == 1
        assert np.allclose(a[0].state, b[0].state, rtol=1e-6)
        assert np.allclose(a[0].state, c[0].state, rtol=1e-6)


def test_negative_roots_rejected_above_kc():
    p = Params()
    pk = p.with_(k=1.05 * critical_k(p))
    eqs, rej = enumerate_equilibria(pk, return_rejected=True)
    assert eqs == []
    assert any(reason == "P<=0" for reason, _ in rej)       # the continuation to P < 0 is rejected as nonphysical


def test_multiple_nonzero_equilibria_are_found_when_present():
    # parameters found in the Phase-2 scan with three nonzero roots (two stable, one saddle)
    p = Params(a=0.3062386609092966, q=0.0015841983962429386, rho=0.24049141728179604, eta=1.6373830589316969, k=0.15)
    eqs = enumerate_equilibria(p)
    assert len(eqs) == 3
    labels = [np.sum(np.linalg.eigvals(jacobian(e.state, p)).real > 0) for e in eqs]
    assert sorted(labels) == [0, 0, 1]


def test_fold_points_have_zero_slope():
    p = Params(eta=0.6)
    fp = fold_P(p)
    assert len(fp) == 1
    assert abs(branch_k_prime(fp[0], p)) < 1e-9
    # numerical derivative check of the analytic slope
    P = 0.03
    fd = (branch_k(P + 1e-7, p) - branch_k(P - 1e-7, p)) / 2e-7
    assert branch_k_prime(P, p) == pytest.approx(fd, rel=1e-6)


def test_branch_state_is_equilibrium_on_the_curve():
    p = Params()
    P = 0.02
    x = branch_state(P, p)
    pk = p.with_(k=float(branch_k(P, p)))
    assert np.max(np.abs(rhs(0, x, pk))) < 1e-14
