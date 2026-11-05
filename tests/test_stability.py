import numpy as np
import pytest

from dormancy import (Params, tumour_free_stable, critical_k, classify, nonzero_equilibria, jacobian)
from dormancy.stability import tumour_free_eigenvalues, tumour_free_threshold


@pytest.mark.parametrize("a,dD", [(0.005, 0.001), (0.05, 0.0), (0.01, 0.02), (0.2, 0.1)])
def test_analytic_boundary_matches_eigenvalues(a, dD):
    p = Params(a=a, delta_D=dD)
    kc = critical_k(p)
    for factor in (0.9, 1.1):
        pk = p.with_(k=kc * factor)
        eigs = tumour_free_eigenvalues(pk)
        assert (classify(eigs) == "stable") == tumour_free_stable(pk)
        assert tumour_free_stable(pk) == (factor > 1)


def test_boundary_eigenvalue_is_zero():
    p = Params()
    pk = p.with_(k=critical_k(p))
    assert np.min(np.abs(tumour_free_eigenvalues(pk))) < 1e-12


def test_no_dormancy_loss_reduces_to_k_s_over_mu():
    p = Params(delta_D=0.0, a=0.01)
    assert tumour_free_threshold(p) == pytest.approx(p.k * p.s / p.mu)


def test_immune_eigenvalue_is_minus_mu():
    p = Params()
    assert np.any(np.isclose(tumour_free_eigenvalues(p), -p.mu))


def test_nonzero_equilibria_are_valid():
    p = Params(k=0.1)
    eqs = nonzero_equilibria(p)
    assert eqs, "expected at least one persisting equilibrium for weak immune killing"
    from dormancy import rhs
    for eq in eqs:
        assert min(eq.P, eq.D, eq.E) >= 0
        assert np.max(np.abs(rhs(0, eq.state, p))) < 1e-9
