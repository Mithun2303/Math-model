import numpy as np
import pytest

from dormancy import Params, critical_k, jacobian, tumour_free
from dormancy.stability import (char_poly, hopf_function, max_real_eigenvalue, n_unstable,
                                routh_hurwitz_stable, stability_label)


def test_stability_label_uses_all_eigenvalues():
    assert stability_label(np.array([-1, -2, -3])) == "stable node"
    assert stability_label(np.array([-0.1 + 1j, -0.1 - 1j, -2])) == "stable focus"
    assert stability_label(np.array([0.3, -1, -2])) == "saddle (n_u=1)"
    assert stability_label(np.array([0.3, 0.1, -2])) == "saddle (n_u=2)"
    assert stability_label(np.array([0.1 + 1j, 0.1 - 1j, 0.2])) == "unstable focus"
    assert stability_label(np.array([-1, 0.0, -2])) == "marginal"
    # one positive eigenvalue hidden behind two very stable ones must not be classed stable
    assert not stability_label(np.array([-100, -100, 1e-3])).startswith("stable")


def test_routh_hurwitz_matches_eigenvalues_on_random_matrices():
    rng = np.random.default_rng(0)
    for _ in range(500):
        J = rng.normal(size=(3, 3))
        assert routh_hurwitz_stable(J) == bool(np.all(np.linalg.eigvals(J).real < 0))


def test_hopf_function_vanishes_for_imaginary_pair_and_signals_crossing():
    w, c = 0.7, 2.0
    for eps, sign in [(-0.05, -1), (0.0, 0), (0.05, +1)]:
        # block matrix with eigenvalues eps +- i w and -c
        J = np.array([[eps, -w, 0], [w, eps, 0], [0, 0, -c]])
        a1, a2, a3 = char_poly(J)
        assert a2 > 0
        H = hopf_function(J)
        if sign == 0:
            assert abs(H) < 1e-12
        else:
            # a1*a2 - a3 > 0 <=> stable; sign flips through the Hopf point
            assert np.sign(H) == -sign


def test_char_poly_matches_numpy():
    J = np.random.default_rng(1).normal(size=(3, 3))
    assert np.allclose(char_poly(J), np.poly(J)[1:])


def test_tumour_free_eigenvalues_split_block_triangular():
    p = Params()
    ev = np.linalg.eigvals(jacobian(tumour_free(p).state, p))
    assert np.any(np.isclose(ev, -p.mu))
    assert n_unstable(ev) == 1 and max_real_eigenvalue(ev) > 0     # k = 0.1 < k_c: tumour-free unstable


def test_critical_k_formula_and_beta_property():
    p = Params()
    assert critical_k(p) == pytest.approx(p.r_eff * p.mu / p.s)
    assert p.beta == pytest.approx(p.k * p.s / p.mu)
    assert p.with_(k=critical_k(p)).beta == pytest.approx(p.r_eff)
    # bracket: sign of max real eigenvalue flips exactly at k_c
    kc = critical_k(p)
    signs = [np.sign(max_real_eigenvalue(jacobian(tumour_free(p.with_(k=kc * f)).state, p.with_(k=kc * f))))
             for f in (0.999, 1.001)]
    assert signs == [1, -1]
