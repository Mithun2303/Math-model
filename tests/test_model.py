import numpy as np
import sympy as sp

from dormancy import Params, rhs, jacobian, tumour_free


def test_tumour_free_is_equilibrium():
    p = Params()
    assert np.allclose(rhs(0, tumour_free(p).state, p), 0.0)


def test_jacobian_matches_symbolic_and_finite_difference():
    P, D, E = sp.symbols("P D E")
    r, K, q, a, dD, k, s, mu, rho, h, eta = sp.symbols("r K q a dD k s mu rho h eta")
    F = sp.Matrix([
        r * P * (1 - (P + D) / K) - q * P + a * D - k * E * P,
        q * P - a * D - dD * D,
        s + rho * P / (h + P) * E - mu * E - eta * (P + D) * E,
    ])
    J_sym = F.jacobian([P, D, E])
    p = Params()
    subs = dict(r=p.r, K=p.K, q=p.q, a=p.a, dD=p.delta_D, k=p.k, s=p.s, mu=p.mu,
                rho=p.rho, h=p.h, eta=p.eta)
    y = np.array([0.03, 0.2, 0.3])
    J_sym_num = np.array(J_sym.subs(subs).subs({P: y[0], D: y[1], E: y[2]}), dtype=float)
    assert np.allclose(jacobian(y, p), J_sym_num, atol=1e-12)

    eps = 1e-7
    J_fd = np.column_stack([
        (np.array(rhs(0, y + eps * e, p)) - np.array(rhs(0, y - eps * e, p))) / (2 * eps)
        for e in np.eye(3)
    ])
    assert np.allclose(jacobian(y, p), J_fd, atol=1e-6)


def test_invalid_params_rejected():
    import pytest
    with pytest.raises(ValueError):
        Params(r=-0.1)
    with pytest.raises(ValueError):
        Params(mu=0.0)
