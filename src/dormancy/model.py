"""Right-hand side and analytic Jacobian of the P-D-E system.

dP/dt = r P (1 - (P+D)/K) - q P + a D - k E P
dD/dt = q P - a D - delta_D D
dE/dt = s + rho P/(h+P) E - mu E - eta (P+D) E
"""
from __future__ import annotations

import numpy as np

from .parameters import Params


def rhs(t, y, p: Params):
    P, D, E = y
    N = P + D
    dP = p.r * P * (1 - N / p.K) - p.q * P + p.a * D - p.k * E * P
    dD = p.q * P - (p.a + p.delta_D) * D
    dE = p.s + p.rho * P / (p.h + P) * E - p.mu * E - p.eta * N * E
    return [dP, dD, dE]


def jacobian(y, p: Params) -> np.ndarray:
    P, D, E = y
    N = P + D
    return np.array([
        [p.r * (1 - N / p.K) - p.r * P / p.K - p.q - p.k * E, -p.r * P / p.K + p.a, -p.k * P],
        [p.q, -(p.a + p.delta_D), 0.0],
        [p.rho * p.h / (p.h + P) ** 2 * E - p.eta * E, -p.eta * E,
         p.rho * P / (p.h + P) - p.mu - p.eta * N],
    ])


def jacobian_k(y, p: Params) -> np.ndarray:
    """Derivative of the Jacobian with respect to the killing parameter k (analytic)."""
    P, D, E = y
    out = np.zeros((3, 3))
    out[0, 0] = -E
    out[0, 2] = -P
    return out


def drhs_dk(y, p: Params) -> np.ndarray:
    """Partial derivative of the vector field with respect to k: (-E P, 0, 0)."""
    P, D, E = y
    return np.array([-E * P, 0.0, 0.0])


def hessian_form(y, w, p: Params) -> np.ndarray:
    """Second derivative D^2 F(y)[w, w] (analytic), a 3-vector.

    Only the logistic term, the mass-action kill term and the immune equation are nonlinear.
    """
    P, D, E = y
    wP, wD, wE = w
    g1 = p.h / (p.h + P) ** 2            # d/dP of P/(h+P)
    g2 = -2.0 * p.h / (p.h + P) ** 3     # d2/dP2
    f1 = -p.r / p.K * (2.0 * wP * wP + 2.0 * wP * wD) - 2.0 * p.k * wE * wP
    f3 = p.rho * (g2 * E * wP * wP + 2.0 * g1 * wP * wE) - 2.0 * p.eta * (wP + wD) * wE
    return np.array([f1, 0.0, f3])

