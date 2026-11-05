"""Normalised P-D-e model (Phase 2, proposal M1): removes the arbitrary absolute scale of the immune effector.

Original model (unchanged, backward compatible):
    dP/dt = r P (1 - (P+D)/K) - q P + a D - k E P
    dD/dt = q P - (a + delta_D) D
    dE/dt = s + rho P/(h+P) E - mu E - eta (P+D) E

Substitute e = E / (s/mu) = mu E / s   (so the tumour-free level is e = 1) and beta = k s / mu.
Then E = (s/mu) e and k E = beta e, and  de/dt = (mu/s) dE/dt:

    dP/dt = r P (1 - (P+D)/K) - q P + a D - beta e P
    dD/dt = q P - (a + delta_D) D
    de/dt = mu (1 - e) + rho P/(h+P) e - eta (P+D) e

The parameters k and s disappear; only beta = k s / mu remains. All of (P, D, N, burden dynamics, k_c up
to the map beta_c = r_eff, ...) are unchanged. The gauge group (E, s, k) -> (cE, c s, k/c) of the
original model acts trivially on (P, D, beta, e).

Equivalent statements used in the tests:
    beta_c = r - q delta_D/(a+delta_D)  = r_eff            (k_c = beta_c mu / s)
    E_original(t) = (s/mu) e(t)  for any s > 0.
The units of e are dimensionless ("fraction of the tumour-free immune level"); beta has units 1/month.
"""
from __future__ import annotations

from dataclasses import dataclass, fields, replace

import numpy as np
from scipy.integrate import solve_ivp

from .parameters import Params


@dataclass(frozen=True)
class NormalizedParams:
    """Parameters of the normalised model. Time in months; burdens normalised to K; e dimensionless."""
    K: float = 1.0
    r: float = 0.12
    q: float = 0.04
    a: float = 0.005
    delta_D: float = 0.001
    beta: float = 0.04     # effective immune killing k s / mu  [1/month]
    mu: float = 0.05
    rho: float = 0.04
    h: float = 0.02
    eta: float = 0.08

    def __post_init__(self):
        for f in fields(self):
            if getattr(self, f.name) < 0:
                raise ValueError(f"parameter {f.name} must be nonnegative")
        for name in ("K", "mu", "h"):
            if getattr(self, name) <= 0:
                raise ValueError(f"parameter {name} must be positive")

    def with_(self, **changes) -> "NormalizedParams":
        return replace(self, **changes)

    @property
    def beta_c(self) -> float:
        """Critical effective killing: tumour-free state stable iff beta > beta_c = r - q delta_D/(a+delta_D)."""
        d = self.a + self.delta_D
        return self.r - (self.q * self.delta_D / d if d > 0 else 0.0)

    def to_params(self, s: float | None = None) -> Params:
        """Map to the original model with immune source s (default s = mu, i.e. E = e) and k = beta mu / s."""
        s = self.mu if s is None else s
        if s <= 0:
            raise ValueError("s must be positive")
        return Params(K=self.K, r=self.r, q=self.q, a=self.a, delta_D=self.delta_D,
                      k=self.beta * self.mu / s, s=s, mu=self.mu, rho=self.rho, h=self.h, eta=self.eta)

    @staticmethod
    def from_params(p: Params) -> "NormalizedParams":
        return NormalizedParams(K=p.K, r=p.r, q=p.q, a=p.a, delta_D=p.delta_D, beta=p.k * p.s / p.mu,
                                mu=p.mu, rho=p.rho, h=p.h, eta=p.eta)


def to_normalized_state(y, p: Params) -> np.ndarray:
    """(P, D, E) -> (P, D, e) with e = mu E / s."""
    y = np.asarray(y, float)
    return np.array([y[0], y[1], y[2] * p.mu / p.s])


def to_original_state(z, p: Params) -> np.ndarray:
    """(P, D, e) -> (P, D, E) with E = (s/mu) e."""
    z = np.asarray(z, float)
    return np.array([z[0], z[1], z[2] * p.s / p.mu])


def gauge_transform(p: Params, c: float) -> Params:
    """Original-model parameters equivalent to p under the E-gauge: s -> c s, k -> k / c (E -> c E)."""
    if c <= 0:
        raise ValueError("gauge factor c must be positive")
    return p.with_(s=c * p.s, k=p.k / c)


def rhs_normalized(t, z, n: NormalizedParams):
    """Right-hand side of the normalised system (independent implementation, derived above)."""
    P, D, e = z
    N = P + D
    dP = n.r * P * (1 - N / n.K) - n.q * P + n.a * D - n.beta * e * P
    dD = n.q * P - (n.a + n.delta_D) * D
    de = n.mu * (1 - e) + n.rho * P / (n.h + P) * e - n.eta * N * e
    return [dP, dD, de]


def jacobian_normalized(z, n: NormalizedParams) -> np.ndarray:
    P, D, e = z
    N = P + D
    return np.array([
        [n.r * (1 - N / n.K) - n.r * P / n.K - n.q - n.beta * e, -n.r * P / n.K + n.a, -n.beta * P],
        [n.q, -(n.a + n.delta_D), 0.0],
        [n.rho * n.h / (n.h + P) ** 2 * e - n.eta * e, -n.eta * e, -n.mu + n.rho * P / (n.h + P) - n.eta * N],
    ])


def simulate_normalized(n: NormalizedParams, z0, N_detect: float = 0.1, horizon: float = 240.0,
                        method: str = "DOP853", rtol: float = 1e-9, atol: float = 1e-12, n_points: int = 2000):
    """First-upward-crossing relapse simulation of the normalised model; returns (t, z, event, time)."""
    z0 = np.asarray(z0, float)
    if z0[0] + z0[1] >= N_detect:
        return np.array([0.0]), z0[:, None], 1, 0.0

    def crossing(t, z, n):
        return z[0] + z[1] - N_detect
    crossing.terminal = True
    crossing.direction = 1
    sol = solve_ivp(rhs_normalized, (0.0, horizon), z0, args=(n,), method=method, rtol=rtol, atol=atol,
                    events=crossing, t_eval=np.linspace(0.0, horizon, n_points))
    if not sol.success:
        raise RuntimeError(sol.message)
    if sol.t_events[0].size:
        t_ev = float(sol.t_events[0][0])
        keep = sol.t < t_ev
        return (np.append(sol.t[keep], t_ev), np.column_stack([sol.y[:, keep], sol.y_events[0][0]]), 1, t_ev)
    return sol.t, sol.y, 0, horizon
