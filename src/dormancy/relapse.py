"""Analytical approximations of the first-passage (relapse) time of N = P + D to N_detect.

Setting: residual disease near the tumour-free state x_TF = (0, 0, s/mu), k below k_c so that
lambda_+ > 0 is the dominant eigenvalue of the (P, D) block of the Jacobian at x_TF.

Estimators (all return months; np.inf means "predicts no relapse"):
  A0  original            T = ln(N_det / N0) / lambda_+
  A1  projected linear    T = ln(N_det / (n_+ c_+)) / lambda_+      c_+ = l_+ . (y0 - x_TF), n_+ = 1.w_+
  A1b two-mode linear     smallest root of N_lin(t) = N_det, N_lin from the exact linear solution
  A2  unstable manifold   second-order (Bernoulli) reduction  u' = lambda u + gamma u^2  on the
                          1-D unstable manifold of x_TF, with the immune response slaved to it:
                          closed form  T = (1/lambda) ln[(1/u0 + gamma/lambda)/(1/u_det + gamma/lambda)]

Derivations and validity conditions: report/phase_2_mathematical_analysis.md, section 6.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .equilibria import branch_k, tumour_free
from .model import hessian_form, jacobian
from .parameters import Params


@dataclass
class UnstableMode:
    lam: float            # dominant real eigenvalue (> 0 required for relapse)
    w: np.ndarray         # right eigenvector, scaled so w_P = 1
    l: np.ndarray         # left eigenvector, l.w = 1
    n1: float             # N-component of w: w_P + w_D
    gamma: float          # quadratic coefficient (1/2) l.D2F(w,w)
    n2: float             # quadratic coefficient of N(u) on the unstable manifold
    valid: bool
    h2: np.ndarray | None = None   # second-order manifold coefficient: x = x_TF + u w + u^2 h2


def unstable_mode(p: Params) -> UnstableMode:
    x0 = tumour_free(p).state
    J = jacobian(x0, p)
    ev, V = np.linalg.eig(J)
    i = int(np.argmax(np.real(ev)))
    lam = float(np.real(ev[i]))
    if lam <= 0 or abs(np.imag(ev[i])) > 1e-12:
        return UnstableMode(lam, np.zeros(3), np.zeros(3), 0.0, 0.0, 0.0, False)
    w = np.real(V[:, i]); w = w / w[0]
    el, U = np.linalg.eig(J.T)
    j = int(np.argmin(np.abs(el - ev[i])))
    l = np.real(U[:, j]); l = l / (l @ w)
    H = hessian_form(x0, w, p)
    gamma = 0.5 * float(l @ H)
    h2 = np.zeros(3)
    try:
        h2 = np.linalg.solve(2 * lam * np.eye(3) - J, 0.5 * H - gamma * w)
        n2 = float(h2[0] + h2[1])
    except np.linalg.LinAlgError:
        n2 = 0.0
    return UnstableMode(lam, w, l, float(w[0] + w[1]), gamma, n2, True, h2)


def _dev(p: Params, y0):
    return np.asarray(y0, float) - tumour_free(p).state


def T_original(p: Params, y0, N_det: float) -> float:
    """A0: ln(N_det/N0)/lambda_+ with N0 = P0 + D0."""
    m = unstable_mode(p)
    N0 = y0[0] + y0[1]
    if not m.valid or N0 <= 0:
        return np.inf
    return max(0.0, np.log(N_det / N0) / m.lam)


def T_projected(p: Params, y0, N_det: float) -> float:
    """A1: single-mode linear estimate with the exact amplitude of the unstable mode."""
    m = unstable_mode(p)
    if not m.valid:
        return np.inf
    amp = m.n1 * float(m.l @ _dev(p, y0))
    if amp <= 0:
        return np.inf
    return max(0.0, np.log(N_det / amp) / m.lam)


def T_two_mode(p: Params, y0, N_det: float, t_max: float = 1e5) -> float:
    """A1b: first crossing of N_det by the exact solution of the linearised system."""
    J = jacobian(tumour_free(p).state, p)
    ev, V = np.linalg.eig(J)
    c = np.linalg.solve(V, _dev(p, y0).astype(complex))

    def N_lin(t):
        x = (V * np.exp(ev * t)) @ c
        return float(np.real(x[0] + x[1]))

    if N_lin(0.0) >= N_det:
        return 0.0
    t_max = min(t_max, 600.0 / max(float(np.max(np.real(ev))), 1e-12))     # avoid exp overflow
    ts = np.concatenate([[0.0], np.geomspace(1e-2, t_max, 4000)])
    vals = np.array([N_lin(t) - N_det for t in ts])
    idx = np.where(vals[:-1] * vals[1:] < 0)[0]
    if len(idx) == 0:
        return np.inf
    i = idx[0]
    return float(brentq(lambda t: N_lin(t) - N_det, ts[i], ts[i + 1], xtol=1e-10))


def T_manifold(p: Params, y0, N_det: float, use_n2: bool = True) -> float:
    """A2: second-order unstable-manifold (Bernoulli) estimate. inf if the manifold saturates below N_det."""
    m = unstable_mode(p)
    if not m.valid:
        return np.inf
    u0 = float(m.l @ _dev(p, y0))
    if u0 <= 0:
        return np.inf
    n2 = m.n2 if use_n2 else 0.0
    # invert N(u) = n1 u + n2 u^2 = N_det for the smallest positive u
    if n2 == 0.0:
        u_det = N_det / m.n1
    else:
        disc = m.n1 ** 2 + 4 * n2 * N_det
        if disc < 0:
            return np.inf
        roots = [(-m.n1 + s * np.sqrt(disc)) / (2 * n2) for s in (+1, -1)]
        roots = [r for r in roots if r > 0]
        if not roots:
            return np.inf
        u_det = min(roots)
    g = m.gamma / m.lam
    num, den = 1.0 / u0 + g, 1.0 / u_det + g
    if den <= 0 or num <= 0:
        return np.inf            # u' = lambda u + gamma u^2 saturates at u* = -lambda/gamma < u_det
    return max(0.0, float(np.log(num / den) / m.lam))


def T_matched_parts(p: Params, y0, N_det: float, N_match: float = 1e-3, **sim_kw):
    """A3 decomposition: (T_total, T_closed_form, tau_numerical_tail).

    T = T_manifold(N_match) [closed form, depends on the initial condition] + tau [numerical integration
    from the second-order-manifold point at N_match to N_det; depends only on the parameters].
    The numerical tail can be most of T when N_det is large (about 2/3 at baseline, k = 0.5 k_c, N_det = 0.1),
    so agreement of A3 with a numerical first-crossing time is NOT independent evidence of a closed form.
    Returns (inf, t1, nan) if the tail never reaches N_det within the horizon.
    """
    from .solver import simulate

    m = unstable_mode(p)
    if not m.valid or N_match >= N_det:
        t = T_manifold(p, y0, N_det)
        return t, t, 0.0
    t1 = T_manifold(p, y0, N_match)
    if not np.isfinite(t1):
        return np.inf, np.inf, np.nan
    if m.n2 == 0:
        u_s = N_match / m.n1
    else:
        u_s = min(r for r in [(-m.n1 + sg * np.sqrt(m.n1 ** 2 + 4 * m.n2 * N_match)) / (2 * m.n2) for sg in (1, -1)] if r > 0)
    xs = tumour_free(p).state + u_s * m.w + u_s ** 2 * m.h2
    sim_kw.setdefault("horizon", 5e4)
    sim = simulate(p, xs, N_det, **sim_kw)
    if not sim.event:
        return np.inf, t1, np.nan
    return t1 + sim.time, t1, sim.time


def T_matched(p: Params, y0, N_det: float, N_match: float = 1e-3, **sim_kw) -> float:
    """A3 (SEMI-analytic): closed-form manifold passage to N_match plus a numerically integrated tail.

    Requires N_match << min(n1 lambda/|gamma|, phi h). See T_matched_parts for the closed-form/numerical split.
    """
    return T_matched_parts(p, y0, N_det, N_match, **sim_kw)[0]


def k_detection_limit(p: Params, N_det: float) -> float:
    """k_det = k(P_det): the equilibrium burden of the nonzero branch equals N_det (P_det = N_det/phi).

    On a forward branch (k'(0) < 0), N*(k) decreases with k, so for k in (k_det, k_c) the stable
    nonzero equilibrium lies BELOW the detection threshold: no event, and every approximation of
    the crossing time must diverge as k -> k_det from below.
    """
    return float(branch_k(N_det / p.phi, p))
