"""Tumour-free and nonzero equilibria.

At equilibrium D = qP/(a+delta_D) and E = s/(mu + eta N - rho P/(h+P)), N = P(1 + q/(a+delta_D)),
so the system reduces to a scalar equation g(P) = 0, solved by bracketing + brentq.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .model import jacobian, rhs
from .parameters import Params


@dataclass(frozen=True)
class Equilibrium:
    P: float
    D: float
    E: float

    @property
    def state(self) -> np.ndarray:
        return np.array([self.P, self.D, self.E])


def tumour_free(p: Params) -> Equilibrium:
    return Equilibrium(0.0, 0.0, p.E_free)


def _reduced(P: float, p: Params):
    """Return (D, E, g) on the dD=0 manifold, or None where E would be undefined/negative."""
    D = p.q * P / (p.a + p.delta_D)
    N = P + D
    denom = p.mu + p.eta * N - p.rho * P / (p.h + P)
    if denom <= 0:
        return None
    E = p.s / denom
    g = p.r * P * (1 - N / p.K) - p.q * P + p.a * D - p.k * E * P
    return D, E, g


def nonzero_equilibria(p: Params, n_grid: int = 20000, tol: float = 1e-9) -> list[Equilibrium]:
    """All admissible equilibria with P > 0, found by sign changes of g(P)/P on a log grid."""
    if p.a + p.delta_D <= 0:
        raise ValueError("a + delta_D must be positive for the reduction")
    # P cannot exceed K (logistic bound), so scan (1e-12, K]
    grid = np.geomspace(1e-12, p.K, n_grid)

    def gp(P):
        out = _reduced(P, p)
        return np.nan if out is None else out[2] / P

    vals = np.array([gp(P) for P in grid])
    roots = []
    for i in range(len(grid) - 1):
        v0, v1 = vals[i], vals[i + 1]
        if np.isnan(v0) or np.isnan(v1) or v0 * v1 > 0:
            continue
        roots.append(brentq(gp, grid[i], grid[i + 1], xtol=1e-14, rtol=1e-12))

    eqs = []
    for P in roots:
        D, E, _ = _reduced(P, p)
        eq = Equilibrium(P, D, E)
        resid = np.max(np.abs(rhs(0.0, eq.state, p)))
        if min(eq.P, eq.D, eq.E) >= 0 and resid < tol:
            eqs.append(eq)
    return eqs


def branch_stability(eq: Equilibrium, p: Params) -> np.ndarray:
    return np.linalg.eigvals(jacobian(eq.state, p))


# ---------------------------------------------------------------------------------------------
# Phase 2: exact branch parametrisation, exhaustive root enumeration, multistart cross-check.
#
# On dD = 0:  D = c P, c = q/(a+delta_D), N = phi P, phi = 1 + c.  dP = 0 gives (for P != 0)
#       r_eff - r phi P / K = k E(P),   E(P) = s / Bf(P),   Bf(P) = mu + eta phi P - rho P/(h+P),
# hence every nonzero equilibrium lies on the explicit curve
#       k(P) = A(P) Bf(P) / s,   A(P) = r_eff - r phi P / K        (k is linear in the equation).
# Clearing denominators gives a cubic in P, so there are at most three nonzero roots.
# ---------------------------------------------------------------------------------------------

def branch_A(P, p: Params):
    return p.r_eff - p.r * p.phi * np.asarray(P, dtype=float) / p.K


def branch_B(P, p: Params):
    P = np.asarray(P, dtype=float)
    return p.mu + p.eta * p.phi * P - p.rho * P / (p.h + P)


def branch_k(P, p: Params):
    """Killing k at which (P, cP, E(P)) is an equilibrium; the explicit branch k(P). Valid for P > -h."""
    return branch_A(P, p) * branch_B(P, p) / p.s


def branch_k_prime(P, p: Params):
    """dk/dP along the branch (analytic)."""
    P = np.asarray(P, dtype=float)
    dA = -p.r * p.phi / p.K
    dB = p.eta * p.phi - p.rho * p.h / (p.h + P) ** 2
    return (dA * branch_B(P, p) + branch_A(P, p) * dB) / p.s


def branch_state(P, p: Params) -> np.ndarray:
    """Point (P, cP, E(P)) on the branch curve (E = s/Bf; nan where Bf <= 0)."""
    P = float(P)
    B = float(branch_B(P, p))
    E = p.s / B if B > 0 else float("nan")
    return np.array([P, p.q / (p.a + p.delta_D) * P, E])


def cubic_coefficients(p: Params) -> np.ndarray:
    """Coefficients (highest power first) of  A(P) Bpoly(P) - k s (h+P),  Bpoly = Bf (h+P)."""
    A = np.poly1d([-p.r * p.phi / p.K, p.r_eff])
    Bp = np.poly1d([p.eta * p.phi, p.mu + p.eta * p.phi * p.h - p.rho, p.mu * p.h])
    return (A * Bp - np.poly1d([p.k * p.s, p.k * p.s * p.h])).coeffs


def fold_P(p: Params) -> np.ndarray:
    """Admissible P > 0 at which dk/dP = 0 (saddle-node points of the nonzero branch). Exact via a cubic."""
    A = np.poly1d([-p.r * p.phi / p.K, p.r_eff])
    Bp = np.poly1d([p.eta * p.phi, p.mu + p.eta * p.phi * p.h - p.rho, p.mu * p.h])
    num = np.polyder(A * Bp) * np.poly1d([1.0, p.h]) - A * Bp   # (A Bp / (h+P))' * (h+P)^2
    roots = [z.real for z in np.roots(num.coeffs) if abs(z.imag) < 1e-12 * (1 + abs(z.real))]
    out = []
    for P in roots:
        if P > 0 and branch_A(P, p) > 0 and branch_B(P, p) > 0:
            out.append(P)
    return np.array(sorted(out))


def _polish(P: float, p: Params, iters: int = 30) -> np.ndarray | None:
    """Newton refinement of the full 3-D system from the algebraic root."""
    x = branch_state(P, p)
    if not np.all(np.isfinite(x)):
        return None
    for _ in range(iters):
        f = np.array(rhs(0.0, x, p))
        if np.max(np.abs(f)) < 1e-15:
            break
        try:
            dx = np.linalg.solve(jacobian(x, p), -f)
        except np.linalg.LinAlgError:
            break
        x = x + dx
    return x


def residual(eq: "Equilibrium | np.ndarray", p: Params) -> float:
    """Max-norm of the vector field at an equilibrium (verification of F(x) = 0)."""
    x = eq.state if isinstance(eq, Equilibrium) else np.asarray(eq, dtype=float)
    return float(np.max(np.abs(rhs(0.0, x, p))))


def enumerate_equilibria(p: Params, tol: float = 1e-10, return_rejected: bool = False):
    """All admissible nonzero equilibria from the exact cubic, Newton-polished and residual-checked.

    Admissible: P > 0, D > 0, E > 0, residual < tol.  Negative or complex roots are rejected
    (listed in the second return value when return_rejected=True).
    """
    if p.a + p.delta_D <= 0:
        raise ValueError("a + delta_D must be positive for the reduction")
    roots = np.roots(cubic_coefficients(p))
    ok, rej = [], []
    for z in roots:
        if abs(z.imag) > 1e-9 * (1.0 + abs(z.real)):
            rej.append(("complex root", z)); continue
        P = float(z.real)
        if not P > 0:
            rej.append(("P<=0", P)); continue
        x = _polish(P, p)
        if x is None:
            rej.append(("E undefined (Bf<=0)", P)); continue
        res = residual(x, p)
        if x.min() <= 0:
            rej.append(("negative component", x)); continue
        if res > tol:
            rej.append((f"residual {res:.2e}", x)); continue
        ok.append(Equilibrium(*x))
    ok.sort(key=lambda e: e.P)
    return (ok, rej) if return_rejected else ok


def multistart_equilibria(p: Params, n_starts: int = 400, seed: int = 0, tol: float = 1e-10,
                          dedupe_rtol: float = 1e-6) -> list[Equilibrium]:
    """Independent search: Newton/hybrid root finding from random (log-uniform) initial guesses.

    Uses no knowledge of the cubic. Negative / nonphysical converged states are rejected.
    Used as a cross-check of enumerate_equilibria.
    """
    from scipy.optimize import root

    rng = np.random.default_rng(seed)
    found: list[np.ndarray] = []
    for _ in range(n_starts):
        P0 = 10 ** rng.uniform(-8, np.log10(p.K))
        D0 = P0 * 10 ** rng.uniform(-1, 1.5) * (p.q / (p.a + p.delta_D) if p.a + p.delta_D > 0 else 1.0)
        E0 = p.E_free * 10 ** rng.uniform(-1, 1)
        sol = root(lambda y: rhs(0.0, y, p), [P0, D0, E0], jac=lambda y: jacobian(y, p),
                   method="hybr", tol=1e-14)
        if not sol.success:
            continue
        x = sol.x
        if x.min() < 0 or residual(x, p) > tol:
            continue
        if x[0] < 1e-12 and x[1] < 1e-12:      # converged onto the tumour-free state
            continue
        if not any(np.allclose(x, y, rtol=dedupe_rtol, atol=1e-12) for y in found):
            found.append(x)
    found.sort(key=lambda v: v[0])
    return [Equilibrium(*x) for x in found]
