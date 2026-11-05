"""Jacobian eigenvalues, classification and the analytic tumour-free stability boundary."""
from __future__ import annotations

import numpy as np

from .equilibria import tumour_free
from .model import jacobian
from .parameters import Params


def eigenvalues(y, p: Params) -> np.ndarray:
    return np.linalg.eigvals(jacobian(y, p))


def classify(eigs, tol: float = 1e-10) -> str:
    re = np.real(eigs)
    if np.all(re < -tol):
        return "stable"
    if np.any(re > tol):
        return "unstable"
    return "marginal"


def tumour_free_threshold(p: Params) -> float:
    """Right-hand side of r < k s/mu + q delta_D/(a+delta_D)."""
    denom = p.a + p.delta_D
    dormancy_term = p.q * p.delta_D / denom if denom > 0 else 0.0
    return p.k * p.E_free + dormancy_term


def tumour_free_stable(p: Params) -> bool:
    """Analytic local stability of the tumour-free state (requires a + delta_D > 0)."""
    return p.r < tumour_free_threshold(p)


def critical_k(p: Params) -> float:
    """Immune killing k at which the tumour-free state loses stability."""
    denom = p.a + p.delta_D
    dormancy_term = p.q * p.delta_D / denom if denom > 0 else 0.0
    return (p.r - dormancy_term) / p.E_free


def tumour_free_eigenvalues(p: Params) -> np.ndarray:
    return eigenvalues(tumour_free(p).state, p)


# ---------------------------------------------------------------------------------------------
# Phase 2: full-eigenvalue stability tools
# ---------------------------------------------------------------------------------------------

def char_poly(J: np.ndarray) -> np.ndarray:
    """(a1, a2, a3) of  lambda^3 + a1 lambda^2 + a2 lambda + a3 = det(lambda I - J)."""
    a1 = -np.trace(J)
    a2 = 0.5 * (np.trace(J) ** 2 - np.trace(J @ J))
    a3 = -np.linalg.det(J)
    return np.array([a1, a2, a3])


def hopf_function(J: np.ndarray) -> float:
    """a1*a2 - a3. With a2 > 0 a sign change of this quantity is a purely imaginary pair (Hopf)."""
    a1, a2, a3 = char_poly(J)
    return float(a1 * a2 - a3)


def routh_hurwitz_stable(J: np.ndarray) -> bool:
    """Necessary and sufficient condition for all eigenvalues of a 3x3 matrix to have Re < 0."""
    a1, a2, a3 = char_poly(J)
    return bool(a1 > 0 and a3 > 0 and a1 * a2 > a3)


def max_real_eigenvalue(J_or_eigs) -> float:
    arr = np.asarray(J_or_eigs)
    eigs = np.linalg.eigvals(arr) if arr.ndim == 2 else arr
    return float(np.max(np.real(eigs)))


def n_unstable(eigs, tol: float = 1e-10) -> int:
    return int(np.sum(np.real(eigs) > tol))


def stability_label(eigs, tol: float = 1e-10) -> str:
    """'stable node/focus', 'saddle (n_u=1/2)', 'unstable node/focus' or 'marginal', from ALL eigenvalues."""
    eigs = np.asarray(eigs)
    nu = n_unstable(eigs, tol)
    nz = int(np.sum(np.abs(np.real(eigs)) <= tol))
    if nz:
        return "marginal"
    cplx = bool(np.any(np.abs(np.imag(eigs)) > 1e-12))
    kind = "focus" if cplx else "node"
    if nu == 0:
        return f"stable {kind}"
    if nu == len(eigs):
        return f"unstable {kind}"
    return f"saddle (n_u={nu})"
