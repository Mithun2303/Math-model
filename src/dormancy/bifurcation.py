"""Local bifurcation analysis of the P-D-E model.

Everything here is *derived from the implemented equations* (see report/phase_2_mathematical_analysis.md).

Key structure
-------------
On dD = 0 every nonzero equilibrium is the graph  k = k(P) = A(P) B(P) / s  (equilibria.branch_k),
with D = cP and E = s/B(P).  Hence

  * transcritical point  : P = 0 on the graph, at k_c = r_eff mu / s; the nonzero branch has slope
                           dk/dP|_0 = k'(0) = [-r phi mu / K + r_eff (eta phi - rho/h)] / s;
                           k'(0) < 0  -> nonzero branch lives at k < k_c (forward, supercritical);
                           k'(0) > 0  -> nonzero branch lives at k > k_c (backward, subcritical);
  * saddle-node (fold)   : dk/dP = 0 at an admissible P > 0 (exact: roots of a cubic);
  * Hopf                 : a1*a2 - a3 = 0 with a2 > 0 along the graph (Routh-Hurwitz).

The Sotomayor conditions for the transcritical bifurcation are evaluated with the analytic
Jacobian and the analytic second derivative D^2F, and compared with the exact k(P) branch.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .equilibria import (branch_k, branch_k_prime, branch_state, enumerate_equilibria, fold_P,
                         tumour_free)
from .model import hessian_form, jacobian, jacobian_k
from .parameters import Params
from .stability import (char_poly, critical_k, routh_hurwitz_stable, stability_label)


# ------------------------------------------------------------------------ transcritical analysis

@dataclass
class TranscriticalReport:
    k_c: float
    eigenvalues_at_kc: np.ndarray
    right_null: np.ndarray        # w, normalised w_P = 1
    left_null: np.ndarray         # l, normalised l.w = 1
    simple_zero_eigenvalue: bool  # exactly one eigenvalue at 0, the other two with Re < 0
    F_k_projection: float         # l . F_k(x_TF, k_c)   -- must be 0 (branch persists)
    alpha: float                  # l . (D_x F_k) w      -- transversality
    beta: float                   # (1/2) l . D^2F(w,w)  -- nondegeneracy
    d_lambda_dk: float            # d(lambda_0)/dk = alpha (eigenvalue speed at k_c)
    dk_dP_exact: float            # k'(0) from the explicit branch
    dP_dk_normal_form: float      # -alpha * w_P / beta  (centre-manifold prediction)
    dP_dk_exact: float            # 1 / k'(0)
    direction: str                # 'forward (supercritical)' | 'backward (subcritical)' | 'degenerate'
    is_transcritical: bool
    notes: str


def critical_eigenpair(p: Params):
    """Right and left null vectors of the Jacobian at the tumour-free state at k = k_c."""
    pc = p.with_(k=critical_k(p))
    J = jacobian(tumour_free(pc).state, pc)
    w_all, V = np.linalg.eig(J)
    i = int(np.argmin(np.abs(w_all)))
    w = np.real(V[:, i])
    wl, U = np.linalg.eig(J.T)
    j = int(np.argmin(np.abs(wl)))
    l = np.real(U[:, j])
    w = w / w[0]
    l = l / (l @ w)
    return pc, J, w_all, w, l


def transcritical_report(p: Params, tol: float = 1e-9) -> TranscriticalReport:
    """Verify the Sotomayor conditions (centre-manifold / normal-form coefficients) at k_c."""
    pc, J, eigs, w, l = critical_eigenpair(p)
    x0 = tumour_free(pc).state
    Fk = np.array([-x0[2] * x0[0], 0.0, 0.0])
    Fk_proj = float(l @ Fk)
    alpha = float(l @ (jacobian_k(x0, pc) @ w))
    beta = 0.5 * float(l @ hessian_form(x0, w, pc))
    others = np.sort(np.real(eigs))
    simple = bool(np.sum(np.abs(eigs) < 1e-9) == 1 and np.sum(np.real(eigs) < -tol) == 2)
    kp0 = float(branch_k_prime(0.0, pc))
    dPdk_nf = -alpha * w[0] / beta if beta != 0 else float("inf")
    dPdk_ex = 1.0 / kp0 if kp0 != 0 else float("inf")
    if abs(kp0) < 1e-12:
        direction = "degenerate (k'(0) = 0)"
    else:
        direction = "forward (supercritical)" if kp0 < 0 else "backward (subcritical)"
    is_tc = bool(simple and abs(Fk_proj) < 1e-12 and abs(alpha) > 1e-12 and abs(beta) > 1e-12)
    notes = ("Sotomayor: F_k = 0 at x_TF; transversality alpha != 0; nondegeneracy beta != 0; "
             "simple zero eigenvalue; both TF and nonzero equilibria exist on the two sides of k_c.")
    return TranscriticalReport(critical_k(p), eigs, w, l, simple, Fk_proj, alpha, beta, alpha, kp0,
                               dPdk_nf, dPdk_ex, direction, is_tc, notes)


def kprime0(p: Params) -> float:
    return float(branch_k_prime(0.0, p))


def degenerate_eta(p: Params) -> float:
    """eta at which k'(0) = 0 (all other parameters fixed): beyond it the bifurcation is backward."""
    if p.r_eff <= 0:
        return float("nan")
    return (p.rho / p.h + p.r * p.phi * p.mu / (p.K * p.r_eff)) / p.phi


def degenerate_rho(p: Params) -> float:
    """rho at which k'(0) = 0: for rho below it the bifurcation is backward."""
    if p.r_eff <= 0:
        return float("nan")
    return p.h * (p.eta * p.phi - p.r * p.phi * p.mu / (p.K * p.r_eff))


# ------------------------------------------------------------------------ branch-wise scans

def _branch_jacobian_arrays(P: np.ndarray, p: Params):
    """Vectorised Jacobian entries along the explicit branch (P array, k = k(P))."""
    c = p.q / (p.a + p.delta_D)
    D = c * P
    N = P + D
    B = p.mu + p.eta * p.phi * P - p.rho * P / (p.h + P)
    E = p.s / B
    k = branch_k(P, p)
    J = np.zeros(P.shape + (3, 3))
    J[..., 0, 0] = p.r * (1 - N / p.K) - p.r * P / p.K - p.q - k * E
    J[..., 0, 1] = -p.r * P / p.K + p.a
    J[..., 0, 2] = -k * P
    J[..., 1, 0] = p.q
    J[..., 1, 1] = -(p.a + p.delta_D)
    J[..., 2, 0] = p.rho * p.h / (p.h + P) ** 2 * E - p.eta * E
    J[..., 2, 1] = -p.eta * E
    J[..., 2, 2] = p.rho * P / (p.h + P) - p.mu - p.eta * N
    return J, k, E


def _charpoly_arrays(J):
    tr = np.trace(J, axis1=-2, axis2=-1)
    tr2 = np.trace(J @ J, axis1=-2, axis2=-1)
    a1 = -tr
    a2 = 0.5 * (tr ** 2 - tr2)
    a3 = -np.linalg.det(J)
    return a1, a2, a3


@dataclass
class BranchScan:
    """Result of an exact scan of the nonzero branch graph of one parameter set."""
    P: np.ndarray
    k: np.ndarray
    E: np.ndarray
    n_unstable: np.ndarray
    hopf_candidates: list     # (k, P, omega) where a1 a2 - a3 changes sign with a2 > 0
    neutral_saddles: list
    folds: list               # (k, P, d2k/dP2 sign)
    stable_k_ranges: list     # list of (k_lo, k_hi) of contiguous stable segments of the nonzero branch
    bistable: bool
    bistable_k_range: tuple | None
    k_c: float
    k_prime0: float
    direction: str


def scan_branch(p: Params, n: int = 801, P_max: float | None = None) -> BranchScan:
    """Exact (graph) scan of the nonzero branch for P in (0, P_max], independent of k.

    Admissible P: A(P) > 0 and B(P) > 0 (so k > 0 and E > 0).  The branch is a graph over P,
    so no root of the equilibrium equation can be missed for ANY k.
    """
    kc = critical_k(p)
    if P_max is None:
        P_max = p.r_eff * p.K / (p.r * p.phi) if p.r_eff > 0 else 0.0
    if P_max <= 0:
        return BranchScan(np.array([]), np.array([]), np.array([]), np.array([], int), [], [], [], [],
                          False, None, kc, kprime0(p), "no nonzero branch (r_eff <= 0)")
    P = np.geomspace(1e-9 * P_max, P_max * (1 - 1e-9), n)
    B = p.mu + p.eta * p.phi * P - p.rho * P / (p.h + P)
    ok = (B > 0) & (branch_k(P, p) > 0)
    P = P[ok]
    J, k, E = _branch_jacobian_arrays(P, p)
    a1, a2, a3 = _charpoly_arrays(J)
    eig = np.linalg.eigvals(J)
    nun = np.sum(np.real(eig) > 0, axis=1)
    H = a1 * a2 - a3
    hopf, neutral = [], []
    for i in range(len(P) - 1):
        if H[i] * H[i + 1] < 0 and (P[i + 1] / P[i] < 1.2):       # contiguous grid points only
            w = H[i] / (H[i] - H[i + 1])
            Pc = P[i] + w * (P[i + 1] - P[i])
            a2c = a2[i] + w * (a2[i + 1] - a2[i])
            if a2c > 0:
                hopf.append((float(branch_k(Pc, p)), float(Pc), float(np.sqrt(a2c))))
            else:
                neutral.append((float(branch_k(Pc, p)), float(Pc)))
    folds = []
    for Pf in fold_P(p):
        d2 = (branch_k_prime(Pf * (1 + 1e-6), p) - branch_k_prime(Pf * (1 - 1e-6), p)) / (2e-6 * Pf)
        folds.append((float(branch_k(Pf, p)), float(Pf), float(np.sign(d2))))
    # contiguous stable segments of the nonzero branch (in P) and their k-ranges
    stable = nun == 0
    ranges = []
    i = 0
    while i < len(P):
        if stable[i]:
            j = i
            while j + 1 < len(P) and stable[j + 1]:
                j += 1
            ranges.append((float(k[i:j + 1].min()), float(k[i:j + 1].max())))
            i = j + 1
        else:
            i += 1
    bi, bik = False, None
    # tumour-free state is stable for k > k_c; stable nonzero equilibrium at such a k = bistability
    for lo, hi in ranges:
        if hi > kc * (1 + 1e-9):
            bi, bik = True, (max(lo, kc), hi)
    kp0 = kprime0(p)
    direction = ("forward (supercritical)" if kp0 < 0 else "backward (subcritical)") if abs(kp0) > 1e-12 else "degenerate"
    return BranchScan(P, k, E, nun, hopf, neutral, folds, ranges, bi, bik, kc, kp0, direction)


def attractor_summary(p: Params) -> dict:
    """Exact equilibria at the given parameters with stability from ALL eigenvalues."""
    eqs = enumerate_equilibria(p)
    rows = []
    tf = tumour_free(p)
    items = [("tumour-free", tf)] + [(f"nonzero {i+1}", e) for i, e in enumerate(eqs)]
    n_stable = 0
    for name, e in items:
        eigs = np.linalg.eigvals(jacobian(e.state, p))
        lab = stability_label(eigs)
        n_stable += lab.startswith("stable")
        rows.append({"name": name, "P": e.P, "D": e.D, "E": e.E, "eigs": eigs, "label": lab})
    return {"rows": rows, "n_stable": int(n_stable), "n_nonzero": len(eqs)}


def extended_root_near_zero(p: Params):
    """Real root of the equilibrium cubic closest to P = 0 (any sign), with Jacobian eigenvalues.

    For k just above k_c on a forward branch this is the P < 0 continuation of the nonzero branch.
    It is NOT an admissible equilibrium (P < 0) but is needed for the local exchange-of-stability
    argument.  Returns (x, eigs, physical) or None.
    """
    from .equilibria import cubic_coefficients
    roots = np.roots(cubic_coefficients(p))
    real = [z.real for z in roots if abs(z.imag) < 1e-9 * (1 + abs(z.real)) and z.real > -p.h * 0.99]
    if not real:
        return None
    P = min(real, key=abs)
    x = branch_state(P, p)
    eigs = np.linalg.eigvals(jacobian(x, p))
    return x, eigs, bool(P > 0)


def stable_overlap_windows(sc: BranchScan, min_width: float = 1e-9) -> list:
    """k-windows in which TWO OR MORE disjoint stable segments of the nonzero branch coexist.

    On such a window at least two stable nonzero equilibria exist at the same k (nonzero-nonzero
    multistability). Independent of the tumour-free state: below k_c the tumour-free state is a saddle,
    so for k < k_c this is bistability between two *nonzero* states (e.g. a low-burden focus and a
    high-burden node), not tumour-free versus nonzero.
    """
    out = []
    r = sc.stable_k_ranges
    for i in range(len(r)):
        for j in range(i + 1, len(r)):
            lo, hi = max(r[i][0], r[j][0]), min(r[i][1], r[j][1])
            if hi - lo > min_width:
                out.append((float(lo), float(hi)))
    return out


def coexisting_stable_nonzero(p: Params, sc: BranchScan | None = None) -> bool:
    """True if two stable nonzero equilibria coexist at some k; verified with the exact cubic at the window midpoint."""
    sc = scan_branch(p) if sc is None else sc
    for lo, hi in stable_overlap_windows(sc):
        pk = p.with_(k=0.5 * (lo + hi))
        n_stable = sum(stability_label(np.linalg.eigvals(jacobian(e.state, pk))).startswith("stable")
                       for e in enumerate_equilibria(pk))
        if n_stable >= 2:
            return True
    return False
