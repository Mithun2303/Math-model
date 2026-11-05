"""Pseudo-arclength continuation of equilibria of the P-D-E system.

Solves F(x, lam) = 0 for x = (P, D, E), tracing one-dimensional solution branches in the
(x, lam) space through folds. The continuation parameter ``lam`` is the immune killing k by
default; any other model parameter can be used (the partial derivative is then taken by central
differences).

Features required by the Phase-2 brief:
  * adaptive step size (Newton iteration count, tangent angle, and proximity to monitored events);
  * the full 3x3 Jacobian and all eigenvalues at every accepted point;
  * event detection with bisection refinement: a zero real eigenvalue (fold or branch-point),
    a Hopf candidate (a1*a2 - a3 changes sign with a2 > 0, i.e. a complex pair), and the branch
    crossing P = 0 (the branch leaves the admissible orthant);
  * explicit residual verification |F| at every point, and a ``physical`` flag (P, D, E >= 0).

Nonphysical points (e.g. the continuation of the nonzero branch to P < 0 for k > k_c) are
*tracked* for the local analysis but flagged and never reported as admissible equilibria.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .equilibria import enumerate_equilibria
from .model import drhs_dk, jacobian, rhs
from .parameters import Params
from .stability import char_poly, hopf_function, stability_label


# --------------------------------------------------------------------------- data structures

@dataclass
class Event:
    kind: str            # 'fold', 'branch point (P=0, zero eigenvalue)', 'zero-eigenvalue (non-fold)',
                         # 'Hopf candidate', 'neutral saddle', 'P=0 crossing'
    lam: float
    x: np.ndarray
    eigs: np.ndarray
    detail: str = ""


@dataclass
class Branch:
    param: str
    lam: np.ndarray            # parameter values (k)
    x: np.ndarray              # (n, 3) states
    eigs: np.ndarray           # (n, 3) complex eigenvalues of the Jacobian
    residual: np.ndarray       # max |F| at each point
    physical: np.ndarray       # bool: P, D, E >= 0
    labels: list[str]          # stability label from all eigenvalues
    events: list[Event] = field(default_factory=list)
    termination: str = ""
    n_steps: int = 0
    n_rejected: int = 0

    @property
    def max_real(self) -> np.ndarray:
        return np.max(np.real(self.eigs), axis=1)

    @property
    def N(self) -> np.ndarray:
        return self.x[:, 0] + self.x[:, 1]


# --------------------------------------------------------------------------- helpers

def _with_param(p: Params, param: str, value: float) -> Params:
    return p.with_(**{param: float(value)}) if value >= 0 else _unchecked(p, param, value)


def _unchecked(p: Params, param: str, value: float) -> Params:
    # allow tiny negative excursions of the continuation parameter without tripping validation
    obj = object.__new__(Params)
    for f in p.__dataclass_fields__:
        object.__setattr__(obj, f, getattr(p, f))
    object.__setattr__(obj, param, float(value))
    return obj


def _dF_dlam(x, p: Params, param: str) -> np.ndarray:
    if param == "k":
        return drhs_dk(x, p)
    v = getattr(p, param)
    h = 1e-6 * max(abs(v), 1e-3)
    fp = np.array(rhs(0.0, x, _with_param(p, param, v + h)))
    fm = np.array(rhs(0.0, x, _with_param(p, param, v - h)))
    return (fp - fm) / (2 * h)


def _scales(p: Params, lam_ref: float):
    return np.array([max(p.K, 1e-12), max(p.K, 1e-12), max(p.E_free, 1e-12)]), max(abs(lam_ref), 1e-3)


class _Problem:
    """Scaled augmented problem u = (x/sx, lam/sl)."""

    def __init__(self, p: Params, param: str, lam_ref: float):
        self.p, self.param = p, param
        self.sx, self.sl = _scales(p, lam_ref)

    def split(self, u):
        return u[:3] * self.sx, u[3] * self.sl

    def join(self, x, lam):
        return np.append(np.asarray(x, float) / self.sx, lam / self.sl)

    def pk(self, lam):
        return _with_param(self.p, self.param, lam)

    def F(self, u):
        x, lam = self.split(u)
        return np.array(rhs(0.0, x, self.pk(lam)))

    def DF(self, u):
        x, lam = self.split(u)
        pk = self.pk(lam)
        Jx = jacobian(x, pk) * self.sx[None, :]
        Jl = _dF_dlam(x, pk, self.param) * self.sl
        return np.column_stack([Jx, Jl])

    def tangent(self, u, t_prev=None, direction=1.0):
        """Unit tangent: null vector of DF (3x4). If the null space is more than one-dimensional
        (a branch point), the direction closest to the previous tangent is taken, so the
        continuation proceeds straight through the intersection on the same branch."""
        _, sv, Vt = np.linalg.svd(self.DF(u))
        null = [Vt[-1]]
        if sv[-1] < 1e-9 * sv[0] and len(sv) >= 2:
            pass
        # rows beyond the rank: Vt has 4 rows; sv has 3 entries; extra null rows are Vt[3]
        ns = [Vt[i] for i in range(3, 4)] + [Vt[i] for i in range(len(sv)) if sv[i] < 1e-10 * sv[0]]
        N = np.array(ns)
        if t_prev is None:
            t = N[0]
            if t[3] * direction < 0:
                t = -t
        else:
            t = N.T @ (N @ t_prev)
            nrm = np.linalg.norm(t)
            t = t / nrm if nrm > 0 else N[0]
            if t @ t_prev < 0:
                t = -t
        return t / np.linalg.norm(t)

    def correct(self, pred, t, max_iter=14, tol=1e-13):
        """Newton on [F(u)=0; t.(u-pred)=0]. Returns (u, n_iter) or (None, n_iter).

        Convergence is judged on the Newton step (scaled, absolute 1e-14) rather than on the
        residual alone: near the transcritical point P ~ 1e-8, so F_P ~ P*(...) is tiny for
        points that are NOT on the branch, and a residual-only test would accept them.
        """
        u = pred.copy()
        for it in range(1, max_iter + 1):
            Fv = self.F(u)
            if not np.all(np.isfinite(Fv)):
                return None, it
            G = np.append(Fv, t @ (u - pred))
            M = np.vstack([self.DF(u), t])
            try:
                du = np.linalg.solve(M, -G)
            except np.linalg.LinAlgError:
                return None, it
            u = u + du
            if not np.all(np.isfinite(u)):
                return None, it
            if np.max(np.abs(du)) < 1e-14 and np.max(np.abs(self.F(u))) < 1e-12:
                return u, it
        if np.max(np.abs(du)) < 1e-11 and np.max(np.abs(self.F(u))) < 1e-12:
            return u, max_iter
        return None, max_iter

    def refine_fixed_lam(self, u, tol=1e-14, max_iter=30):
        u = u.copy()
        for _ in range(max_iter):
            Fv = self.F(u)
            if np.max(np.abs(Fv)) < tol:
                return u
            A = self.DF(u)[:, :3]
            u[:3] += np.linalg.solve(A, -Fv)
        return u


def _monitors(prob: _Problem, u):
    x, lam = prob.split(u)
    J = jacobian(x, prob.pk(lam))
    a1, a2, a3 = char_poly(J)
    return {
        "det": float(np.linalg.det(J)),
        "hopf": float(a1 * a2 - a3),
        "a2": float(a2),
        "P": float(x[0]),
        "eigs": np.linalg.eigvals(J),
        "J": J,
    }


# --------------------------------------------------------------------------- main routine

def continue_branch(p: Params, x0, param: str = "k", direction: int = +1,
                    lam_bounds: tuple[float, float] = (0.0, 1.0),
                    ds0: float = 0.02, ds_min: float = 1e-9, ds_max: float = 0.1,
                    max_steps: int = 4000, newton_tol: float = 1e-13,
                    P_floor: float | None = None, E_ceiling_factor: float = 50.0,
                    detect_events: bool = True, lam0: float | None = None) -> Branch:
    """Trace the branch through (x0, lam0) in the given direction (+1 / -1 along the parameter)."""
    lam0 = getattr(p, param) if lam0 is None else lam0
    lo, hi = lam_bounds
    if P_floor is None:
        P_floor = -0.5 * p.h
    prob = _Problem(p, param, max(abs(hi), abs(lo)))
    u = prob.refine_fixed_lam(prob.join(x0, lam0))
    t = prob.tangent(u, None, direction)

    lam_l, xs, eig_l, res_l, phys_l, lab_l = [], [], [], [], [], []
    events: list[Event] = []

    def record(uu):
        x, lam = prob.split(uu)
        pk = prob.pk(lam)
        J = jacobian(x, pk)
        e = np.linalg.eigvals(J)
        lam_l.append(lam); xs.append(x); eig_l.append(e)
        res_l.append(float(np.max(np.abs(rhs(0.0, x, pk)))))
        phys_l.append(bool(x.min() >= -1e-14))
        lab_l.append(stability_label(e))

    record(u)
    m_prev = _monitors(prob, u)
    ds = ds0
    n_steps = n_rej = 0
    term = "max_steps"
    prev_dfds = {"det": None, "hopf": None}

    while n_steps < max_steps:
        # land exactly on the parameter bound
        lam_now = prob.split(u)[1]
        dlam = t[3] * prob.sl
        step = ds
        if abs(dlam) > 1e-14:
            bound = hi if dlam > 0 else lo
            remaining = (bound - lam_now) / (dlam / 1.0) * 1.0     # arclength to bound along t
            if remaining <= 0:
                term = f"reached {param} bound"; break
            if remaining < step:
                step = remaining
        pred = u + step * t
        u_new, nit = prob.correct(pred, t, tol=newton_tol)
        if u_new is not None and abs(dlam) > 1e-14 and step == remaining:
            u_fix = u_new.copy(); u_fix[3] = bound / prob.sl      # land exactly on the bound
            try:
                u_fix = prob.refine_fixed_lam(u_fix)
                if np.max(np.abs(prob.F(u_fix))) < newton_tol and np.linalg.norm(u_fix - u_new) < 1e-6:
                    u_new = u_fix
            except np.linalg.LinAlgError:
                pass
        if u_new is None:
            ds *= 0.5; n_rej += 1
            if ds < ds_min:
                term = "step collapsed (no convergence)"; break
            continue
        t_new = prob.tangent(u_new, t)
        cosang = float(t @ t_new)
        if cosang < 0.98 and step > 4 * ds_min:
            ds = step * 0.5; n_rej += 1           # too much curvature: shrink the step
            continue
        x_new, lam_new = prob.split(u_new)
        m_new = _monitors(prob, u_new)

        # ----- event detection and bisection refinement
        found = []
        if detect_events:
            for name in ("det", "hopf", "P"):
                fa, fb = m_prev[name], m_new[name]
                if name == "P" and (abs(fa) < 1e-13 or abs(fb) < 1e-13):
                    continue                      # tumour-free branch: P == 0 identically
                if fa * fb < 0:
                    found.append(_bisect_event(prob, u, t, step, name, m_prev[name], newton_tol))
        found = [f for f in found if f is not None]
        found.sort(key=lambda f: f[0])
        ev_here: list[Event] = []
        for s_ev, u_ev, name in found:
            x_ev, lam_ev = prob.split(u_ev)
            ev = _classify_event(name, x_ev, lam_ev, _monitors(prob, u_ev), prob.tangent(u_ev, t), p)
            if ev.kind == "P=0 crossing" and any(abs(e.lam - ev.lam) < 1e-7 * (1 + abs(ev.lam)) for e in ev_here):
                continue                          # same point as a zero-eigenvalue event: merged
            ev_here.append(ev)
            record(u_ev)
        events.extend(ev_here)

        record(u_new)
        n_steps += 1
        f_old, f_new = dict((n, m_prev[n]) for n in ("det", "hopf")), dict((n, m_new[n]) for n in ("det", "hopf"))
        u, t, m_prev = u_new, t_new, m_new

        # ----- adaptive step: iteration count, angle, and linear-extrapolation proximity to an event
        if nit <= 3 and cosang > 0.999:
            ds = min(ds * 1.6, ds_max)
        elif nit >= 6:
            ds = max(ds * 0.6, ds_min)
        for name in ("det", "hopf"):
            delta = f_new[name] - f_old[name]
            if delta != 0 and f_new[name] * (f_new[name] + delta) < 0:      # crossing predicted next step
                ds = min(ds, max(0.5 * step * abs(f_new[name] / delta), 10 * ds_min))

        # ----- terminate on leaving the region of interest
        if x_new[0] < P_floor:
            term = "P below floor (branch left the extended region)"; break
        if x_new[2] > E_ceiling_factor * max(p.E_free, 1e-9) or x_new[2] < 0:
            term = "E outside range (denominator singularity)"; break
        if (lam_new <= lo + 1e-14 and t_new[3] < 0) or (lam_new >= hi - 1e-14 and t_new[3] > 0):
            term = f"reached {param} bound"; break

    return Branch(param, np.array(lam_l), np.array(xs), np.array(eig_l), np.array(res_l),
                  np.array(phys_l), lab_l, events, term, n_steps, n_rej)


def _bisect_event(prob, u, t, step, name, m_a, tol):
    """Locate the sign change of monitor ``name`` on the corrector curve between u and u+step*t.

    Special case: when the nonzero branch runs into P = 0 (transcritical point) the corrector
    becomes singular, because every parameter value solves F_P = 0 on the line P = 0. The
    bisection then stops once |P| < 2e-7 and the point at P = 0 is obtained by a quadratic
    extrapolation (in P) of the last reliable branch points; the error is O(P^3).
    """
    lo_s, hi_s, u_ev = 0.0, step, None
    nonzero_start = abs(u[0]) > 1e-13
    trail = []
    for _ in range(100):
        mid = 0.5 * (lo_s + hi_s)
        u_try, _ = prob.correct(u + mid * t, t, tol=tol)
        if u_try is None:
            hi_s = mid; continue
        mv = _monitors(prob, u_try)[name]
        if mv * m_a < 0:
            hi_s = mid
        else:
            lo_s = mid
        u_ev = u_try
        if nonzero_start and abs(u_try[0] * prob.sx[0]) > 2e-7:
            trail.append(u_try)
        if (hi_s - lo_s) < 1e-14 or abs(mv) < 1e-18:
            break
        if nonzero_start and abs(u_try[0] * prob.sx[0]) < 2e-7:
            break
    if u_ev is None:
        return None
    if nonzero_start and abs(u_ev[0] * prob.sx[0]) < 2e-7 and len(trail) >= 3:
        pts = sorted(trail, key=lambda v: abs(v[0]))[:3]
        Ps = np.array([v[0] for v in pts])
        u_ev = np.array([np.polyval(np.polyfit(Ps, [v[i] for v in pts], 2), 0.0) for i in range(4)])
        u_ev[0] = 0.0
    return (0.5 * (lo_s + hi_s), u_ev, name)


def _classify_event(name, x, lam, m, t, p) -> Event:
    eigs = m["eigs"]
    if name == "det":
        t_lam = abs(t[3])
        if t_lam < 1e-5:
            kind, det = "fold", "tangent d(lam)/ds = 0: saddle-node"
        elif abs(x[0]) < 1e-7 * p.K:
            kind, det = "branch point (P=0, zero eigenvalue)", "branch meets the tumour-free line"
        else:
            kind, det = "zero-eigenvalue (non-fold)", "det J = 0 with nonzero d(lam)/ds"
        return Event(kind, lam, x, eigs, det)
    if name == "hopf":
        if m["a2"] > 0:
            return Event("Hopf candidate", lam, x, eigs,
                         f"a1*a2-a3 = 0, a2 = {m['a2']:.3e} > 0; max Im(eig) = {np.max(np.abs(np.imag(eigs))):.3e}")
        return Event("neutral saddle", lam, x, eigs, f"a1*a2-a3 = 0 but a2 = {m['a2']:.3e} <= 0 (real pair, no Hopf)")
    return Event("P=0 crossing", lam, x, eigs, "branch crosses P = 0")


# --------------------------------------------------------------------------- multi-seed tracing

def _on_branch(lam, x, br: Branch, rtol: float = 0.1) -> bool:
    """True if the point (x, lam) lies on the traced polyline of ``br`` (interpolated at equal lam)."""
    for i in range(len(br.lam) - 1):
        l0, l1 = br.lam[i], br.lam[i + 1]
        if min(l0, l1) - 1e-12 <= lam <= max(l0, l1) + 1e-12 and l0 != l1:
            w = (lam - l0) / (l1 - l0)
            xi = (1 - w) * br.x[i] + w * br.x[i + 1]
            if abs(xi[0] - x[0]) <= rtol * abs(x[0]) + 1e-12 and abs(xi[2] - x[2]) <= rtol * abs(x[2]) + 1e-12:
                return True
    return False


def trace_all_branches(p: Params, k_range=(0.0, 1.0), seed_ks=None, n_seed: int = 15,
                       param: str = "k", **kw) -> list[Branch]:
    """Trace the tumour-free branch and every nonzero branch reachable from the seed parameters.

    Seeds are the exact admissible equilibria (cubic roots) at several parameter values,
    including values just below/above the analytic k_c and around any fold. A seed that already
    lies on a traced branch is skipped; this discovers disconnected branches (isolas), which a
    single continuation from one starting point would miss.
    """
    lo, hi = k_range
    branches: list[Branch] = []

    def run_pair(x0, lam0, name):
        pk = _with_param(p, param, lam0)
        out = []
        for d in (+1, -1):
            b = continue_branch(pk, x0, param, d, (lo, hi), lam0=lam0, **kw)
            out.append(b)
        return out

    # tumour-free branch (exact, but traced so that its events come from the same machinery)
    mid = 0.5 * (lo + hi)
    pk = _with_param(p, param, mid)
    tf = run_pair(np.array([0.0, 0.0, pk.E_free]), mid, "tumour-free")
    for b in tf:
        branches.append(b)

    if seed_ks is None:
        seed_ks = list(np.linspace(lo, hi, n_seed))
    seeds = []
    for kk in seed_ks:
        pk = _with_param(p, param, kk)
        try:
            for e in enumerate_equilibria(pk):
                seeds.append((kk, e.state))
        except ValueError:
            pass

    nonzero: list[Branch] = []
    for kk, x in seeds:
        pk0 = _with_param(p, param, kk)
        if any(_on_branch(kk, x, b) for b in nonzero):
            continue
        nonzero.extend(run_pair(x, kk, "nonzero"))
    branches.extend(nonzero)
    return branches
