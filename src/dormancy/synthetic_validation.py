"""Phase 4: synthetic parameter recovery, profile likelihoods, practical identifiability.

All data are SYNTHETIC (known ground truth, hypothetical observation model, see observation.py). The fitter uses a
different solver and tolerance (DOP853 rtol 1e-9 with analytic forward sensitivities) than the data generator (Radau rtol 1e-12):
no inverse crime. Fitting is constrained (box bounds on log parameters) multistart least squares; the detection limit enters as a
Tobit (left-censored) likelihood; estimates ON a bound are reported as bound hits (a sign of non-identifiability, not success).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares
from scipy.special import log_ndtr

from .identifiability import PARAMS, forward_sensitivities, numerical_spectrum
from .observation import CHANNELS, Dataset, ObservationAssumptions, theta_to_model

MODEL_FREE = ["r", "q", "a", "delta_D", "beta", "rho", "eta"]
ALL_FIT = MODEL_FREE + ["N0", "fD", "c", "ce"]


@dataclass
class FitProblem:
    data: Dataset
    free: list
    fixed: dict                      # values of all non-free quantities (incl. h, mu, K and anything not in `free`)
    lo: np.ndarray = None
    hi: np.ndarray = None
    sigma: float = None              # noise sd (natural log); default: the data's noise level (assumed known)
    mu: float = 0.05
    h: float = 0.02
    K: float = 1.0
    fit_solver: dict = field(default_factory=lambda: {"method": "DOP853", "rtol": 1e-9, "atol": 1e-12})

    def __post_init__(self):
        for nm in self.free:
            if nm not in ALL_FIT:
                raise ValueError(f"cannot fit {nm}")
        if self.sigma is None:
            self.sigma = self.data.assumptions.noise_sigma_ln
        if self.lo is None:
            raise ValueError("bounds required")

    # ---- transforms: positive quantities -> log, dormant fraction -> logit
    def to_x(self, theta: dict) -> np.ndarray:
        return np.array([_fwd(nm, theta[nm]) for nm in self.free])

    def to_theta(self, x: np.ndarray) -> dict:
        th = dict(self.fixed)
        for nm, xi in zip(self.free, x):
            th[nm] = _inv(nm, xi)
        return th


def _fwd(nm, v):
    return float(np.log(v / (1 - v))) if nm == "fD" else float(np.log(v))


def _inv(nm, x):
    return float(1 / (1 + np.exp(-x))) if nm == "fD" else float(np.exp(x))


def make_bounds(truth: dict, free: list, log_factor: float = 30.0):
    """Box bounds: truth x / factor ... truth x factor on the log scale (logit scale +-ln(factor) for the dormant fraction)."""
    lf = np.log(log_factor)
    c = np.array([_fwd(nm, truth[nm]) for nm in free])
    return c - lf, c + lf


def residuals_jacobian(prob: FitProblem, x: np.ndarray, want_jac: bool = True):
    """Weighted residual vector and Jacobian d r / d x (x = log / logit parameters)."""
    th = prob.to_theta(x)
    d = prob.data
    a = d.assumptions
    p_free_model = [nm for nm in prob.free if nm in MODEL_FREE]
    ic_names = [nm for nm in prob.free if nm in ("N0", "fD")]
    P0 = th["N0"] * (1 - th["fD"]); D0 = th["N0"] * th["fD"]
    x0 = np.array([P0, D0, th.get("e0", 1.0) * prob.mu / prob.mu])        # e0 relative to tumour-free level (= 1 in the s = mu gauge)
    thm = {k: th[k] for k in ("r", "q", "a", "delta_D", "beta", "rho", "eta")}
    thm.update(h=prob.h, mu=prob.mu, K=prob.K)
    sens, ic_sens = list(p_free_model), {}
    if "N0" in ic_names:
        sens.append("N0"); ic_sens["N0"] = np.array([P0, D0, 0.0])
    if "fD" in ic_names:
        sens.append("fD"); ic_sens["fD"] = np.array([-1.0, 1.0, 0.0]) * th["N0"] * th["fD"] * (1 - th["fD"])
    try:
        X, S = forward_sensitivities(thm, x0, sens, d.times, rtol=prob.fit_solver["rtol"], atol=prob.fit_solver["atol"],
                                     ic_sens=ic_sens, method=prob.fit_solver["method"])
    except Exception:                                      # noqa: BLE001
        big = np.full(sum(len(d.times) for _ in CHANNELS[a.observed]), 1e3)
        return big, np.zeros((len(big), len(x)))
    if X.shape[0] != len(d.times) or not np.all(np.isfinite(X)) or np.any(X[:, :2] <= 0):
        big = np.full(sum(len(d.times) for _ in CHANNELS[a.observed]), 1e3)
        return big, np.zeros((len(big), len(x)))
    ncol = len(prob.free)
    # d state / d x  (columns follow prob.free); model params: chain rule d/dln(theta) = theta d/dtheta
    dX = np.zeros((len(d.times), 3, ncol))
    for j, nm in enumerate(prob.free):
        if nm in MODEL_FREE:
            dX[:, :, j] = S[:, :, sens.index(nm)] * thm[nm]
        elif nm in ("N0", "fD"):
            dX[:, :, j] = S[:, :, sens.index(nm)]
    res, jac = [], []
    for ch in CHANNELS[a.observed]:
        if ch == "N":
            val = X[:, 0] + X[:, 1]; dval = dX[:, 0, :] + dX[:, 1, :]
        elif ch == "P":
            val, dval = X[:, 0], dX[:, 0, :]
        elif ch == "D":
            val, dval = X[:, 1], dX[:, 1, :]
        else:
            val, dval = X[:, 2], dX[:, 2, :]
        scale_name, scale_val = ("ce", th.get("ce", a.immune_scale_ce)) if ch == "e" else ("c", th.get("c", a.signal_scale_c))
        m = np.log(scale_val) + np.log(val)
        dm = dval / val[:, None]
        if scale_name in prob.free:
            dm = dm.copy(); dm[:, prob.free.index(scale_name)] += 1.0
        yv, cen = d.y[ch], d.censored[ch]
        r = np.empty(len(m)); J = np.zeros((len(m), ncol))
        unc = ~cen
        r[unc] = (np.log(yv[unc]) - m[unc]) / prob.sigma
        J[unc] = -dm[unc] / prob.sigma
        if cen.any():                                      # Tobit: left-censored below the limit of detection
            z = (np.log(a.lod) - m[cen]) / prob.sigma
            nll = -log_ndtr(z)
            rc = np.sqrt(2 * np.maximum(nll, 1e-300))
            ratio = np.exp(-0.5 * z * z - 0.5 * np.log(2 * np.pi) - log_ndtr(z))      # phi(z)/Phi(z)
            dnll_dm = ratio / prob.sigma
            r[cen] = rc
            J[cen] = (dnll_dm / np.maximum(rc, 1e-12))[:, None] * dm[cen]
        res.append(r); jac.append(J)
    return np.concatenate(res), (np.vstack(jac) if want_jac else None)


@dataclass
class FitResult:
    x: np.ndarray
    cost: float
    success: bool
    nfev: int
    bound_hit: np.ndarray
    x0: np.ndarray
    J: np.ndarray = None


def fit_from(prob: FitProblem, x0: np.ndarray, max_nfev: int = 80) -> FitResult:
    f = lambda x: residuals_jacobian(prob, x, False)[0]
    j = lambda x: residuals_jacobian(prob, x, True)[1]
    sol = least_squares(f, np.clip(x0, prob.lo, prob.hi), jac=j, bounds=(prob.lo, prob.hi), method="trf", x_scale="jac",
                        max_nfev=max_nfev, xtol=1e-11, ftol=1e-11, gtol=1e-9)
    hit = (np.abs(sol.x - prob.lo) < 1e-3) | (np.abs(sol.x - prob.hi) < 1e-3)
    return FitResult(sol.x, 2 * sol.cost, bool(sol.success), int(sol.nfev), hit, np.asarray(x0), sol.jac)


def random_starts(prob: FitProblem, n: int, seed: int) -> np.ndarray:
    from scipy.stats import qmc
    U = qmc.LatinHypercube(d=len(prob.free), seed=seed).random(n)
    return prob.lo + U * (prob.hi - prob.lo)


def multistart(prob: FitProblem, n_starts: int, seed: int, max_nfev: int = 80, extra_starts=()) -> list[FitResult]:
    starts = list(random_starts(prob, n_starts, seed)) + list(extra_starts)
    out = [fit_from(prob, s, max_nfev) for s in starts]
    out.sort(key=lambda r: r.cost)
    return out


def wald_and_spectrum(prob: FitProblem, best: FitResult):
    """Gauss-Newton Fisher information at the optimum: covariance (pseudo-inverse), 95% Wald CI, singular spectrum, correlations."""
    _, J = residuals_jacobian(prob, best.x, True)
    sp = numerical_spectrum(J, 1.0)
    JTJ = J.T @ J
    cov = np.linalg.pinv(JTJ, rcond=1e-12)
    sd = np.sqrt(np.maximum(np.diag(cov), 0))
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = cov / np.outer(sd, sd)
    return {"cov": cov, "sd": sd, "corr": corr, "spectrum": sp, "J": J}


def profile_likelihood(prob: FitProblem, best: FitResult, name: str, half_width_log: float, n_points: int, seed: int = 0,
                       n_restarts: int = 2, max_nfev: int = 60):
    """Profile of ``name`` over best +- half_width_log (log/logit scale). Warm-started outward from the optimum in both directions.

    Returns dict(grid, chi2 (relative to the best fit), ci95 (lo, hi) or None per side, finite_both). A side is 'open' if the
    profile has not exceeded the 95% threshold within the grid or the box bounds (Raue et al. 2009: not identifiable on that side).
    """
    j = prob.free.index(name)
    others = [i for i in range(len(prob.free)) if i != j]
    sub_free = [prob.free[i] for i in others]
    grid = np.linspace(best.x[j] - half_width_log, best.x[j] + half_width_log, n_points)
    grid = np.clip(grid, prob.lo[j], prob.hi[j])
    mid = int(np.argmin(np.abs(grid - best.x[j])))
    chi2 = np.full(n_points, np.nan)
    rng = np.random.default_rng(seed)

    def solve(val, x_start):
        sub = FitProblem(prob.data, sub_free, {**prob.fixed, name: _inv(name, val)}, prob.lo[others], prob.hi[others], prob.sigma,
                         prob.mu, prob.h, prob.K, prob.fit_solver)
        cands = [fit_from(sub, np.clip(x_start[others], sub.lo, sub.hi), max_nfev)]
        for _ in range(n_restarts):
            xs = x_start[others] + rng.standard_normal(len(others)) * 0.3
            cands.append(fit_from(sub, np.clip(xs, sub.lo, sub.hi), max_nfev))
        b = min(cands, key=lambda r: r.cost)
        full = np.array(x_start, float); full[others] = b.x; full[j] = val
        return b.cost, full

    # outward sweeps with warm start
    for direction in (+1, -1):
        x_prev = best.x.copy()
        idxs = range(mid, n_points) if direction > 0 else range(mid, -1, -1)
        for i in idxs:
            c, x_prev = solve(grid[i], x_prev)
            chi2[i] = c
    chi2 = chi2 - best.cost
    chi2 = np.maximum(chi2, 0.0)
    thr = 3.841

    def crossing(idx_range):
        prev = None
        for i in idx_range:
            if chi2[i] > thr:
                if prev is None:
                    return grid[i]
                lo_, hi_ = grid[prev], grid[i]
                f0, f1 = chi2[prev], chi2[i]
                return lo_ + (thr - f0) / (f1 - f0) * (hi_ - lo_)
            prev = i
        return None

    hi_c, lo_c = crossing(range(mid, n_points)), crossing(range(mid, -1, -1))
    return {"name": name, "grid": grid, "chi2": chi2, "ci_lo": lo_c, "ci_hi": hi_c, "finite_both": (lo_c is not None and hi_c is not None),
            "at_bound_lo": bool(grid[0] <= prob.lo[j] + 1e-9), "at_bound_hi": bool(grid[-1] >= prob.hi[j] - 1e-9)}


def recovery_metrics(prob: FitProblem, best: FitResult, truth: dict) -> dict:
    """Parameter recovery error (log ratio), prediction error in the window, extrapolation to 240 months and crossing time."""
    est = prob.to_theta(best.x)
    out = {}
    for nm in prob.free:
        out[f"err_{nm}"] = float(np.log(est[nm] / truth[nm])) if nm != "fD" else float(_fwd("fD", est["fD"]) - _fwd("fD", truth["fD"]))
    return out


def predicted_curve(theta: dict, times, mu: float = 0.05, h: float = 0.02, K: float = 1.0, method="DOP853", rtol=1e-9, atol=1e-12):
    from scipy.integrate import solve_ivp
    from .model import rhs
    p, x0 = theta_to_model(theta, mu, h, K)
    sol = solve_ivp(rhs, (0, float(max(times))), x0, args=(p,), method=method, rtol=rtol, atol=atol, t_eval=np.asarray(times, float))
    return sol.y.T


def prediction_errors(prob: FitProblem, best: FitResult, truth: dict, horizon: float = 240.0, N_detect: float = 0.1) -> dict:
    """Held-out error of ln N inside the window (monthly grid), extrapolation error at 120/240 months, and crossing-time error."""
    from .observation import detection_time, simulate_hidden
    est = prob.to_theta(best.x)
    est_full = {**truth, **est}
    t_in = np.arange(1.0, max(prob.data.times.max(), 2.0), 1.0)
    hid = simulate_hidden(truth, t_in)
    pred = predicted_curve(est_full, t_in)
    e_in = np.sqrt(np.mean((np.log(pred[:, 0] + pred[:, 1]) - np.log(hid.N)) ** 2))
    out = {"rmse_lnN_in_window": float(e_in)}
    for T in (120.0, 240.0):
        a = simulate_hidden(truth, [0.0, T]).N[-1]
        b = predicted_curve(est_full, [0.0, T]); b = b[-1, 0] + b[-1, 1]
        out[f"abs_lnN_error_at_{int(T)}"] = float(abs(np.log(b) - np.log(a)))
    try:
        Tt = detection_time(truth, N_detect, horizon * 10)
        Te = detection_time(est_full, N_detect, horizon * 10, method="DOP853", rtol=1e-9, atol=1e-12)
    except Exception:                                      # noqa: BLE001
        Tt = Te = float("inf")
    out["T_true"], out["T_pred"] = Tt, Te
    out["T_rel_error"] = float(abs(Te - Tt) / Tt) if (np.isfinite(Tt) and np.isfinite(Te)) else float("nan")
    return out


# ----------------------------------------------------------------------------------------------- replicate runner
FREE_SETS = {"ic_unknown": MODEL_FREE + ["N0", "fD"], "ic_known": list(MODEL_FREE),
             "ic_scale_unknown": MODEL_FREE + ["N0", "fD", "c"]}


def run_recovery_replicate(args) -> dict:
    """Top-level worker: one synthetic dataset -> multistart fit -> recovery / identifiability / prediction metrics.

    args = (truth, assumptions_dict, free_set_name, seed, n_starts, bounds_log_factor, condition_name)
    """
    from .observation import ObservationAssumptions, generate_dataset
    truth, adict, fset, seed, n_starts, lf, cname = args
    a = ObservationAssumptions(**adict)
    free = FREE_SETS[fset]
    d = generate_dataset(truth, a, seed)
    fixed = {k: truth[k] for k in ("e0",) if k in truth}
    fixed.update({k: truth[k] for k in MODEL_FREE + ["N0", "fD"] if k not in free})
    if "c" not in free:
        fixed["c"] = a.signal_scale_c
    truth_fit = dict(truth); truth_fit["c"] = a.signal_scale_c
    lo, hi = make_bounds(truth_fit, free, lf)
    prob = FitProblem(d, free, fixed, lo, hi)
    xt = prob.to_x(truth_fit)
    res = multistart(prob, n_starts, seed + 12345)
    best = res[0]
    out = {"condition": cname, "seed": seed, "n_obs_times": len(d.times), "T_detect_true": d.T_detect, "obs_end": float(d.times[-1]),
           "best_cost": best.cost, "n_dof": int(sum(len(d.times) * 1 for _ in CHANNELS[a.observed]) - len(free)),
           "n_starts_within_chi2_3.84": int(sum(r.cost <= best.cost + 3.84 for r in res)), "n_starts": len(res),
           "any_bound_hit": bool(best.bound_hit.any()), "n_bound_hits": int(best.bound_hit.sum())}
    near = np.array([r.x for r in res if r.cost <= best.cost + 3.84])
    wal = wald_and_spectrum(prob, best)
    _, Jt = residuals_jacobian(prob, xt, True)
    spt = numerical_spectrum(Jt, 1.0)
    out["n_dirs_s_gt_3_at_truth"] = spt["n_above_3"]; out["n_dirs_s_gt_1_at_truth"] = spt["n_above_1"]
    out["sv_at_truth"] = ";".join(f"{v:.4g}" for v in spt["singular_values"])
    err = best.x - xt
    for j, nm in enumerate(free):
        out[f"err_{nm}"] = float(err[j])
        out[f"sd_wald_{nm}"] = float(wal["sd"][j])
        out[f"covered95_{nm}"] = bool(abs(err[j]) <= 1.96 * wal["sd"][j])
        out[f"spread_near_best_{nm}"] = float(np.std(near[:, j])) if len(near) > 1 else 0.0
        out[f"at_bound_{nm}"] = bool(best.bound_hit[j])
    # well-determined combinations (right singular vectors at truth with s > 3): estimation error in units of the predicted sd
    z = []
    for k in range(spt["n_above_3"]):
        v = spt["Vt"][k]
        z.append(float(v @ err * spt["singular_values"][k]))
    out["z_combos_s_gt_3"] = ";".join(f"{v:.3g}" for v in z)
    out["combos_covered95"] = int(sum(abs(v) <= 1.96 for v in z)); out["n_combos"] = len(z)
    pe = prediction_errors(prob, best, truth_fit)
    out.update(pe)
    out["rmse_ln_at_truth_residuals"] = float(np.sqrt(np.mean(residuals_jacobian(prob, xt, False)[0] ** 2)))
    return out


# ----------------------------------------------------------------------------------------------- survival estimators (Task 7)
def kaplan_meier(time: np.ndarray, event: np.ndarray, grid=None):
    """Kaplan-Meier estimate with Greenwood variance. Returns (t_unique, S, se). Ties handled; right-censoring only."""
    time = np.asarray(time, float); event = np.asarray(event, int)
    order = np.argsort(time, kind="stable")
    t, e = time[order], event[order]
    ut = np.unique(t[e == 1])
    S, var_sum, out_S, out_se = 1.0, 0.0, [], []
    n = len(t)
    for u in ut:
        at_risk = np.sum(t >= u)
        d = np.sum((t == u) & (e == 1))
        S *= 1 - d / at_risk
        if at_risk > d:
            var_sum += d / (at_risk * (at_risk - d))
        out_S.append(S); out_se.append(S * np.sqrt(var_sum))
    return ut, np.array(out_S), np.array(out_se)


def km_at(ut, S, se, t: float):
    """S(t) and se at a time (step function; 1 before the first event)."""
    j = np.searchsorted(ut, t, side="right") - 1
    return (1.0, 0.0) if j < 0 else (float(S[j]), float(se[j]))


def km_cumulative_incidence_ci(time, event, t: float, conf: float = 0.95):
    """1 - S(t) with log-log (Kalbfleisch-Prentice) CI."""
    from scipy import stats
    ut, S, se = kaplan_meier(time, event)
    s, sd = km_at(ut, S, se, t)
    z = stats.norm.ppf(0.5 + conf / 2)
    if s <= 0 or s >= 1 or sd == 0:
        return 1 - s, max(0.0, 1 - s - z * sd), min(1.0, 1 - s + z * sd)
    se_ll = sd / (s * abs(np.log(s)))
    lo = s ** np.exp(z * se_ll); hi = s ** np.exp(-z * se_ll)
    return 1 - s, 1 - hi, 1 - lo


def turnbull(L: np.ndarray, R: np.ndarray, max_iter: int = 2000, tol: float = 1e-10):
    """Turnbull NPMLE of the survival function for interval-censored data (L, R]; R = inf means right-censored.

    Returns (support points = right ends of the innermost intervals, mass, function F(t)). Standard self-consistency EM on the
    Turnbull intervals (Turnbull 1976); L may be 0 for left-unbounded intervals.
    """
    L = np.asarray(L, float); R = np.asarray(R, float)
    pts = np.unique(np.concatenate([L, R[np.isfinite(R)]]))
    # innermost intervals (q_j, p_j] with q a left end, p a right end and no other endpoint in between
    lefts = np.unique(L); rights = np.unique(R[np.isfinite(R)])
    allp = np.unique(np.concatenate([lefts, rights]))
    typ = {x: ("L" if x in set(lefts) else "") + ("R" if x in set(rights) else "") for x in allp}
    inner = []
    for i in range(len(allp) - 1):
        if "L" in typ[allp[i]] and "R" in typ[allp[i + 1]]:
            inner.append((allp[i], allp[i + 1]))
    if not inner:
        return np.array([]), np.array([]), lambda t: np.zeros_like(np.asarray(t, float))
    A = np.zeros((len(L), len(inner)))
    for j, (q, p) in enumerate(inner):
        A[:, j] = ((L <= q) & (p <= R)).astype(float)
    mass = np.full(len(inner), 1.0 / len(inner))
    for _ in range(max_iter):
        denom = A @ mass
        denom = np.where(denom > 0, denom, 1.0)
        new = (A * mass[None, :] / denom[:, None]).sum(0) / len(L)
        if np.max(np.abs(new - mass)) < tol:
            mass = new
            break
        mass = new
    right_ends = np.array([p for (_, p) in inner])

    def F(t):
        t = np.asarray(t, float)
        return np.array([mass[right_ends <= x].sum() for x in np.atleast_1d(t)]).reshape(t.shape)

    return right_ends, mass, F


def population_T_table(T_of_beta, betas):
    """Monotone tabulation of first-crossing time vs beta (event only for betas below beta_det). Returns (betas_ev, T_ev)."""
    T = np.array([T_of_beta(b) for b in betas])
    ok = np.isfinite(T)
    return np.asarray(betas)[ok], T[ok]


def lognormal_beta_nll(params, betas_ev, T_ev, L, R, beta_det):
    """Negative log-likelihood of ln(beta) ~ N(m, s^2) for interval-censored events (L, R]; T(beta) is increasing in beta.

    F(t) = P(T <= t) = Phi((ln beta*(t) - m)/s) where beta*(t) inverts the tabulated T(beta); beyond the largest event time
    F stays at Phi((ln beta_det - m)/s) (no event for beta above beta_det: the 'cure' mass). R = inf: right-censored.
    """
    from scipy.stats import norm
    m, ls = params
    s = np.exp(ls)
    lt = np.log(betas_ev)

    def F(t):
        t = np.atleast_1d(np.asarray(t, float))
        out = np.empty(len(t))
        for i, x in enumerate(t):
            if not np.isfinite(x):
                out[i] = norm.cdf((np.log(beta_det) - m) / s)
            elif x <= T_ev[0]:
                out[i] = 0.0
            elif x >= T_ev[-1]:
                out[i] = norm.cdf((np.log(beta_det) - m) / s)
            else:
                out[i] = norm.cdf((np.interp(x, T_ev, lt) - m) / s)
        return out

    FL, FR = F(L), F(np.where(np.isfinite(R), R, np.inf))
    ev = np.isfinite(R)
    p = np.where(ev, FR - FL, 1.0 - FL)
    return float(-np.sum(np.log(np.maximum(p, 1e-300))))
