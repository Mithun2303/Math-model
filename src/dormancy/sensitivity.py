"""Phase 3: local and global sensitivity estimators.

Sensitivity is a property of the model over the chosen *ranges*; it is not uncertainty about the parameters
(see uncertainty.py for propagated uncertainty). All ranges are hypothetical.

Local:   central differences in log-parameter at two steps (h, h/2); derivative reported only where both agree
         and every stencil point defines the output; event indicators are never differentiated (flip margin).
Global:  Spearman / partial rank correlations (+ bootstrap), Morris elementary effects, Sobol indices with the
         Saltelli (2010) first-order and Jansen total-order estimators and Saltelli (2002) second-order indices.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import stats
from scipy.stats import qmc

# ----------------------------------------------------------------------------------------------- local

def perturb_row(row: dict, name: str, log_step: float) -> dict:
    """Multiply input ``name`` by exp(log_step). Supports the pseudo-inputs P0 and D0 (initial burdens)."""
    r = dict(row)
    fac = math.exp(log_step)
    if name == "P0":
        P0 = row["N0"] * (1 - row["fD"]) * fac
        D0 = row["N0"] * row["fD"]
        r["N0"], r["fD"] = P0 + D0, D0 / (P0 + D0)
    elif name == "D0":
        P0 = row["N0"] * (1 - row["fD"])
        D0 = row["N0"] * row["fD"] * fac
        r["N0"], r["fD"] = P0 + D0, D0 / (P0 + D0)
    else:
        r[name] = row[name] * fac
    return r


def log_derivative(values: dict, h: float, log_output: bool, tol: float = 0.2, tiny: float = 1e-12):
    """Central difference d out / d ln(theta) at steps h and h/2.

    values: {+h, -h, +h/2, -h/2 -> output}. Returns (derivative, status). Status:
      'ok'                    both estimates defined and agree within ``tol`` (relative);
      'undefined'             some stencil output is missing/NaN (e.g. event absent at a stencil point);
      'unstable'              the two step sizes disagree (near a discontinuity or a switching boundary);
      'zero'                  all stencil values identical (no dependence).
    For ``log_output`` the output is log-transformed (an elasticity d ln out / d ln theta).
    """
    keys = [h, -h, h / 2, -h / 2]
    v = [values.get(k, float("nan")) for k in keys]
    if not all(np.isfinite(x) for x in v):
        return float("nan"), "undefined"
    if log_output:
        if any(x <= 0 for x in v):
            return float("nan"), "undefined"
        v = [math.log(x) for x in v]
    if max(v) - min(v) <= tiny * max(1.0, max(abs(x) for x in v)):
        return 0.0, "zero"
    d1 = (v[0] - v[1]) / (2 * h)
    d2 = (v[2] - v[3]) / h
    if abs(d1 - d2) > tol * max(abs(d1), abs(d2), tiny):
        return float("nan"), "unstable"
    return d2, "ok"                      # report the finer step (h/2)


def flip_margin(indicator, row: dict, name: str, base_value: int, max_log: float = math.log(10.0),
                n_scan: int = 120, n_bisect: int = 60) -> dict:
    """Smallest change of input ``name`` (multiplicative) that flips a 0/1 indicator, searched up and down.

    indicator(row) -> 0/1 (or NaN if undefined). Returns, per direction, the factor at which the indicator first
    differs from ``base_value`` (refined by bisection), or inf if none within a factor ``exp(max_log)``.
    """
    out = {}
    for sgn, label in ((+1, "up"), (-1, "down")):
        prev = 0.0
        found = None
        for s in np.linspace(max_log / n_scan, max_log, n_scan):
            val = indicator(perturb_row(row, name, sgn * s))
            if np.isfinite(val) and int(val) != base_value:
                lo, hi = prev, s
                for _ in range(n_bisect):
                    mid = 0.5 * (lo + hi)
                    vm = indicator(perturb_row(row, name, sgn * mid))
                    if np.isfinite(vm) and int(vm) != base_value:
                        hi = mid
                    else:
                        lo = mid
                found = hi
                break
            prev = s
        out[f"{label}_factor"] = math.exp(sgn * found) if found is not None else float("inf")
    return out


# ----------------------------------------------------------------------------------------------- rank statistics

def _rank(X: np.ndarray) -> np.ndarray:
    return np.apply_along_axis(stats.rankdata, 0, X)


def spearman_columns(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Spearman correlation of each column of X with y (ties get average ranks)."""
    Rx = _rank(X)
    ry = stats.rankdata(y)
    Rx = Rx - Rx.mean(0)
    ry = ry - ry.mean()
    den = np.sqrt((Rx ** 2).sum(0) * (ry ** 2).sum())
    return (Rx * ry[:, None]).sum(0) / np.where(den > 0, den, np.nan)


def prcc(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Partial rank correlation of y with each column of X, controlling for the other columns.

    Uses the inverse of the rank correlation matrix of [X, y]: r_{iy.rest} = -P_iy / sqrt(P_ii P_yy).
    Assumes monotone relations; it is not meaningful for strongly non-monotone (e.g. bistable) responses.
    """
    Z = np.column_stack([_rank(X), stats.rankdata(y)])
    C = np.corrcoef(Z, rowvar=False)
    if not np.all(np.isfinite(C)):
        return np.full(X.shape[1], np.nan)
    P = np.linalg.pinv(C)
    d = X.shape[1]
    return -P[:d, d] / np.sqrt(P[:d, :d].diagonal() * P[d, d])


def bootstrap_ci(stat_fn, X: np.ndarray, y: np.ndarray, B: int = 500, seed: int = 0, conf: float = 0.95):
    """Percentile bootstrap (resampling rows) for a vector-valued statistic stat_fn(X, y)."""
    rng = np.random.default_rng(seed)
    n = len(y)
    reps = np.empty((B, X.shape[1]))
    for b in range(B):
        idx = rng.integers(0, n, n)
        reps[b] = stat_fn(X[idx], y[idx])
    a = (1 - conf) / 2
    return np.nanquantile(reps, a, axis=0), np.nanquantile(reps, 1 - a, axis=0)


def kendall_tau(a: np.ndarray, b: np.ndarray) -> float:
    return float(stats.kendalltau(a, b)[0])


# ----------------------------------------------------------------------------------------------- Morris

def morris_design(d: int, r: int, levels: int = 4, seed: int = 0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Morris trajectories in the unit cube. Returns (X, order, sign): X has shape (r, d+1, d)."""
    rng = np.random.default_rng(seed)
    delta = levels / (2.0 * (levels - 1))
    grid = np.arange(levels) / (levels - 1)
    base_levels = grid[grid <= 1 - delta + 1e-12]
    X = np.empty((r, d + 1, d))
    orders = np.empty((r, d), int)
    signs = np.empty((r, d), int)
    for t in range(r):
        x = rng.choice(base_levels, d)
        order = rng.permutation(d)
        sign = np.where(x + delta <= 1 + 1e-12, 1, -1)
        sign = np.where(rng.random(d) < 0.5, sign, np.where(x - delta >= -1e-12, -1, 1))
        X[t, 0] = x
        for k, i in enumerate(order):
            x = x.copy()
            x[i] = x[i] + sign[i] * delta
            X[t, k + 1] = np.clip(x, 0, 1)
        orders[t], signs[t] = order, sign
    return X, orders, signs


def morris_effects(Y: np.ndarray, orders: np.ndarray, signs: np.ndarray, delta: float) -> dict:
    """Elementary effects from outputs Y of shape (r, d+1); returns mu, mu_star, sigma per input."""
    r, d1 = Y.shape
    d = d1 - 1
    ee = np.full((r, d), np.nan)
    for t in range(r):
        for k, i in enumerate(orders[t]):
            ee[t, i] = (Y[t, k + 1] - Y[t, k]) / (signs[t, i] * delta)
    return {"mu": np.nanmean(ee, 0), "mu_star": np.nanmean(np.abs(ee), 0), "sigma": np.nanstd(ee, 0, ddof=1), "ee": ee}


# ----------------------------------------------------------------------------------------------- Sobol / Saltelli

def saltelli_unit_matrices(d: int, N: int, seed: int, second_order: bool = True):
    """Scrambled-Sobol base matrices A, B (N, d) from one 2d-dimensional sequence, and AB_i, BA_i lists.

    N should be a power of 2 (balance properties). The first n rows (n a power of 2) form a valid smaller design.
    """
    U = qmc.Sobol(d=2 * d, scramble=True, seed=seed).random(N)
    A, B = U[:, :d], U[:, d:]
    AB = []
    BA = []
    for i in range(d):
        M = A.copy(); M[:, i] = B[:, i]; AB.append(M)
        if second_order:
            M2 = B.copy(); M2[:, i] = A[:, i]; BA.append(M2)
    return A, B, AB, BA


def sobol_from_outputs(fA, fB, fAB, fBA=None) -> dict:
    """Sobol indices from model outputs: first-order (Saltelli 2010), total (Jansen 1999), second-order (Saltelli 2002)."""
    fA, fB = np.asarray(fA, float), np.asarray(fB, float)
    fAB = np.asarray(fAB, float)            # shape (d, N)
    d = fAB.shape[0]
    V = np.var(np.concatenate([fA, fB]), ddof=1)
    if not V > 0:
        return {"S1": np.full(d, np.nan), "ST": np.full(d, np.nan), "S2": np.full((d, d), np.nan), "V": V}
    S1 = np.array([np.mean(fB * (fAB[i] - fA)) / V for i in range(d)])
    ST = np.array([0.5 * np.mean((fA - fAB[i]) ** 2) / V for i in range(d)])
    S2 = np.full((d, d), np.nan)
    if fBA is not None:
        fBA = np.asarray(fBA, float)
        for j in range(d):
            for k in range(j + 1, d):
                Vjk = np.mean(fBA[j] * fAB[k] - fA * fB)
                S2[j, k] = S2[k, j] = Vjk / V - S1[j] - S1[k]
    return {"S1": S1, "ST": ST, "S2": S2, "V": V}


def sobol_bootstrap(fA, fB, fAB, fBA=None, B: int = 300, seed: int = 0, conf: float = 0.95) -> dict:
    """Bootstrap (resampling the N base rows jointly) percentile intervals for S1, ST (and S2)."""
    rng = np.random.default_rng(seed)
    fA, fB, fAB = np.asarray(fA), np.asarray(fB), np.asarray(fAB)
    fBA = None if fBA is None else np.asarray(fBA)
    N = len(fA)
    S1, ST, S2 = [], [], []
    for _ in range(B):
        idx = rng.integers(0, N, N)
        res = sobol_from_outputs(fA[idx], fB[idx], fAB[:, idx], None if fBA is None else fBA[:, idx])
        S1.append(res["S1"]); ST.append(res["ST"]); S2.append(res["S2"])
    a = (1 - conf) / 2
    q = lambda arr: (np.nanquantile(arr, a, axis=0), np.nanquantile(arr, 1 - a, axis=0))
    return {"S1": q(np.array(S1)), "ST": q(np.array(ST)), "S2": q(np.array(S2)) if fBA is not None else None}


def ishigami(x: np.ndarray, a: float = 7.0, b: float = 0.1) -> np.ndarray:
    """Ishigami function on [-pi, pi]^3 (x shape (n, 3)); analytic Sobol indices:
    S1 = (0.3139, 0.4424, 0), ST = (0.5576, 0.4424, 0.2437), S13 = 0.2437, S12 = S23 = 0."""
    return np.sin(x[:, 0]) + a * np.sin(x[:, 1]) ** 2 + b * x[:, 2] ** 4 * np.sin(x[:, 0])
