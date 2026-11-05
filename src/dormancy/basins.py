"""Phase 3: basins of attraction of the stable equilibria and their link to simulated relapse events.

A locally stable tumour-free state does NOT guarantee control from every initial condition when other stable
states coexist (Phase 2: tumour-free + nonzero, or two nonzero). This module classifies the attractor reached from
initial conditions (long-horizon integration matched to the enumerated stable equilibria), maps basins on 2-D
slices of the (P, D, E) initial-condition space, and records the 240-month event outcome from the SAME initial
conditions so that attractor and event can be compared (they differ: an initial condition can cross N_detect
and then fall back, or lie in a nonzero basin whose attractor is still below N_detect).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .parameters import Params
from .uncertainty import long_run, parallel_map, regime_of, run_trajectory


def classify_ic(args) -> dict:
    """Top-level worker. args = (p, y0, N_detect, opts, regime). Returns attractor label and 240-month event info."""
    p, y0, N_detect, opts, reg = args
    y0 = np.asarray(y0, float)
    out = {"P0": y0[0], "D0": y0[1], "E0": y0[2], "N0": y0[0] + y0[1]}
    start_above = bool(y0[0] + y0[1] >= N_detect)
    out["start_above"] = start_above
    extra = list(opts.get("extra_thresholds", []))
    traj = run_trajectory(p, y0, sorted(set([N_detect] + extra)), opts.get("horizon", 240.0), method=opts.get("method", "DOP853"),
                          rtol=opts.get("rtol", 1e-9), atol=opts.get("atol", 1e-12), main_threshold=N_detect)
    k = f"{N_detect:g}".replace(".", "p")
    for th in extra:
        kk = f"{th:g}".replace(".", "p")
        out[f"event240_{kk}"] = traj.get(f"event_{kk}", float("nan"))
        out[f"transient240_{kk}"] = traj.get(f"transient_{kk}", float("nan"))
        out[f"start_above_{kk}"] = bool(y0[0] + y0[1] >= th)
    out["event240"] = traj.get(f"event_{k}", float("nan"))
    out["T240"] = traj.get(f"T_{k}", float("nan"))
    out["transient240"] = traj.get(f"transient_{k}", float("nan"))
    out["N240"] = traj.get("N_240", float("nan"))
    lr = long_run(p, y0, N_detect, opts.get("long_floor", 1e5), opts.get("long_tc", 50.0),
                  rtol=opts.get("lr_rtol", 1e-9), atol=opts.get("lr_atol", 1e-13), method=opts.get("lr_method", "LSODA"),
                  regime=reg, with_event=opts.get("with_event_long", False))
    out.update({"attractor": lr["attractor"], "N_inf": lr["N_inf"], "event_long": lr["event_long"],
                "T_event_long": lr["T_event_long"], "T_long": lr["T_long"], "solver_ok": traj["solver_ok"] and lr["long_ok"]})
    return out


def slice_initial_conditions(kind: str, p: Params, n_a: int, n_b: int, lo: float = 1e-6, hi: float = 0.6,
                             fD: float | None = None, e_range=(0.05, 1.6), e_fixed: float = 1.0):
    """Initial-condition grids.

    kind 'PD': axes P0 (a) and D0 (b), both log-spaced in [lo, hi], E0 = e_fixed * E_free.
    kind 'NE': axes N0 (a, log-spaced) and relative immune level e0 (b, linear), dormant fraction fD fixed.
    Returns (a_values, b_values, list_of_y0) with y0 ordered a-major (a varies slowest).
    """
    a = np.geomspace(lo, hi, n_a)
    if kind == "PD":
        b = np.geomspace(lo, hi, n_b)
        ys = [np.array([pa, pb, e_fixed * p.E_free]) for pa in a for pb in b]
    elif kind == "NE":
        if fD is None:
            fD = p.q / (p.q + p.a + p.delta_D)       # equilibrium dormant fraction
        b = np.linspace(*e_range, n_b)
        ys = [np.array([na * (1 - fD), na * fD, eb * p.E_free]) for na in a for eb in b]
    else:
        raise ValueError(kind)
    return a, b, ys


def basin_map(p: Params, ys: list, N_detect: float, opts: dict, workers: int | None = None) -> pd.DataFrame:
    reg = regime_of(p)
    res = parallel_map(classify_ic, [(p, y, N_detect, opts, reg) for y in ys], workers=workers, chunksize=16)
    return pd.DataFrame(res)


def label_codes(labels: pd.Series, order: list[str] | None = None):
    order = order or sorted(set(labels))
    return labels.map({l: i for i, l in enumerate(order)}).values, order


def boundary_mask(L: np.ndarray) -> np.ndarray:
    """True for grid cells having a 4-neighbour with a different label."""
    m = np.zeros_like(L, dtype=bool)
    m[:-1, :] |= L[:-1, :] != L[1:, :]
    m[1:, :] |= L[1:, :] != L[:-1, :]
    m[:, :-1] |= L[:, :-1] != L[:, 1:]
    m[:, 1:] |= L[:, 1:] != L[:, :-1]
    return m


def bisect_boundary(p: Params, ya, yb, N_detect: float, opts: dict, n_iter: int = 40) -> dict:
    """Locate the basin boundary on the path from ya to yb (geometric in P, D; linear in E).

    The two end points must reach different attractors. Returns the path parameter t* at which the label
    changes (bracket width 2^-n_iter), the labels, and the initial state at t*.
    """
    reg = regime_of(p)
    ya, yb = np.asarray(ya, float), np.asarray(yb, float)

    def y_at(t):
        return np.array([ya[0] ** (1 - t) * yb[0] ** t, ya[1] ** (1 - t) * yb[1] ** t, (1 - t) * ya[2] + t * yb[2]])

    la = classify_ic((p, ya, N_detect, opts, reg))["attractor"]
    lb = classify_ic((p, yb, N_detect, opts, reg))["attractor"]
    if la == lb:
        return {"found": False, "label_a": la, "label_b": lb}
    lo, hi = 0.0, 1.0
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        lm = classify_ic((p, y_at(mid), N_detect, opts, reg))["attractor"]
        if lm == la:
            lo = mid
        else:
            hi = mid
    return {"found": True, "label_a": la, "label_b": lb, "t_star": 0.5 * (lo + hi), "bracket": hi - lo, "y_star": y_at(0.5 * (lo + hi))}


def jitter_flip_rate(p: Params, ys: list, base_labels: list, N_detect: float, opts: dict, rel: float = 1e-6,
                     seed: int = 0, workers: int | None = None) -> float:
    """Fraction of initial conditions whose attractor label changes under a relative jitter of size ``rel``."""
    rng = np.random.default_rng(seed)
    ysj = [y * (1.0 + rel * rng.standard_normal(3)) for y in ys]
    df = basin_map(p, ysj, N_detect, opts, workers)
    return float(np.mean(df["attractor"].values != np.asarray(base_labels)))


def random_initial_conditions(p: Params, n: int, seed: int, lo: float = 1e-6, hi: float = 0.6, e_range=(0.05, 1.6)):
    """n random initial states: P0, D0 log-uniform in [lo, hi], E0 = e0 E_free with e0 uniform in e_range."""
    rng = np.random.default_rng(seed)
    P0 = np.exp(rng.uniform(np.log(lo), np.log(hi), n))
    D0 = np.exp(rng.uniform(np.log(lo), np.log(hi), n))
    e0 = rng.uniform(*e_range, n)
    return [np.array([a, b, c * p.E_free]) for a, b, c in zip(P0, D0, e0)]


def tight_options(opts: dict) -> dict:
    """Stricter numerics used to verify labels: Radau, tighter tolerances, 5x longer horizons."""
    t = dict(opts)
    t.update(method="DOP853", rtol=1e-12, atol=1e-15, lr_method="Radau", lr_rtol=1e-11, lr_atol=1e-14,
             long_floor=5 * opts.get("long_floor", 1e5), long_tc=2 * opts.get("long_tc", 50.0))
    return t
