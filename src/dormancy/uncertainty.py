"""Phase 3: uncertainty distributions, sampling, a dense trajectory runner and aggregation of simulated relapse events.

Everything here is a HYPOTHETICAL mathematical construction. Priors are not calibrated to patients; a simulated
event probability is a model-defined quantity, not a patient risk.

Design (see configs/phase3.yaml and report/phase_3_sensitivity_and_uncertainty.md):
  * normalised gauge: NormalizedParams -> Params(s = mu), so the tumour-free immune level is E = e = 1;
  * ONE trajectory runner (`run_trajectory`): a dense, non-terminal integration on [0, H] from which the first
    upward crossing at any number of thresholds, N(60/120/240), downward crossings (transient vs persistent
    exceedance), the dormancy-duration metric and solver diagnostics are derived. It is cross-checked against the
    verified terminal-event solver `dormancy.solver.simulate`;
  * `long_run`: separate long-horizon classifier (first crossing by T_long and the attractor reached), never mixed
    with the 240-month outcome;
  * samples are drawn up-front from a seed and mapped over a process pool with top-level functions and plain
    dicts, so results do not depend on scheduling.
"""
from __future__ import annotations

import math
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import yaml
from scipy import stats
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.stats import qmc

from .equilibria import enumerate_equilibria, tumour_free
from .model import jacobian, rhs
from .normalized import NormalizedParams
from .parameters import Params
from .solver import simulate
from .stability import stability_label

ROOT = Path(__file__).resolve().parents[2]
PARAM_NAMES = ["r", "q", "a", "delta_D", "beta", "rho", "eta", "N0", "fD", "N_detect"]


def load_phase3_config(path: str | Path | None = None) -> dict:
    with open(path or ROOT / "configs" / "phase3.yaml") as fh:
        return yaml.safe_load(fh)


# ----------------------------------------------------------------------------------------------- priors

class Prior:
    """Independent marginals (optionally coupled by a Gaussian copula) with explicit inverse CDFs.

    marginals: dict name -> spec, spec['family'] in {loguniform, uniform, logit, lognormal, beta, fixed}.
    Parameters not listed take their value from ``fixed`` (a dict).
    """

    def __init__(self, marginals: dict, fixed: dict | None = None, copula: list | None = None, name: str = "prior"):
        self.name = name
        self.marginals = dict(marginals)
        self.fixed = dict(fixed or {})
        self.names = [n for n, s in self.marginals.items() if s["family"] != "fixed"]
        for n, s in self.marginals.items():
            if s["family"] == "fixed":
                self.fixed[n] = float(s["value"])
        self.copula = list(copula or [])
        self._validate()
        self.chol = self._cholesky()

    # -- validation ---------------------------------------------------------------------------------
    def _validate(self):
        for n, s in self.marginals.items():
            fam = s["family"]
            if n not in PARAM_NAMES:
                raise ValueError(f"unknown parameter {n!r}")
            if fam in ("loguniform", "uniform", "logit"):
                lo, hi = s["lo"], s["hi"]
                if not (lo < hi):
                    raise ValueError(f"{n}: need lo < hi")
                if fam == "loguniform" and lo <= 0:
                    raise ValueError(f"{n}: log-uniform needs lo > 0")
                if fam == "logit" and not (0 < lo < hi < 1):
                    raise ValueError(f"{n}: logit range must lie in (0,1)")
            elif fam == "lognormal":
                if s["median"] <= 0 or s["factor95"] <= 1:
                    raise ValueError(f"{n}: lognormal needs median > 0 and factor95 > 1")
                if "lo" in s and "hi" in s and not (0 < s["lo"] < s["hi"]):
                    raise ValueError(f"{n}: bad truncation")
            elif fam == "beta":
                if s["alpha"] <= 0 or s["beta"] <= 0:
                    raise ValueError(f"{n}: beta shape parameters must be positive")
            elif fam == "fixed":
                pass
            else:
                raise ValueError(f"{n}: unknown family {fam!r}")
        for a, b, rho in self.copula:
            if a not in self.names or b not in self.names:
                raise ValueError("copula variables must be sampled (non-fixed) parameters")
            if not (-1 < rho < 1):
                raise ValueError("copula correlation must lie in (-1, 1)")

    def _cholesky(self):
        if not self.copula:
            return None
        d = len(self.names)
        C = np.eye(d)
        for a, b, rho in self.copula:
            i, j = self.names.index(a), self.names.index(b)
            C[i, j] = C[j, i] = rho
        return np.linalg.cholesky(C)        # raises if not positive definite

    # -- transforms -----------------------------------------------------------------------------------
    def ppf(self, name: str, u: np.ndarray) -> np.ndarray:
        s = self.marginals[name]
        fam = s["family"]
        u = np.clip(u, 1e-12, 1 - 1e-12)
        if fam == "loguniform":
            return np.exp(np.log(s["lo"]) + u * (np.log(s["hi"]) - np.log(s["lo"])))
        if fam == "uniform":
            return s["lo"] + u * (s["hi"] - s["lo"])
        if fam == "logit":
            L = lambda x: np.log(x / (1 - x))
            return 1.0 / (1.0 + np.exp(-(L(s["lo"]) + u * (L(s["hi"]) - L(s["lo"])))))
        if fam == "lognormal":
            sig = math.log(s["factor95"]) / 1.96 * s.get("sigma_scale", 1.0)
            if "lo" in s and "hi" in s:      # truncated: map u into the CDF interval of the truncation
                cdf = lambda x: stats.norm.cdf(math.log(x / s["median"]) / sig)
                u = cdf(s["lo"]) + u * (cdf(s["hi"]) - cdf(s["lo"]))
            return s["median"] * np.exp(sig * stats.norm.ppf(u))
        if fam == "beta":
            return stats.beta.ppf(u, s["alpha"], s["beta"])
        raise ValueError(fam)

    def from_unit(self, U: np.ndarray) -> list[dict]:
        """Map an (n, d) unit-cube array (columns ordered as ``self.names``) to parameter dicts."""
        U = np.asarray(U, dtype=float)
        if U.shape[1] != len(self.names):
            raise ValueError("unit array has wrong number of columns")
        if self.chol is not None:       # Gaussian copula on the supplied uniforms
            Z = stats.norm.ppf(np.clip(U, 1e-12, 1 - 1e-12)) @ self.chol.T
            U = stats.norm.cdf(Z)
        cols = {n: self.ppf(n, U[:, i]) for i, n in enumerate(self.names)}
        out = []
        for k in range(U.shape[0]):
            row = dict(self.fixed)
            row.update({n: float(cols[n][k]) for n in self.names})
            out.append(row)
        return out

    def sample(self, n: int, seed: int, method: str = "mc") -> list[dict]:
        d = len(self.names)
        if method == "mc":
            U = np.random.default_rng(seed).random((n, d))
        elif method == "lhs":
            U = qmc.LatinHypercube(d=d, seed=seed).random(n)
        elif method == "sobol":
            U = qmc.Sobol(d=d, scramble=True, seed=seed).random(n)
        else:
            raise ValueError(method)
        return self.from_unit(U)


def build_prior(cfg: dict, name: str) -> Prior:
    """Build a named prior from configs/phase3.yaml (supports base / scale_sigma / copula / overrides)."""
    spec = dict(cfg["priors"][name])
    fixed = {k: v for k, v in cfg["fixed"].items() if k in PARAM_NAMES}
    fixed.update({"N_detect": cfg["baseline"]["N_detect"]})
    if "from_global_ranges" in spec:
        marg = {}
        for n in spec["from_global_ranges"]:
            g = cfg["global_ranges"][n]
            marg[n] = {"family": {"log": "loguniform", "logit": "logit", "lin": "uniform"}[g["scale"]], "lo": g["lo"], "hi": g["hi"]}
        return Prior(marg, fixed, name=name)
    marg: dict = {}
    if "base" in spec:
        base = cfg["priors"][spec["base"]]
        marg = {k: dict(v) for k, v in base.items() if isinstance(v, dict) and "family" in v}
    for k, v in spec.items():
        if isinstance(v, dict) and "family" in v:
            marg[k] = dict(v)
    if "scale_sigma" in spec:
        for k, v in marg.items():
            if v["family"] == "lognormal":
                v["sigma_scale"] = spec["scale_sigma"]
    return Prior(marg, fixed, copula=spec.get("copula"), name=name)


def build_global_prior(cfg: dict) -> Prior:
    """Log/logit-uniform ranges over ALL ten inputs (the screening design of the global sensitivity analysis)."""
    marg = {}
    for n, g in cfg["global_ranges"].items():
        marg[n] = {"family": {"log": "loguniform", "logit": "logit", "lin": "uniform"}[g["scale"]], "lo": g["lo"], "hi": g["hi"]}
    fixed = {k: v for k, v in cfg["fixed"].items() if k in PARAM_NAMES}
    return Prior(marg, fixed, name="global_ranges")


# ----------------------------------------------------------------------------------------------- model construction

def row_to_model(row: dict, cfg_fixed: dict | None = None):
    """(Params in the s = mu gauge, y0, N_detect) from a parameter dict."""
    f = {"K": 1.0, "mu": 0.05, "h": 0.02, "e0": 1.0}
    f.update(cfg_fixed or {})
    n = NormalizedParams(K=f["K"], r=row["r"], q=row["q"], a=row["a"], delta_D=row["delta_D"], beta=row["beta"],
                         mu=f["mu"], rho=row["rho"], h=f["h"], eta=row["eta"])
    p = n.to_params()
    e0 = row.get("e0", f["e0"])
    N0, fD = row["N0"], row["fD"]
    y0 = np.array([N0 * (1.0 - fD), N0 * fD, e0 * p.E_free])
    return p, y0, row.get("N_detect", 0.1)


# ----------------------------------------------------------------------------------------------- trajectory runner

def run_trajectory(p: Params, y0, thresholds, horizon: float = 240.0, report_times=(60.0, 120.0, 240.0),
                   method: str = "DOP853", rtol: float = 1e-9, atol: float = 1e-12, dt_grid: float = 0.05,
                   dorm_eps: float = 0.01, dorm_min_frac: float = 0.5, main_threshold: float | None = None) -> dict:
    """Dense, non-terminal integration on [0, horizon]; first upward crossings at all ``thresholds``.

    Event definition = the project's: N = P + D crosses a threshold *upward* for the first time, starting below it.
    A start at or above a threshold has no defined event (``start_above_<i>`` flag, event/time = NaN).
    Exceedance after the first crossing is *transient* if N later crosses the threshold downward within the
    horizon, *persistent* otherwise (N stays >= threshold to the horizon).
    The grid used to bracket crossings has step ``dt_grid`` months (an up-and-down excursion shorter than that
    would be missed; the root is then refined by Brent's method on the dense interpolant).
    """
    y0 = np.asarray(y0, float)
    thresholds = list(thresholds)
    out: dict = {"status": 0, "nfev": 0, "nonfinite": False, "min_state": float("nan"), "E_max": float("nan"),
                 "solver_ok": True}
    try:
        sol = solve_ivp(rhs, (0.0, horizon), y0, args=(p,), method=method, rtol=rtol, atol=atol, dense_output=True)
    except Exception as exc:                                    # noqa: BLE001 - recorded, never silent
        out.update(solver_ok=False, status=-99, error=str(exc))
        return _empty_outputs(out, thresholds, report_times)
    out["status"], out["nfev"] = int(sol.status), int(sol.nfev)
    if sol.status != 0:
        out["solver_ok"] = False
        return _empty_outputs(out, thresholds, report_times)
    tg = np.linspace(0.0, horizon, int(round(horizon / dt_grid)) + 1)
    Y = sol.sol(tg)
    if not np.all(np.isfinite(Y)):
        out.update(nonfinite=True, solver_ok=False)
        return _empty_outputs(out, thresholds, report_times)
    N = Y[0] + Y[1]
    out["min_state"] = float(Y.min())
    out["E_max"] = float(Y[2].max())
    Nfun = lambda t: float(sum(sol.sol(t)[:2]))
    for T in report_times:
        out[f"N_{int(T)}"] = float(Nfun(min(T, horizon)))
    out["N_max"] = float(N.max())
    out["DN_final"] = float(Y[1, -1] / N[-1]) if N[-1] > 0 else float("nan")
    main = main_threshold if main_threshold is not None else thresholds[0]
    for th in thresholds:
        key = _tkey(th)
        out[f"start_above_{key}"] = bool(N[0] >= th)
        if N[0] >= th:
            out[f"event_{key}"] = float("nan"); out[f"T_{key}"] = float("nan")
            out[f"transient_{key}"] = float("nan"); out[f"persistent_{key}"] = float("nan"); out[f"T_down_{key}"] = float("nan")
            continue
        up = np.where((N[:-1] < th) & (N[1:] >= th))[0]
        if up.size == 0:
            out[f"event_{key}"] = 0; out[f"T_{key}"] = float("nan")
            out[f"transient_{key}"] = 0; out[f"persistent_{key}"] = 0; out[f"T_down_{key}"] = float("nan")
            continue
        i = int(up[0])
        t_up = brentq(lambda t: Nfun(t) - th, tg[i], tg[i + 1], xtol=1e-10, rtol=1e-12)
        out[f"event_{key}"] = 1; out[f"T_{key}"] = float(t_up)
        down = np.where((N[:-1] >= th) & (N[1:] < th) & (np.arange(len(N) - 1) > i))[0]
        if down.size:
            j = int(down[0])
            out[f"T_down_{key}"] = float(brentq(lambda t: Nfun(t) - th, tg[j], tg[j + 1], xtol=1e-10, rtol=1e-12))
            out[f"transient_{key}"] = 1; out[f"persistent_{key}"] = 0
        else:
            out[f"T_down_{key}"] = float("nan")
            out[f"transient_{key}"] = 0; out[f"persistent_{key}"] = int(N[-1] >= th)
    # operational dormancy duration (hypothetical definition): time before first crossing (or horizon) with
    # N < threshold, |d ln N/dt| < dorm_eps per month and dormant fraction >= dorm_min_frac
    key = _tkey(main)
    t_end = out.get(f"T_{key}")
    t_end = horizon if (t_end is None or not np.isfinite(t_end)) else t_end
    with np.errstate(divide="ignore", invalid="ignore"):
        dln = np.gradient(np.log(np.maximum(N, 1e-300)), tg)
        dfrac = Y[1] / np.maximum(N, 1e-300)
    mask = (tg <= t_end) & (N < main) & (np.abs(dln) < dorm_eps) & (dfrac >= dorm_min_frac)
    out["dormancy_duration"] = float(mask.sum() * dt_grid)
    out["dormancy_ever_exited"] = bool(np.isfinite(out.get(f"T_{key}", np.nan)))
    return out


def _tkey(th: float) -> str:
    return f"{th:g}".replace(".", "p")


def _empty_outputs(out: dict, thresholds, report_times) -> dict:
    for T in report_times:
        out[f"N_{int(T)}"] = float("nan")
    out.update(N_max=float("nan"), DN_final=float("nan"), dormancy_duration=float("nan"), dormancy_ever_exited=False)
    for th in thresholds:
        k = _tkey(th)
        for pre in ("event_", "T_", "transient_", "persistent_", "T_down_"):
            out[pre + k] = float("nan")
        out[f"start_above_{k}"] = False
    return out


# ----------------------------------------------------------------------------------------------- regimes and long run

def regime_of(p: Params, tol: float = 1e-10) -> dict:
    """Equilibrium structure of a parameter vector (stability from ALL eigenvalues)."""
    tf_eig = np.linalg.eigvals(jacobian(tumour_free(p).state, p))
    tf_stable = bool(np.max(tf_eig.real) < -tol)
    try:
        eqs = enumerate_equilibria(p)
    except ValueError:
        eqs = []
    stable = []
    slowest = [abs(float(np.max(tf_eig.real)))] if tf_stable else []
    for e in eqs:
        ev = np.linalg.eigvals(jacobian(e.state, p))
        if stability_label(ev).startswith("stable"):
            stable.append(e)
            slowest.append(abs(float(np.max(ev.real))))
    n_nz = len(stable)
    if tf_stable and n_nz == 0:
        label = "TF only"
    elif tf_stable and n_nz == 1:
        label = "TF + NZ"
    elif tf_stable and n_nz >= 2:
        label = "TF + NZ + NZ"
    elif (not tf_stable) and n_nz == 1:
        label = "NZ only"
    elif (not tf_stable) and n_nz >= 2:
        label = "NZ + NZ"
    else:
        label = "no stable equilibrium"
    Ns = sorted(e.P + e.D for e in stable)
    return {"regime": label, "tf_stable": tf_stable, "lambda_tf": float(np.max(tf_eig.real)), "n_stable_nonzero": n_nz,
            "stable_N": Ns, "max_stable_N": float(max(Ns)) if Ns else 0.0, "min_stable_N": float(min(Ns)) if Ns else 0.0,
            "n_nonzero_equilibria": len(eqs), "slowest_rate": float(min(slowest)) if slowest else float("nan"),
            "r_eff": p.r_eff}


def long_run(p: Params, y0, N_detect: float, floor: float = 1e5, n_tc: float = 50.0, rtol: float = 1e-9,
             atol: float = 1e-13, method: str = "LSODA", regime: dict | None = None, with_event: bool = True) -> dict:
    """Long-horizon outcome, kept separate from the 240-month outcome.

    * ``T_long = max(floor, n_tc / (slowest stable-equilibrium rate))``;
    * first crossing of N_detect by T_long from the verified ``simulate`` (terminal event);
    * attractor reached: end state of an LSODA run (RHS clipped to the nonnegative orthant ONLY to avoid the
      pole of P/(h+P) at P = -h; the true flow keeps the orthant invariant) matched to the enumerated stable
      equilibria. Label is 'unresolved' unless the final state lies within 1e-3 (relative) of one of them with a
      residual below 1e-8.
    """
    y0 = np.asarray(y0, float)
    reg = regime or regime_of(p)
    rate = reg["slowest_rate"]
    T = max(floor, n_tc / rate) if np.isfinite(rate) and rate > 0 else floor
    T = min(T, 5e6)
    out = {"T_long": float(T), "event_long": float("nan"), "T_event_long": float("nan")}
    if with_event:
        sim = simulate(p, y0, N_detect, T, method="DOP853", rtol=1e-9, atol=1e-12, n_points=10)
        out.update(event_long=int(sim.event), T_event_long=float(sim.time) if sim.event else float("nan"))
    try:
        sol = solve_ivp(lambda t, y, p_: rhs(t, np.maximum(y, 0.0), p_), (0.0, T), y0, args=(p,), method=method,
                        rtol=rtol, atol=atol, dense_output=True)
        yf = sol.y[:, -1]
        ok = sol.status == 0 and np.all(np.isfinite(yf))
        if ok:     # maximum burden along the whole long run (lets any threshold's long-horizon event be read off)
            tt = np.concatenate([np.linspace(0.0, min(240.0, T), 481), np.geomspace(min(240.0, T) + 1e-9, T, 3000)])
            YY = np.maximum(sol.sol(tt), 0.0)
            Nn = YY[0] + YY[1]
            out["N_max_long"] = float(Nn.max())
            out["t_N_max_long"] = float(tt[int(np.argmax(Nn))])
    except Exception:                                           # noqa: BLE001
        ok, yf = False, np.full(3, np.nan)
    label, Ninf = "unresolved", float("nan")
    if ok:
        res = float(np.max(np.abs(rhs(0.0, np.maximum(yf, 0.0), p))))
        cands = []
        if reg["tf_stable"]:
            cands.append(("TF", tumour_free(p).state))
        stab = [e for e in (enumerate_equilibria(p) if reg["n_stable_nonzero"] else [])
                if stability_label(np.linalg.eigvals(jacobian(e.state, p))).startswith("stable")]
        stab.sort(key=lambda e: e.P + e.D)
        names = ["NZ"] if len(stab) == 1 else ["NZ-low", "NZ-high", "NZ-3"][:len(stab)]
        cands += [(nm, e.state) for nm, e in zip(names, stab)]
        for nm, x in cands:
            if np.linalg.norm(yf - x) <= 1e-3 * max(np.linalg.norm(x), 1e-3) + 1e-9 and res < 1e-8:
                label, Ninf = nm, float(x[0] + x[1])
                break
        out["residual_end"] = res
    out.update(attractor=label, N_inf=Ninf, long_ok=bool(ok), end_N=float(yf[0] + yf[1]) if ok else float("nan"))
    return out


# ----------------------------------------------------------------------------------------------- one-sample pipeline

def evaluate_row(args) -> dict:
    """Top-level (picklable) worker: one parameter dict -> all outputs. ``args`` = (row, options)."""
    row, opt = args
    thresholds = opt["thresholds"]
    main = row.get("N_detect", opt.get("main_threshold", 0.1))
    ths = sorted(set(thresholds) | {main}) if opt.get("add_main", True) else list(thresholds)
    p, y0, _ = row_to_model(row, opt.get("fixed"))
    res = dict(row)
    t = run_trajectory(p, y0, ths, opt["horizon"], opt.get("report_times", (60.0, 120.0, 240.0)),
                       opt.get("method", "DOP853"), opt.get("rtol", 1e-9), opt.get("atol", 1e-12),
                       main_threshold=main, dorm_eps=opt.get("dorm_eps", 0.01), dorm_min_frac=opt.get("dorm_min_frac", 0.5))
    res.update(t)
    mk = _tkey(main)
    res["event_main"] = t.get(f"event_{mk}", float("nan"))
    res["T_main"] = t.get(f"T_{mk}", float("nan"))
    res["transient_main"] = t.get(f"transient_{mk}", float("nan"))
    res["persistent_main"] = t.get(f"persistent_{mk}", float("nan"))
    res["start_above_main"] = t.get(f"start_above_{mk}", False)
    if opt.get("check_against_simulate", True) and not res["start_above_main"] and t["solver_ok"]:
        s = simulate(p, y0, main, opt["horizon"], opt.get("method", "DOP853"), opt.get("rtol", 1e-9), opt.get("atol", 1e-12), n_points=10)
        res["simulate_event"] = int(s.event)
        res["simulate_dT"] = abs(s.time - res["T_main"]) if (s.event and res["event_main"] == 1) else (0.0 if s.event == res["event_main"] else float("nan"))
        res["simulate_agrees"] = bool(s.event == res["event_main"] and (not s.event or abs(s.time - res["T_main"]) < 1e-4))
    reg = regime_of(p)
    res.update({f"reg_{k}": v for k, v in reg.items() if k not in ("stable_N",)})
    res["stable_N_list"] = ";".join(f"{x:.6g}" for x in reg["stable_N"])
    res["beta_over_betac"] = row["beta"] / max(p.r_eff, 1e-12) if p.r_eff > 0 else float("inf")
    res["r_eff"] = p.r_eff
    if opt.get("long_run", False):
        lr = long_run(p, y0, main, opt.get("long_floor", 1e5), opt.get("long_tc", 50.0), regime=reg)
        res.update(lr)
    return res


def parallel_map(fn, items, workers: int | None = None, chunksize: int = 32):
    """Ordered process-pool map; results are independent of scheduling. workers=1 runs in-process."""
    items = list(items)
    workers = workers if workers is not None else min(8, os.cpu_count() or 1)
    if workers <= 1 or len(items) < 2 * workers:
        return [fn(x) for x in items]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(fn, items, chunksize=chunksize))


def default_options(cfg: dict, **over) -> dict:
    opt = {"thresholds": list(cfg["thresholds"]), "horizon": cfg["horizon_months"], "fixed": cfg["fixed"],
           "method": cfg["solver"]["method"], "rtol": cfg["solver"]["rtol"], "atol": cfg["solver"]["atol"],
           "dorm_eps": cfg["dormancy_definition"]["max_abs_dlnN_dt_per_month"],
           "dorm_min_frac": cfg["dormancy_definition"]["min_dormant_fraction"],
           "long_floor": cfg["long_run"]["horizon_floor_months"], "long_tc": cfg["long_run"]["n_time_constants"]}
    opt.update(over)
    return opt


# ----------------------------------------------------------------------------------------------- aggregation

def wilson_interval(k: int, n: int, conf: float = 0.95) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    z = stats.norm.ppf(0.5 + conf / 2)
    ph = k / n
    den = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / den
    h = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return (max(0.0, min(ph, c - h)), min(1.0, max(ph, c + h)))      # clamp: interval always contains the estimate


def clopper_pearson(k: int, n: int, conf: float = 0.95) -> tuple[float, float]:
    if n == 0:
        return (float("nan"), float("nan"))
    a = 1 - conf
    lo = 0.0 if k == 0 else stats.beta.ppf(a / 2, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(1 - a / 2, k + 1, n - k)
    return (float(lo), float(hi))


def upper_bound_rare(k: int, n: int, conf: float = 0.95) -> float:
    """One-sided Clopper-Pearson upper bound (use when k is 0 or very small)."""
    return 1.0 if k == n else float(stats.beta.ppf(conf, k + 1, n - k))


def event_summary(T: np.ndarray, event: np.ndarray, horizons=(60, 120, 240), horizon: float = 240.0) -> dict:
    """Simulated cumulative event probabilities and summaries.

    T: first-crossing times (NaN when no event); event: 0/1. Only horizon censoring exists, so the Kaplan-Meier
    estimate equals the empirical CDF and RMST(tau) = mean(min(T, tau)) exactly.
    The unconditional median is returned only if the probability of an event by the horizon is >= 0.5.
    """
    T = np.asarray(T, float)
    event = np.asarray(event, float)
    ok = np.isfinite(event)
    T, event = T[ok], event[ok]
    n = len(event)
    out = {"n": int(n)}
    for h in horizons:
        k = int(np.sum((event == 1) & (T <= h)))
        lo, hi = wilson_interval(k, n)
        clo, chi = clopper_pearson(k, n)
        out[f"P_by_{int(h)}"] = k / n if n else float("nan")
        out[f"P_by_{int(h)}_wilson_lo"], out[f"P_by_{int(h)}_wilson_hi"] = lo, hi
        out[f"P_by_{int(h)}_cp_lo"], out[f"P_by_{int(h)}_cp_hi"] = clo, chi
        out[f"P_by_{int(h)}_upper95_onesided"] = upper_bound_rare(k, n)
        out[f"count_by_{int(h)}"] = k
        tau = np.where(event == 1, np.minimum(T, h), h)
        out[f"RMST_{int(h)}"] = float(np.mean(tau))
        out[f"RMST_{int(h)}_se"] = float(np.std(tau, ddof=1) / math.sqrt(n)) if n > 1 else float("nan")
    ev_t = T[(event == 1) & (T <= horizon)]
    out["n_events"] = int(len(ev_t))
    out["conditional_median_T"] = float(np.median(ev_t)) if len(ev_t) else float("nan")
    out["conditional_q25_T"] = float(np.quantile(ev_t, 0.25)) if len(ev_t) else float("nan")
    out["conditional_q75_T"] = float(np.quantile(ev_t, 0.75)) if len(ev_t) else float("nan")
    pH = out[f"P_by_{int(horizons[-1])}"]
    out["unconditional_median_T"] = float(np.median(np.where(event == 1, T, np.inf))) if pH >= 0.5 else float("nan")
    return out
