"""Phase 5 (descriptive): model-generated event-time distributions versus the METABRIC RFS Kaplan-Meier curve.

Nothing here fits, optimises or selects a model; every function is a pure transformation or a descriptive discrepancy metric.
The design (scenario grid, horizons, metrics, consistency criterion) is pre-registered in configs/phase5.yaml and its sha256 is
recorded in results/tables/p5_prereg_hash.txt BEFORE any comparison is computed (`verify_prereg`).

Time unit: months everywhere. Mapping (Phase 4 `endpoint_mapping`): t_clin = tau + t_model + d (recurrence),
t_clin = tau + t_model + Lag (death-first). Model events are first upward crossings of a model-defined burden threshold.
"""
from __future__ import annotations

import hashlib
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .endpoint_mapping import ClinicalMapping, cumulative_incidence, dkw_epsilon, map_event_times
from .km_descriptive import km_point
from .synthetic_validation import kaplan_meier


# ------------------------------------------------------------------------------------------------ pre-registration
def sha256_file(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_prereg(path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def parse_hash_record(path) -> dict:
    rec = {}
    for line in Path(path).read_text().splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            rec[k.strip()] = v.strip()
    return rec


def verify_prereg(cfg_path, hash_path) -> str:
    """Raise unless the design file still matches the hash recorded before the comparison. Returns the hash."""
    rec = parse_hash_record(hash_path)
    h = sha256_file(cfg_path)
    if rec.get("sha256") != h:
        raise RuntimeError(f"pre-registered design changed: recorded {rec.get('sha256')} != current {h}")
    return h


def prereg_precedes(hash_path, *output_paths) -> bool:
    """True if the hash record file was last written no later than every given comparison output (ordering logic)."""
    t0 = Path(hash_path).stat().st_mtime
    return all(Path(p).stat().st_mtime >= t0 for p in output_paths)


# ------------------------------------------------------------------------------------------------ simulated populations
def tkey(th: float) -> str:
    return f"{th:g}".replace(".", "p")


def load_prior_times(table_dir, files, threshold: float) -> np.ndarray:
    """Model first-crossing times (months) for one prior at one threshold; np.inf = no event within 240 months (censored).

    Samples that start above the threshold have no defined event in Phase 3 and are excluded (count was 0 in every run)."""
    k = tkey(threshold)
    out = []
    for f in files:
        d = pd.read_csv(Path(table_dir) / f, comment="#")
        keep = d[f"start_above_{k}"].astype(int).to_numpy() == 0
        ev = d[f"event_{k}"].to_numpy().astype(int) == 1
        T = np.where(ev, d[f"T_{k}"].to_numpy(float), np.inf)
        out.append(T[keep])
    return np.concatenate(out)


def map_times(T, tau=0.0, d=0.0, death_first=0.0, lag=24.0, seed=0) -> np.ndarray:
    m = ClinicalMapping(tau=float(tau), detect_delay=float(d), death_first_fraction=float(death_first), death_lag_mean=float(lag))
    return map_event_times(T, m, np.random.default_rng(seed))


def ecdf_at(T_clin, t) -> np.ndarray:
    """Empirical event-free complement: F(t) = P(T_clin <= t); inf = no event (what Kaplan-Meier equals with only administrative censoring)."""
    return cumulative_incidence(T_clin, np.atleast_1d(np.asarray(t, float)))


def shifted_cdf(T_model, grid, shift) -> np.ndarray:
    """F_model(t - shift) on a clinical grid for any real shift (negative = earlier); t - shift < 0 gives 0."""
    return cumulative_incidence(T_model, np.asarray(grid, float) - shift)


def wilson(k, n, z=1.959964):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, c - h), min(1.0, c + h)


# ------------------------------------------------------------------------------------------------ METABRIC reference
def load_metabric(path):
    pat = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    pat["RFS_MONTHS"] = pd.to_numeric(pat["RFS_MONTHS"], errors="coerce")
    ev = pat.RFS_STATUS.str.startswith("1").astype(float).where(pat.RFS_STATUS.notna())
    d = pat.assign(event=ev).dropna(subset=["RFS_MONTHS", "event"]).copy()
    d["event"] = d["event"].astype(int)
    return d["RFS_MONTHS"].to_numpy(float), d["event"].to_numpy(int)


def km_cuminc_grid(time, event, grid) -> np.ndarray:
    ut, S, _ = kaplan_meier(np.asarray(time, float), np.asarray(event, int))
    idx = np.searchsorted(ut, np.asarray(grid, float), side="right")
    Sg = np.concatenate([[1.0], S])[idx]
    return 1.0 - Sg


def km_reference(time, event, horizons) -> dict:
    return {int(h): km_point(time, event, float(h)) for h in horizons}


# ------------------------------------------------------------------------------------------------ discrepancy metrics
def signed_difference(F_sim, F_km) -> float:
    return float(F_sim - F_km)


def inside_ci(F_sim, lo, hi) -> bool:
    return bool(lo <= F_sim <= hi)


def sup_norm(F_a, F_b) -> float:
    return float(np.max(np.abs(np.asarray(F_a) - np.asarray(F_b))))


def late_increment(F120, F240) -> float:
    return float(F240 - F120)


def conditional_late_risk(F120, F240) -> float:
    return float((F240 - F120) / (1.0 - F120)) if F120 < 1 else float("nan")


def ratio_240_120(F120, F240) -> float:
    return float(F240 / F120) if F120 > 0 else float("inf")


def scenario_grid(cfg) -> list[dict]:
    g = cfg["grid"]
    rows = []
    for i, (p, th, tau, d, df) in enumerate(itertools.product(g["priors"], g["thresholds"], g["tau_months"],
                                                              g["detect_delay_months"], g["death_first_fraction"])):
        rows.append({"scenario_id": i, "prior": p, "threshold": th, "tau": tau, "d": d, "death_first": df})
    return rows


def evaluate_scenario(sc, T_model, km, km_grid_F, cfg, grid_t) -> dict:
    """All pre-registered metrics for one scenario. `km` = km_reference dict, `km_grid_F` = KM F on `grid_t` (1-month grid)."""
    seed = int(cfg["simulation_source"]["seed_for_stochastic_mapping"]) + int(sc["scenario_id"])
    Tc = map_times(T_model, sc["tau"], sc["d"], sc["death_first"], cfg["grid"]["death_lag_mean_months"], seed)
    hs = cfg["horizons_months"]["primary"] + cfg["horizons_months"]["secondary"]
    row = dict(sc)
    row["n_sim"] = len(Tc)
    Fs = {}
    for h in hs:
        F = float(ecdf_at(Tc, h)[0])
        Fs[h] = F
        k = km[h]
        row[f"F_sim_{h}"] = F
        row[f"F_km_{h}"] = k["cum_inc"]
        row[f"D_{h}"] = signed_difference(F, k["cum_inc"])
        row[f"inside_ci_{h}"] = inside_ci(F, k["ci_lo"], k["ci_hi"])
        row[f"below_ci_{h}"] = F < k["ci_lo"]
        row[f"above_ci_{h}"] = F > k["ci_hi"]
    Fg = ecdf_at(Tc, grid_t)
    for h in (120, 240):
        m = grid_t <= h
        row[f"sup_norm_{h}"] = sup_norm(Fg[m], km_grid_F[m])
    li_s, li_k = late_increment(Fs[120], Fs[240]), late_increment(km[120]["cum_inc"], km[240]["cum_inc"])
    row["late_inc_sim"], row["late_inc_km"] = li_s, li_k
    row["cond_late_risk_sim"] = conditional_late_risk(Fs[120], Fs[240])
    row["cond_late_risk_km"] = conditional_late_risk(km[120]["cum_inc"], km[240]["cum_inc"])
    row["ratio_240_120_sim"] = ratio_240_120(Fs[120], Fs[240])
    row["ratio_240_120_km"] = ratio_240_120(km[120]["cum_inc"], km[240]["cum_inc"])
    row["frac_after_120_sim"] = (Fs[240] - Fs[120]) / Fs[240] if Fs[240] > 0 else float("nan")
    row["frac_after_120_km"] = (km[240]["cum_inc"] - km[120]["cum_inc"]) / km[240]["cum_inc"]
    # pre-specified statements
    row["S1"] = bool(row["below_ci_60"])
    row["S2"] = bool(row["below_ci_120"])
    row["S3"] = bool(row["inside_ci_240"])
    row["S4"] = bool(li_s > li_k)
    row["S5"] = bool(row["frac_after_120_sim"] > row["frac_after_120_km"])
    row["S6"] = bool(row["ratio_240_120_sim"] > 2 and row["ratio_240_120_km"] < 2)
    row["exploratory"] = False
    row["_Fg"] = Fg      # curve on the 1-month grid (popped by the caller before writing tables)
    return row
