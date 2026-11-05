"""Phase 3, Task 7: robustness of simulated relapse outcomes to the (model-defined) detection threshold N_detect.

Run AFTER monte_carlo_relapse.py:  .venv/bin/python experiments/threshold_robustness.py
Uses common random numbers: the same 5000 prior-P1 samples and the same dense trajectories for every threshold.
The thresholds are model-defined burden fractions of K and are NOT clinical imaging limits.
Outputs: results/tables/p3_threshold_*.csv, results/figures/fig_p3_09_detection_threshold_sensitivity.png
"""
import numpy as np

from _p3common import *  # noqa: F401,F403
import matplotlib.pyplot as plt

from dormancy.equilibria import enumerate_equilibria
from dormancy.model import jacobian
from dormancy.normalized import NormalizedParams
from dormancy.parameters import Params
from dormancy.stability import stability_label
from dormancy.uncertainty import (event_summary, long_run, parallel_map, regime_of, row_to_model, wilson_interval)

THS = CFG["thresholds"]
HZ = CFG["report_horizons_months"]
GRID = np.geomspace(0.05, 20.0, 80)                    # multiplicative factors applied to beta in the boundary scan


def tkey(th):
    return f"{th:g}".replace(".", "p")


def scan_beta(args):
    """For one sample: long-horizon max burden and equilibrium predictions as beta is scaled by GRID (all else fixed)."""
    row, thetas = args
    out = {"Nmax": [], "min_stable_N": [], "max_stable_N": [], "tf_stable": []}
    for g in GRID:
        r = dict(row); r["beta"] = row["beta"] * g
        p, y0, _ = row_to_model(r, CFG["fixed"])
        reg = regime_of(p)
        lr = long_run(p, y0, 0.1, with_event=False, regime=reg)
        out["Nmax"].append(lr.get("N_max_long", float("nan")))
        out["min_stable_N"].append(reg["min_stable_N"] if reg["n_stable_nonzero"] else 0.0)
        out["max_stable_N"].append(reg["max_stable_N"])
        out["tf_stable"].append(reg["tf_stable"])
    out = {k: np.array(v) for k, v in out.items()}
    res = {"r_eff": Params().r_eff}
    p0, _, _ = row_to_model(row, CFG["fixed"])
    res["r_eff"] = p0.r_eff
    res["beta"] = row["beta"]
    for th in thetas:
        k = tkey(th)
        sim_ev = out["Nmax"] >= th
        for nm, pred in (("sim", sim_ev), ("eq_min_stable", out["min_stable_N"] >= th), ("eq_max_stable", out["max_stable_N"] >= th)):
            # boundary = first beta (increasing) at which the indicator is False after being True at the smallest beta
            res[f"beta_det_{nm}_{k}"] = _boundary(row["beta"] * GRID, pred)
        res[f"reentrant_sim_{k}"] = bool(_reentrant(sim_ev))
        res[f"agree_eq_min_{k}"] = float(np.mean(sim_ev == (out["min_stable_N"] >= th)))
        res[f"agree_eq_max_{k}"] = float(np.mean(sim_ev == (out["max_stable_N"] >= th)))
        res[f"agree_kc_{k}"] = float(np.mean(sim_ev == ~out["tf_stable"]))
    return res


def _boundary(betas, pred):
    pred = np.asarray(pred, bool)
    if not pred[0]:
        return float("nan")                       # no event even at the smallest beta of the scan
    idx = np.where(~pred)[0]
    if len(idx) == 0:
        return float("inf")                       # event at every beta of the scan
    i = idx[0]
    return float(np.sqrt(betas[i - 1] * betas[i]))


def _reentrant(sim_ev):
    sim_ev = np.asarray(sim_ev, bool)
    d = np.diff(sim_ev.astype(int))
    return bool((d == 1).any())                   # False -> True as beta increases: event re-appears (non-monotone in beta)


def main():
    D = pd.read_csv(TAB / "p3_mc_P1_samples.csv", comment="#")
    n = len(D)
    log(f"loaded {n} samples (prior P1_base); start-above counts per threshold: " +
        str({th: int(D[f'start_above_{tkey(th)}'].sum()) for th in THS if f'start_above_{tkey(th)}' in D}))

    # ----------------------------------------------------------------- 1. event probabilities and exceedance types vs threshold
    rows = []
    for th in THS:
        k = tkey(th)
        valid = ~D[f"start_above_{k}"].astype(bool)
        Dv = D[valid]
        es = event_summary(Dv[f"T_{k}"].values, Dv[f"event_{k}"].values, HZ, CFG["horizon_months"])
        ev = Dv[f"event_{k}"] == 1
        row = {"threshold": th, "n_valid": len(Dv), **{f: es[f] for f in
               ("P_by_60", "P_by_60_wilson_lo", "P_by_60_wilson_hi", "P_by_120", "P_by_120_wilson_lo", "P_by_120_wilson_hi",
                "P_by_240", "P_by_240_wilson_lo", "P_by_240_wilson_hi", "conditional_median_T", "conditional_q25_T",
                "conditional_q75_T", "RMST_240", "n_events", "count_by_60", "P_by_60_upper95_onesided")}}
        row["frac_events_transient_within_240"] = float(Dv.loc[ev, f"transient_{k}"].mean()) if ev.any() else float("nan")
        row["frac_events_persistent_at_240"] = float(Dv.loc[ev, f"persistent_{k}"].mean()) if ev.any() else float("nan")
        # long-horizon event: the trajectory ever reaches theta (max burden along the long run)
        row["P_long_horizon_event"] = float((Dv.N_max_long >= th).mean())
        row["frac_long_event_not_seen_by_240"] = float(((Dv.N_max_long >= th) & ~ev).mean())
        # equilibrium-based predictions at this threshold (from the SAME samples, N_inf = attractor reached)
        obs24, obsL = ev.values, (Dv.N_max_long >= th).values
        for nm, pv in (("TF_unstable", ~Dv.reg_tf_stable.astype(bool)), ("attractor_N_ge_theta", Dv.N_inf >= th),
                       ("max_stable_N_ge_theta", Dv.reg_max_stable_N >= th), ("min_stable_N_ge_theta", Dv.reg_min_stable_N >= th)):
            pv = pv.values
            row[f"acc_{nm}_vs_event240"] = float(np.mean(pv == obs24))
            row[f"acc_{nm}_vs_long_event"] = float(np.mean(pv == obsL))
        # transient crossing consistent with a low attractor
        row["events_with_attractor_below_theta"] = float(((Dv.N_inf < th) & ev).sum() / max(ev.sum(), 1))
        rows.append(row)
        log(f"theta={th:5.2f}: P60={row['P_by_60']:.4f} P120={row['P_by_120']:.4f} P240={row['P_by_240']:.4f}; cond.median={row['conditional_median_T']:.1f}; "
            f"transient={row['frac_events_transient_within_240']:.3f}; P(long)={row['P_long_horizon_event']:.3f}; "
            f"acc kc/event240={row['acc_TF_unstable_vs_event240']:.3f}, attractor/long={row['acc_attractor_N_ge_theta_vs_long_event']:.3f}")
    R = pd.DataFrame(rows)
    save_table(R, "p3_threshold_event_probabilities")

    # ----------------------------------------------------------------- 2. per-sample boundary scan in beta (simulated vs equilibrium)
    sub = D[(D.r_eff > 0)].iloc[:300]
    items = [({k: float(r[k]) for k in ("r", "q", "a", "delta_D", "beta", "rho", "eta", "N0", "fD", "N_detect")}, THS) for _, r in sub.iterrows()]
    sc = pd.DataFrame(parallel_map(scan_beta, items, chunksize=4))
    log(f"beta-boundary scan done for {len(sc)} samples")
    brows = []
    for th in THS:
        k = tkey(th)
        beta_c = sc.r_eff.values
        for nm in ("sim", "eq_min_stable", "eq_max_stable"):
            v = sc[f"beta_det_{nm}_{k}"].values / beta_c
            fin = np.isfinite(v)
            brows.append({"threshold": th, "method": nm, "n_finite": int(fin.sum()), "n_inf_no_boundary": int(np.isinf(v).sum()),
                          "n_nan_no_event_at_min_beta": int(np.isnan(v).sum()), "median_beta_det_over_beta_c": float(np.median(v[fin])) if fin.any() else float("nan"),
                          "q25": float(np.quantile(v[fin], 0.25)) if fin.any() else float("nan"), "q75": float(np.quantile(v[fin], 0.75)) if fin.any() else float("nan"),
                          "min": float(np.min(v[fin])) if fin.any() else float("nan"), "max": float(np.max(v[fin])) if fin.any() else float("nan")})
        ratio = sc[f"beta_det_sim_{k}"].values / sc[f"beta_det_eq_min_stable_{k}"].values
        ratio = ratio[np.isfinite(ratio)]
        brows.append({"threshold": th, "method": "ratio sim/eq_min_stable", "n_finite": len(ratio),
                      "median_beta_det_over_beta_c": float(np.median(ratio)) if len(ratio) else float("nan"),
                      "q25": float(np.quantile(ratio, 0.25)) if len(ratio) else float("nan"), "q75": float(np.quantile(ratio, 0.75)) if len(ratio) else float("nan"),
                      "min": float(ratio.min()) if len(ratio) else float("nan"), "max": float(ratio.max()) if len(ratio) else float("nan")})
        brows.append({"threshold": th, "method": "fraction of samples with re-entrant (non-monotone in beta) events",
                      "median_beta_det_over_beta_c": float(sc[f"reentrant_sim_{k}"].mean()),
                      "n_finite": int(sc[f"reentrant_sim_{k}"].sum())})
        brows.append({"threshold": th, "method": "agreement over beta grid: sim event vs eq (lowest stable N>=theta) / (highest) / k_c-type (TF unstable)",
                      "median_beta_det_over_beta_c": float(sc[f"agree_eq_min_{k}"].mean()), "q25": float(sc[f"agree_eq_max_{k}"].mean()), "q75": float(sc[f"agree_kc_{k}"].mean())})
    B = pd.DataFrame(brows)
    save_table(B, "p3_threshold_beta_boundary")
    save_table(sc, "p3_threshold_beta_boundary_per_sample")

    # ----------------------------------------------------------------- 3. baseline k_det: formula vs simulation per threshold
    from dormancy import Params as P0
    from dormancy.relapse import k_detection_limit
    p0 = P0()
    kc = p0.r_eff * p0.mu / p0.s
    chk = []
    for th in THS:
        kd = k_detection_limit(p0, th)
        res = {"threshold": th, "k_det_formula_over_kc": kd / kc, "beta_det_over_beta_c_formula": kd / kc}
        for f in (0.98, 1.02):
            pk = P0(k=f * kd)
            nz = NormalizedParams.from_params(pk).to_params()
            y0 = [1.01e-4 * 0.01, 1.01e-4 * 0.99, nz.E_free]
            lr = long_run(nz, y0, th, with_event=False)
            res[f"Nmax_long_at_{f}kdet"] = lr.get("N_max_long", float("nan"))
            res[f"event_at_{f}kdet"] = int(lr.get("N_max_long", 0) >= th)
        chk.append(res)
    save_table(pd.DataFrame(chk), "p3_threshold_baseline_kdet_check")
    log("baseline k_det check: " + str([(c["threshold"], c["event_at_0.98kdet"], c["event_at_1.02kdet"]) for c in chk]))

    # ----------------------------------------------------------------- figure
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 9))
    ax = axes[0, 0]
    for h, c in zip(HZ, [BLUE, ORANGE, VIOLET]):
        ax.plot(R.threshold, R[f"P_by_{h}"], "o-", color=c, label=f"by {h} months ({h//12} y)")
        ax.fill_between(R.threshold, R[f"P_by_{h}_wilson_lo"], R[f"P_by_{h}_wilson_hi"], color=c, alpha=0.18)
    ax.plot(R.threshold, R.P_long_horizon_event, "s--", color=INK, ms=4, label="long horizon (>=1e5 mo), reaches threshold")
    ax.axvline(0.1, color=GREY, ls=":"); ax.set_xscale("log")
    ax.set_xlabel("model-defined detection threshold $N_{detect}$ [fraction of K]"); ax.set_ylabel("simulated P(first crossing)")
    ax.set_title("Event probability vs threshold (prior P1, n=5000, 95% Wilson)", fontsize=9); ax.legend(fontsize=7)
    ax = axes[0, 1]
    ax.errorbar(R.threshold, R.conditional_median_T, yerr=[R.conditional_median_T - R.conditional_q25_T, R.conditional_q75_T - R.conditional_median_T],
                fmt="o-", color=BLUE, capsize=3, label="conditional median (IQR) among events")
    ax.plot(R.threshold, R.RMST_240, "s-", color=ORANGE, label="restricted mean event-free time (240 mo)")
    ax.set_xscale("log"); ax.set_xlabel("$N_{detect}$"); ax.set_ylabel("months"); ax.legend(fontsize=7); ax.set_ylim(0, 245)
    ax.set_title("Event-time summaries vs threshold", fontsize=9)
    ax = axes[1, 0]
    ax.plot(R.threshold, R.frac_events_transient_within_240, "o-", color=MAGENTA, label="events later followed by a downward crossing (transient)")
    ax.plot(R.threshold, R.frac_events_persistent_at_240, "o-", color=AQUA, label="events still above threshold at 240 mo")
    ax.plot(R.threshold, R.events_with_attractor_below_theta, "s--", color=INK, label="events whose long-run attractor is below threshold")
    ax.plot(R.threshold, R.acc_TF_unstable_vs_event240, "^:", color=GREY, label="accuracy: 'TF unstable' vs event by 240")
    ax.plot(R.threshold, R.acc_attractor_N_ge_theta_vs_long_event, "v:", color=BLUE, label="accuracy: 'attractor $N\\geq\\theta$' vs long-horizon event")
    ax.set_xscale("log"); ax.set_xlabel("$N_{detect}$"); ax.set_ylabel("fraction"); ax.set_ylim(-0.02, 1.05); ax.legend(fontsize=6.3, loc="center left")
    ax.set_title("Transient vs persistent exceedance; equilibrium predictors", fontsize=9)
    ax = axes[1, 1]
    for nm, c, lab in [("sim", INK, "simulated (long run)"), ("eq_min_stable", BLUE, "equilibrium: lowest stable $N^*\\geq\\theta$"),
                       ("eq_max_stable", ORANGE, "equilibrium: highest stable $N^*\\geq\\theta$")]:
        s = B[B.method == nm]
        ax.plot(s.threshold, s.median_beta_det_over_beta_c, "o-", color=c, label=lab)
        ax.fill_between(s.threshold, s.q25, s.q75, color=c, alpha=0.12)
    kdc = [c["k_det_formula_over_kc"] for c in chk]
    ax.plot(THS, kdc, "k*--", ms=9, mfc="white", label="baseline-parameter formula $k_{det}/k_c$ (Phase 2)")
    ax.axhline(1, color=GREY, ls=":"); ax.set_xscale("log")
    ax.set_xlabel("$N_{detect}$"); ax.set_ylabel("$\\beta_{det}/\\beta_c$ (per-sample boundary; median, IQR)")
    ax.set_ylim(0, 1.15); ax.legend(fontsize=6.3, loc="lower right")
    ax.set_title("Apparent detection boundary: killing above which no relapse occurs", fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 1)); p3_footer(fig, "Thresholds are model-defined burdens, not imaging limits. Same 5000 samples/dense paths at every threshold (common random numbers).")
    fig.savefig(FIG / "fig_p3_09_detection_threshold_sensitivity.png"); plt.close(fig)
    log("done")


if __name__ == "__main__":
    main()
