"""Phase 4, Task 6: time-origin, detection-delay, threshold, residual-burden and death-first NUISANCE sensitivity (SYNTHETIC).

Run: .venv/bin/python -u experiments/time_origin_sensitivity.py
Population: 5000 prior-P1 samples (Phase 3 run A). All mappings are hypothetical nuisance assumptions; no clinical delay is
estimated and no model curve is compared with METABRIC. Outputs: p4_time_origin_*.csv, p4_threshold_*.csv, p4_equivalence_map.csv;
figures fig_p4_09, fig_p4_10 (SYNTHETIC) and fig_p4_12 (CLINICAL, descriptive).
"""
from _p4common import *  # noqa: F401,F403
import matplotlib.pyplot as plt

from dormancy.endpoint_mapping import (ClinicalMapping, cumulative_incidence, dkw_epsilon, map_event_times, sup_distance, time_to_level)
from dormancy.relapse import unstable_mode
from dormancy.uncertainty import load_phase3_config, parallel_map, row_to_model, run_trajectory

CFG3 = load_phase3_config()
TH_DENSE = sorted(set(list(np.round(np.geomspace(0.01, 0.5, 13), 5)) + [0.0333, 0.1, 0.3]))
TH_SCALED = [0.1]
SCALES = [0.1, 0.2, 1 / 3, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0]
LEVELS = [0.05, 0.10, 0.20, 0.30, 0.40]
TGRID = np.arange(0, 301, 1.0)
HZ = [60, 120, 240]
KEYS = ["r", "q", "a", "delta_D", "beta", "rho", "eta", "N0", "fD", "N_detect"]


def tk(th):
    return f"{th:g}".replace(".", "p")


def sim_worker(args):
    """One sample: first-crossing times at dense thresholds (scale 1) and at TH_SCALED for each N0 scale; also lambda_plus."""
    row = args
    out = {"lam": float("nan")}
    p, y0, _ = row_to_model(row, CFG3["fixed"])
    try:
        out["lam"] = float(unstable_mode(p).lam) if unstable_mode(p).valid else float("nan")
    except Exception:
        pass
    for s in SCALES:
        r2 = dict(row); r2["N0"] = row["N0"] * s
        p2, y02, _ = row_to_model(r2, CFG3["fixed"])
        ths = TH_DENSE if s == 1.0 else TH_SCALED
        ths = [t for t in ths if y02[0] + y02[1] < t]
        t = run_trajectory(p2, y02, ths, 240.0, report_times=(240.0,), main_threshold=ths[0] if ths else 0.1) if ths else {}
        for th in (TH_DENSE if s == 1.0 else TH_SCALED):
            if th not in ths:
                out[f"s{s:.4f}_T{tk(th)}"] = float("nan"); out[f"s{s:.4f}_above{tk(th)}"] = 1
            else:
                T = t.get(f"T_{tk(th)}", float("nan"))
                out[f"s{s:.4f}_T{tk(th)}"] = float(T) if np.isfinite(T) else float("inf")
                out[f"s{s:.4f}_above{tk(th)}"] = 0
    return out


def F_of(T, shift=0.0, mapping=None, rng=None):
    T = np.asarray(T, float)
    T = T[~np.isnan(T)]
    if mapping is not None:
        T = map_event_times(T, mapping, rng)
    else:
        T = np.where(np.isfinite(T), T + shift, np.inf)
    return cumulative_incidence(T, TGRID)


def levels_row(F):
    return {f"t_to_{int(l * 100)}pct": time_to_level(F, TGRID, l) for l in LEVELS}


def main():
    base = pd.read_csv(TAB / "p3_mc_P1_samples.csv", comment="#")
    rows = [{k: float(r[k]) for k in KEYS} for _, r in base.iterrows()]
    log(f"re-simulating {len(rows)} samples x {len(SCALES)} N0 scales (workers=3)")
    sims = pd.DataFrame(parallel_map(sim_worker, rows, workers=3, chunksize=25))
    save_table(sims.describe().T.reset_index().rename(columns={"index": "column"}), "p4_time_origin_resim_columns_summary")
    log("re-simulation done")
    n = len(sims)
    T1 = sims[f"s1.0000_T{tk(0.1)}"].values

    # ------------------------------------------------------------ 1. tau x d
    F0 = F_of(T1)
    out = []
    for tau in (0, 3, 6, 12, 18, 24):
        for d in (0, 3, 6, 12):
            F = F_of(T1, tau + d)
            row = {"tau": tau, "detect_delay": d, "total_shift": tau + d, **levels_row(F)}
            for h in HZ:
                row[f"F_clin_{h}mo"] = float(F[int(h)])
                row[f"dF_vs_unmapped_{h}mo"] = float(F[int(h)] - F0[int(h)])
            for l in LEVELS:
                a, b = time_to_level(F0, TGRID, l), row[f"t_to_{int(l * 100)}pct"]
                row[f"shift_to_{int(l * 100)}pct"] = b - a if np.isfinite(a) and np.isfinite(b) else float("nan")
            out.append(row)
    t1 = pd.DataFrame(out)
    save_table(t1, "p4_time_origin_tau_delay")
    print(t1[t1.detect_delay == 0][["tau", "t_to_10pct", "t_to_20pct", "t_to_30pct", "F_clin_120mo", "F_clin_240mo", "dF_vs_unmapped_240mo"]].round(3).to_string())

    # ------------------------------------------------------------ 2. thresholds
    out = []
    Fth = {}
    for th in TH_DENSE:
        col = f"s1.0000_T{tk(th)}"
        valid = ~sims[f"s1.0000_above{tk(th)}"].astype(bool)
        Fth[th] = F_of(sims.loc[valid, col].values)
        row = {"threshold": th, "n_valid": int(valid.sum()), "n_start_above_excluded": int((~valid).sum()), **levels_row(Fth[th])}
        for h in HZ:
            row[f"F_{h}mo"] = float(Fth[th][h])
        out.append(row)
    t2 = pd.DataFrame(out)
    save_table(t2, "p4_threshold_effects")
    print(t2[["threshold", "n_start_above_excluded", "t_to_10pct", "t_to_20pct", "t_to_30pct", "F_120mo", "F_240mo"]].round(3).to_string())

    # ------------------------------------------------------------ 3. N0 scale
    out, Fs = [], {}
    for s in SCALES:
        col = f"s{s:.4f}_T{tk(0.1)}"
        valid = ~sims[f"s{s:.4f}_above{tk(0.1)}"].astype(bool)
        Fs[s] = F_of(sims.loc[valid, col].values)
        both = valid & np.isfinite(sims[col]) & np.isfinite(T1) & np.isfinite(sims["lam"])
        dT = (sims.loc[both, col] - T1[both.values]).values
        pred = (-np.log(s) / sims.loc[both, "lam"]).values
        row = {"N0_scale": s, "n_valid": int(valid.sum()), "n_start_above_excluded": int((~valid).sum()), "n_events_both_scales": int(both.sum()),
               "median_dT_numeric": float(np.median(dT)) if both.sum() else float("nan"),
               "median_dT_linear_prediction (-ln s / lambda_plus)": float(np.median(pred)) if both.sum() else float("nan"),
               "median_abs_diff_numeric_vs_prediction": float(np.median(np.abs(dT - pred))) if both.sum() else float("nan"), **levels_row(Fs[s])}
        for h in HZ:
            row[f"F_{h}mo"] = float(Fs[s][h])
        out.append(row)
    t3 = pd.DataFrame(out)
    save_table(t3, "p4_time_origin_N0_scale")
    print(t3[["N0_scale", "n_events_both_scales", "median_dT_numeric", "median_dT_linear_prediction (-ln s / lambda_plus)", "t_to_10pct", "F_240mo"]].round(2).to_string())

    # ------------------------------------------------------------ 4. death-first
    out, Fd = [], {}
    rng = np.random.default_rng(CFG["seeds"]["endpoint_mapping"])
    for fr in (0.0, 0.105, 0.272):
        m = ClinicalMapping(0.0, 0.0, fr, 24.0)
        Fd[fr] = F_of(T1, mapping=m, rng=np.random.default_rng(CFG["seeds"]["endpoint_mapping"]))
        row = {"death_first_fraction": fr, "death_lag_mean": 24.0, **levels_row(Fd[fr])}
        for h in HZ:
            row[f"F_{h}mo"] = float(Fd[fr][h]); row[f"dF_vs_0_{h}mo"] = float(Fd[fr][h] - Fd[0.0][h])
        out.append(row)
    save_table(pd.DataFrame(out), "p4_time_origin_death_first")
    print(pd.DataFrame(out)[["death_first_fraction", "t_to_10pct", "t_to_20pct", "F_120mo", "F_240mo"]].round(3).to_string())

    # ------------------------------------------------------------ 5. equivalence map
    cand = {}
    for th in TH_DENSE:
        cand[("threshold", th)] = Fth[th]
    for s in SCALES:
        cand[("N0_scale", s)] = Fs[s]
    eps = {"n=2380 (METABRIC usable)": dkw_epsilon(2380), "n=5000": dkw_epsilon(5000)}
    out = []
    W = TGRID <= 240
    for tau in (3, 6, 12, 24):
        ref = F_of(T1, tau)
        for kind in ("threshold", "N0_scale"):
            best = min(((sup_distance(ref[W], F[W]), v) for (k, v), F in cand.items() if k == kind), key=lambda x: x[0])
            out.append({"tau": tau, "equivalent_to": kind, "best_value": best[1], "min_sup_distance": best[0],
                        "dkw_eps_n2380": eps["n=2380 (METABRIC usable)"], "dkw_eps_n5000": eps["n=5000"],
                        "distance_below_dkw_n2380 (event times cannot separate)": bool(best[0] < eps["n=2380 (METABRIC usable)"]),
                        "distance_below_dkw_n5000": bool(best[0] < eps["n=5000"]),
                        "distance_without_any_change (tau vs tau=0, theta=0.1)": sup_distance(ref[W], F0[W])})
    # tau equivalent to a factor-3 change of threshold or N0 (shift tau' applied to the unchanged reference threshold 0.1 / s = 1)
    for lab, F in (("threshold x3 (0.3)", Fth[0.3]), ("threshold /3 (0.0333)", Fth[0.0333]), ("N0 x3", Fs[3.0]), ("N0 /3", Fs[1 / 3])):
        best = min(((sup_distance(F[W], F_of(T1, tp)[W]), tp) for tp in np.arange(-60, 61, 1.0)), key=lambda x: x[0])
        out.append({"tau": float("nan"), "equivalent_to": f"tau equivalent of {lab}", "best_value": best[1], "min_sup_distance": best[0],
                    "dkw_eps_n2380": eps["n=2380 (METABRIC usable)"], "dkw_eps_n5000": eps["n=5000"],
                    "distance_below_dkw_n2380 (event times cannot separate)": bool(best[0] < eps["n=2380 (METABRIC usable)"]),
                    "distance_below_dkw_n5000": bool(best[0] < eps["n=5000"]), "distance_without_any_change (tau vs tau=0, theta=0.1)": sup_distance(F[W], F0[W])})
    eq = pd.DataFrame(out)
    save_table(eq, "p4_equivalence_map")
    print(eq.round(4).to_string())

    # ------------------------------------------------------------ figures
    cm = plt.cm.viridis
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    ax = axes[0]
    for i, th in enumerate(TH_DENSE):
        ax.plot(TGRID[W], Fth[th][W], color=cm(i / (len(TH_DENSE) - 1)), lw=1.3)
    sm = plt.cm.ScalarMappable(cmap=cm, norm=matplotlib.colors.LogNorm(0.01, 0.5)); fig.colorbar(sm, ax=ax, label="detection threshold $N_{detect}$")
    ax.set_xlabel("model time since post-treatment origin [months]"); ax.set_ylabel("model-defined cumulative first-crossing probability")
    ax.set_title("F(t) for thresholds 0.01-0.5 (prior P1, n=5000)", fontsize=9)
    ax = axes[1]
    for l, col in zip((0.10, 0.20, 0.30), (BLUE, ORANGE, VIOLET)):
        ax.plot(t2.threshold, t2[f"t_to_{int(l * 100)}pct"], "o-", color=col, label=f"time to {int(l * 100)}% cumulative probability")
    ax.set_xscale("log"); ax.set_xlabel("$N_{detect}$"); ax.set_ylabel("months"); ax.legend(fontsize=7); ax.set_title("Horizontal shift caused by the threshold", fontsize=9)
    ax = axes[2]
    for h, col in zip(HZ, (BLUE, ORANGE, VIOLET)):
        ax.plot(t2.threshold, t2[f"F_{h}mo"], "s-", color=col, label=f"F at {h} months")
    ax.set_xscale("log"); ax.set_xlabel("$N_{detect}$"); ax.set_ylabel("cumulative probability"); ax.legend(fontsize=7); ax.set_title("Vertical effect of the threshold", fontsize=9)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97)); p4_footer(fig, "Thresholds are model-defined burdens (not imaging limits); samples starting above a threshold are excluded (0 for all thresholds >= 0.01 in this prior).")
    fig.savefig(FIG / "fig_p4_09_detection_threshold_uncertainty.png"); plt.close(fig)

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.6))
    ax = axes[0, 0]
    for i, tau in enumerate((0, 3, 6, 12, 18, 24)):
        ax.plot(TGRID, F_of(T1, tau), color=cm(i / 5), label=f"tau = {tau}")
    ax.plot(TGRID, F_of(T1, 24 + 12), "k--", lw=1, label="tau + d = 36")
    ax.set_xlim(0, 300); ax.set_xlabel("clinical time [months]"); ax.set_ylabel("cumulative probability"); ax.legend(fontsize=7)
    ax.set_title("Time-origin offset tau (and delay d: only tau + d matters)", fontsize=9)
    ax = axes[0, 1]
    tt = t1[t1.detect_delay == 0]
    for l, col in zip((0.05, 0.10, 0.20, 0.30), (AQUA, BLUE, ORANGE, VIOLET)):
        ax.plot(tt.tau, tt[f"shift_to_{int(l * 100)}pct"], "o-", color=col, label=f"{int(l * 100)}% level")
    ax.set_xlabel("tau [months]"); ax.set_ylabel("horizontal shift of the curve [months]"); ax.legend(fontsize=7); ax.set_title("Horizontal shift = tau (+ d) exactly", fontsize=9)
    ax = axes[0, 2]
    for i, s in enumerate((0.1, 1 / 3, 1.0, 3.0, 10.0)):
        ax.plot(TGRID[W], Fs[s][W], color=cm(i / 4), label=f"N0 x {s:.2g}")
    ax.set_xlabel("model time [months]"); ax.set_ylabel("cumulative probability"); ax.legend(fontsize=7); ax.set_title("Residual-burden scale (threshold 0.1)", fontsize=9)
    ax = axes[1, 0]
    for fr, col in zip((0.0, 0.105, 0.272), (BLUE, ORANGE, VIOLET)):
        ax.plot(TGRID, Fd[fr], color=col, label=f"death-first share {fr:g}")
    ax.set_xlim(0, 300); ax.set_xlabel("clinical time [months] (tau = 0)"); ax.set_ylabel("cumulative probability"); ax.legend(fontsize=7)
    ax.set_title("Death-first mapping (lag mean 24 mo; shares = METABRIC bounds)", fontsize=9)
    ax = axes[1, 1]
    for tau in (6, 12, 24):
        ax.plot(TGRID[W], F_of(T1, tau)[W], color=BLUE, lw=0.8, alpha=0.4)
    e12 = eq[(eq.tau == 12)]
    th12 = float(e12[e12.equivalent_to == "threshold"].best_value.iloc[0]); s12 = float(e12[e12.equivalent_to == "N0_scale"].best_value.iloc[0])
    ax.plot(TGRID[W], F_of(T1, 12)[W], color=INK, lw=2, label="reference: tau = 12, threshold 0.1")
    ax.plot(TGRID[W], Fth[th12][W], color=ORANGE, label=f"best threshold alone: {th12:g} (tau 0)")
    ax.plot(TGRID[W], Fs[s12][W], color=VIOLET, ls="--", label=f"best N0 scale alone: x{s12:.2g} (tau 0)")
    ax.set_xlabel("time [months]"); ax.set_ylabel("cumulative probability"); ax.legend(fontsize=7); ax.set_title("Equivalence: a tau = 12 shift reproduced by threshold / N0", fontsize=9)
    ax = axes[1, 2]
    xs = np.arange(4); w = 0.38
    for k, (kind, col) in enumerate((("threshold", ORANGE), ("N0_scale", VIOLET))):
        v = [eq[(eq.tau == t) & (eq.equivalent_to == kind)].min_sup_distance.iloc[0] for t in (3, 6, 12, 24)]
        ax.bar(xs + (k - 0.5) * w, v, w, color=col, label=f"best {kind} alone")
    ax.axhline(dkw_epsilon(2380), color=MAGENTA, ls="--", label="DKW noise, n = 2380"); ax.axhline(dkw_epsilon(5000), color=GREY, ls=":", label="DKW noise, n = 5000")
    ax.set_xticks(xs); ax.set_xticklabels([f"tau = {t}" for t in (3, 6, 12, 24)]); ax.set_ylabel("min sup-norm distance of F over [0, 240]"); ax.legend(fontsize=7)
    ax.set_title("Is the nuisance separable from threshold / N0?", fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97)); p4_footer(fig, "Population: prior-P1 events (Phase 3 run A, n = 5000). All mappings are hypothetical nuisance assumptions; no delay distribution is estimated.")
    fig.savefig(FIG / "fig_p4_10_time_origin_uncertainty.png"); plt.close(fig)

    # fig 12 (CLINICAL, descriptive)
    km = pd.read_csv(TAB / "metabric_km_overall.csv")
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 5), gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    ax.set_xlim(-6, 64); ax.set_ylim(-1.2, 10); ax.axis("off"); ax.grid(False)
    ev = [("initial diagnosis\n(origin of OS_MONTHS: documented)", 0, 8.4, INK), ("surgery\n(start of 'disease-free after surgery', Rueda 2019)", 4, 7.2, INK),
          ("end of treatment\n(not recorded; no dates)", 18, 6.0, INK2), ("model t = 0\n(post-treatment residual disease)", 24, 4.8, BLUE),
          ("model threshold crossing\n(N = N_detect)", 44, 3.6, ORANGE), ("recorded recurrence\n(RFS event; composite)", 52, 2.4, MAGENTA), ("disease-specific death\n(also an RFS event)", 60, 1.2, VIOLET)]
    for lab, x, y, col in ev:
        ax.plot([x], [y], "o", color=col, ms=8); ax.text(x + 1.0, y, lab, fontsize=7.5, va="center", color=col)
    ax.annotate("", xy=(0, 9.3), xytext=(24, 9.3), arrowprops=dict(arrowstyle="<->", color=GREY)); ax.text(12, 9.45, "tau: unknown, no treatment dates", ha="center", fontsize=8, color=INK2)
    ax.annotate("", xy=(44, 4.2), xytext=(52, 4.2), arrowprops=dict(arrowstyle="<->", color=GREY)); ax.text(48, 4.4, "d: unknown", ha="center", fontsize=8, color=INK2)
    ax.annotate("", xy=(52, 0.1), xytext=(60, 0.1), arrowprops=dict(arrowstyle="<->", color=GREY)); ax.text(56, 0.3, "L: unknown", ha="center", fontsize=8, color=INK2)
    ax.text(-5.5, -1.1, "Schematic order only (not to scale). VERIFIED: OS origin = diagnosis; RFS origin NOT documented; no treatment dates.\nUNVERIFIED nuisance windows: tau, d, L. Event-time curves cannot separate tau from threshold / N0 (equivalence map).", fontsize=7.2, va="bottom")
    ax.set_title("Clinical endpoint mapping: what is documented and what is a nuisance", fontsize=9.5)
    ax = axes[1]
    yrs = km.horizon_years.values * 12
    ax.errorbar(yrs, km.cumulative_incidence_1_minus_S, yerr=[km.cumulative_incidence_1_minus_S - (1 - km.ci95_high), (1 - km.ci95_low) - km.cumulative_incidence_1_minus_S],
                fmt="o-", color=BLUE, capsize=3)
    for x, y, nr in zip(yrs, km.cumulative_incidence_1_minus_S, km.n_at_risk):
        ax.text(x, y + 0.06, f"at risk {int(nr)}", ha="center", fontsize=7.5, color=INK2)
    ax.set_xlabel("months from the (undocumented) RFS origin"); ax.set_ylabel("Kaplan-Meier cumulative RFS event probability")
    ax.set_xlim(0, 260); ax.set_ylim(0, 0.75); ax.set_title("METABRIC export: RFS (composite; non-disease deaths censored), n = 2,380, 964 events", fontsize=8.5)
    fig.tight_layout(rect=(0, 0.07, 1, 1)); p4_footer(fig, "Descriptive only: no model fitted or overlaid; 20-year estimate rests on 162 patients at risk.", kind="CLIN")
    fig.savefig(FIG / "fig_p4_12_clinical_endpoint_mapping_limitations.png"); plt.close(fig)
    log("done")


import matplotlib  # noqa: E402

if __name__ == "__main__":
    main()
