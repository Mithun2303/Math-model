"""Phase 4, Task 2: synthetic parameter recovery (one-factor-at-a-time design around a base condition).

Run: .venv/bin/python experiments/synthetic_recovery.py      (about 5-10 minutes on 8 cores)
Outputs: results/tables/p4_recovery_*.csv; figures fig_p4_01, 05, 06, 07.
SYNTHETIC data; observation model is an assumption (see src/dormancy/observation.py and configs/phase4.yaml).
"""
from _p4common import *  # noqa: F401,F403
import matplotlib.pyplot as plt

from dormancy.observation import ObservationAssumptions, generate_dataset
from dormancy.synthetic_validation import (FREE_SETS, FitProblem, MODEL_FREE, make_bounds, multistart, predicted_curve,
                                           run_recovery_replicate, wald_and_spectrum)
from dormancy.uncertainty import parallel_map

BASE = dict(observed="N", noise_sigma_ln=0.10, lod=1e-6, schedule_months=6.0, follow_up_months=240.0, signal_scale_c=1.0)
NREP = CFG["fitting"]["n_replicates"] + 2
NST = CFG["fitting"]["n_starts"]
LF = CFG["fitting"]["bounds_log_factor"]


def conditions():
    cs = []
    add = lambda name, gt, fset="ic_unknown", **kw: cs.append((name, gt, fset, {**BASE, **kw}))
    add("base: G1, N, sigma 0.10, 6 mo, 240 mo", "G1_baseline")
    for s in (0.001, 0.01, 0.05, 0.25, 0.5):
        add(f"noise sigma {s}", "G1_baseline", noise_sigma_ln=s)
    for st in (1.0, 3.0, 12.0, 24.0):
        add(f"schedule {st:g} mo", "G1_baseline", schedule_months=st)
    for fu in (60.0, 120.0):
        add(f"follow-up {fu:g} mo", "G1_baseline", follow_up_months=fu)
    for ob in ("PD", "Ne", "PDe"):
        add(f"observed {ob}", "G1_baseline", observed=ob)
    for g in ("G2_fast", "G3_slow", "G4_controlled", "G5_lowdormant", "G6_stimulated"):
        add(f"truth {g}", g)
    add("initial state known", "G1_baseline", "ic_known")
    add("signal scale c unknown (true c = 2)", "G1_baseline", "ic_scale_unknown", signal_scale_c=2.0)
    add("high LOD 3e-4 (early data censored)", "G1_baseline", lod=3e-4)
    return cs


def main():
    seeds = CFG["seeds"]
    cs = conditions()
    items, tags = [], []
    for ci, (name, gt, fset, adict) in enumerate(cs):
        th = truth(gt)
        for rep in range(NREP):
            items.append((th, adict, fset, seeds["synthetic_data"] + 1000 * ci + rep, NST, LF, name))
            tags.append((name, gt, fset, adict["observed"], rep))
    import sys
    figs_only = "--figures-only" in sys.argv
    if figs_only:
        df = pd.read_csv(TAB / "p4_recovery_replicates.csv", comment="#")
    else:
        log(f"{len(cs)} conditions x {NREP} replicates = {len(items)} fits")
        out = parallel_map(run_recovery_replicate, items, chunksize=1)
        df = pd.DataFrame(out)
        df["truth"] = [t[1] for t in tags]; df["free_set"] = [t[2] for t in tags]; df["observed"] = [t[3] for t in tags]; df["rep"] = [t[4] for t in tags]
        for k, v in BASE.items():
            pass
        save_table(df.drop(columns=["sv_at_truth", "z_combos_s_gt_3"]), "p4_recovery_replicates")
        save_table(df[["condition", "seed", "sv_at_truth", "z_combos_s_gt_3"]], "p4_recovery_spectra")

        # summary per condition
        rows = []
        for name, g in df.groupby("condition", sort=False):
            fs = FREE_SETS[g.free_set.iloc[0]]
            r = {"condition": name, "truth": g.truth.iloc[0], "observed": g.observed.iloc[0], "n_rep": len(g),
                 "median_n_dirs_s_gt_3": float(g.n_dirs_s_gt_3_at_truth.median()), "frac_multimodal_starts(>1 distinct within chi2 3.84)": float((g["n_starts_within_chi2_3.84"] > 1).mean()),
                 "frac_bound_hit_any": float(g.any_bound_hit.mean()), "median_T_rel_error": float(g.T_rel_error.median()),
                 "median_rmse_lnN_in_window": float(g.rmse_lnN_in_window.median()), "median_abs_lnN_err_240": float(g.abs_lnN_error_at_240.median()),
                 "frac_T_pred_within_10pct": float((g.T_rel_error < 0.1).mean())}
            for nm in fs:
                e = g[f"err_{nm}"].abs()
                r[f"med_abs_err_{nm}"] = float(e.median())
                r[f"frac_within_x1.5_{nm}"] = float((e <= np.log(1.5)).mean())
                r[f"wald_coverage_{nm}"] = float(g[f"covered95_{nm}"].mean())
            r["median_abs_err_all_params"] = float(np.median([g[f"err_{nm}"].abs().median() for nm in fs]))
            rows.append(r)
        summ = pd.DataFrame(rows)
        save_table(summ, "p4_recovery_summary")
        pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
        print(summ[["condition", "median_n_dirs_s_gt_3", "frac_bound_hit_any", "median_abs_err_all_params", "median_T_rel_error", "frac_T_pred_within_10pct", "median_rmse_lnN_in_window"]].round(3).to_string())

    # ------------------------------------------------------------------ figures
    colors = {"N": BLUE, "PD": AQUA, "Ne": ORANGE, "PDe": VIOLET}
    # Fig 1: true vs estimated (sigma 0.10 conditions across truths and observed sets)
    sel = df[(df.condition.str.startswith("truth ")) | (df.condition.str.startswith("observed ")) | (df.condition.str.startswith("base"))].copy()
    params = MODEL_FREE + ["N0", "fD"]
    fig, axes = plt.subplots(3, 3, figsize=(12.5, 11))
    th_by_gt = {gt: truth(gt) for gt in CFG["ground_truths"]}
    for ax, nm in zip(axes.flat, params):
        for obs, col in colors.items():
            g = sel[sel.observed == obs]
            if g.empty:
                continue
            tv = np.array([th_by_gt[t][nm] for t in g.truth])
            ev = np.exp(np.log(tv) + g[f"err_{nm}"].values) if nm != "fD" else 1 / (1 + np.exp(-(np.log(tv / (1 - tv)) + g[f"err_{nm}"].values)))
            ax.scatter(tv, ev, s=14, alpha=0.55, color=col, label=f"observed: {obs}")
        lo, hi = ax.get_xlim()
        xs = np.geomspace(max(1e-6, min(ax.get_xlim()[0], ax.get_ylim()[0])), 2, 50)
        ax.plot(xs, xs, color=INK, lw=0.8)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(f"true {nm}"); ax.set_ylabel(f"estimated {nm}")
        ax.set_title(nm if nm not in ("delta_D",) else "$\\delta_D$", fontsize=9)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Parameter recovery: true vs estimated (6 ground truths x observed-state sets; each point a noisy synthetic data set)", y=0.995, fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97)); p4_footer(fig, "sigma_ln = 0.10, 6-month schedule, truncated at detection; multistart 8 random starts, box bounds truth x/30...x30. Line = identity.")
    fig.savefig(FIG / "fig_p4_01_parameter_recovery_true_vs_estimated.png"); plt.close(fig)

    # Fig 5/6: error vs noise, vs schedule, vs follow-up
    def err_panel(ax, conds, xs, xlabel, logx=True):
        for nm in params:
            ys = []
            for c in conds:
                g = df[df.condition == c]
                ys.append(g[f"err_{nm}"].abs().median() if f"err_{nm}" in g else np.nan)
            ax.plot(xs, ys, "o-", ms=3, label=nm if nm != "delta_D" else "$\\delta_D$")
        ax.axhline(np.log(1.5), color=INK2, ls=":", lw=1)
        ax.text(xs[0], np.log(1.5) * 1.05, "x1.5", fontsize=7, color=INK2)
        if logx:
            ax.set_xscale("log")
        ax.set_yscale("log"); ax.set_xlabel(xlabel); ax.set_ylabel("median |log recovery error|")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    sig = [0.001, 0.01, 0.05, 0.10, 0.25, 0.5]
    cond_s = [f"noise sigma {s}" if s != 0.10 else cs[0][0] for s in sig]
    err_panel(axes[0], cond_s, sig, "measurement noise sd (natural log units)")
    axes[0].set_title("Recovery error vs noise (observed N, G1)", fontsize=9); axes[0].legend(fontsize=6.5, ncol=3)
    ax = axes[1]
    ax.plot(sig, [df[df.condition == c].T_rel_error.median() for c in cond_s], "s-", color=INK, label="crossing-time relative error")
    ax.plot(sig, [df[df.condition == c].median_rmse if False else df[df.condition == c].abs_lnN_error_at_240.median() for c in cond_s], "^-", color=ORANGE, label="|ln N(240) error|")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("measurement noise sd (natural log units)"); ax.set_ylabel("median prediction error")
    ax.set_title("Prediction error vs noise", fontsize=9); ax.legend(fontsize=7)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p4_footer(fig, "G1 baseline, observed total burden N, 6-month visits, truncated at detection; errors that persist at sigma=0.001 reflect weakly determined parameter directions (Fisher spectrum, Fig 4), not noise.")
    fig.savefig(FIG / "fig_p4_05_recovery_error_vs_noise.png"); plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    sch = [1.0, 3.0, 6.0, 12.0, 24.0]
    cond_sch = [f"schedule {s:g} mo" if s != 6.0 else cs[0][0] for s in sch]
    err_panel(axes[0], cond_sch, sch, "interval between visits [months]")
    axes[0].set_title("Recovery error vs observation frequency", fontsize=9); axes[0].legend(fontsize=6, ncol=3)
    fu = [60.0, 120.0, 240.0]
    cond_fu = [f"follow-up {f:g} mo" if f != 240.0 else cs[0][0] for f in fu]
    err_panel(axes[1], cond_fu, fu, "follow-up limit [months] (observation also stops at detection)", logx=False)
    axes[1].set_title("Recovery error vs follow-up", fontsize=9)
    ax = axes[2]
    ax.plot(sch, [df[df.condition == c].median_nobs if False else df[df.condition == c].n_obs_times.median() for c in cond_sch], "o-", color=INK, label="observations per data set")
    ax2 = ax.twinx()
    ax2.plot(sch, [df[df.condition == c].T_rel_error.median() for c in cond_sch], "s-", color=ORANGE, label="crossing-time rel. error")
    ax.set_xscale("log"); ax.set_xlabel("interval between visits [months]"); ax.set_ylabel("observations"); ax2.set_ylabel("median T relative error")
    ax2.grid(False); ax.set_title("Information vs prediction", fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p4_footer(fig, "G1 baseline, observed N, sigma_ln=0.10.")
    fig.savefig(FIG / "fig_p4_06_recovery_error_vs_observation_frequency.png"); plt.close(fig)

    # Fig 7: observed vs reconstructed trajectories (G1 base and G6), with near-best multistart solutions
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    for ax, (gt, obs, ttl) in zip(axes, [("G1_baseline", "N", "G1, observed N"), ("G1_baseline", "PDe", "G1, observed P, D, e"), ("G6_stimulated", "N", "G6 (nonlinear, transient; never reaches 0.1: full 240-month window)")]):
        th = truth(gt)
        adict = {**BASE, "observed": obs, "detection_threshold": 0.1}
        a = ObservationAssumptions(**adict)
        d = generate_dataset(th, a, seed=seeds["synthetic_data"] + 7)
        free = FREE_SETS["ic_unknown"]
        fixed = {"e0": th["e0"], "c": 1.0}
        lo, hi = make_bounds({**th, "c": 1.0}, free, LF)
        prob = FitProblem(d, free, fixed, lo, hi)
        res = multistart(prob, 8, seeds["multistart"])
        tt = np.linspace(0, 240, 241)
        from dormancy.observation import simulate_hidden
        hid = simulate_hidden(th, tt)
        for r_ in res:
            if r_.cost <= res[0].cost + 3.84:
                c = predicted_curve({**th, **prob.to_theta(r_.x)}, tt)
                ax.semilogy(tt, c[:, 0] + c[:, 1], color=ORANGE, lw=0.8, alpha=0.5)
        ax.semilogy(tt, hid.N, color=INK, lw=2.0, label="true N(t) (hidden)")
        for ch, mk in zip(CHANNELS_KEYS[obs], ["o", "s", "^"]):
            yv = d.y[ch]
            ax.semilogy(d.times, yv, mk, ms=3.5, alpha=0.8, label=f"data: {ch}")
        ax.axvline(d.times[-1], color=GREY, ls=":", lw=1); ax.text(d.times[-1], ax.get_ylim()[0] * 3, " observation stops\n (detection/treatment)", fontsize=7, color=INK2)
        ax.axhline(a.detection_threshold, color=MAGENTA, ls="--", lw=0.9)
        ax.plot([], [], color=ORANGE, label="near-best fits (chi2 <= best+3.84), extrapolated")
        ax.set_xlabel("time since start of simulation [months]"); ax.set_ylabel("burden N (fraction of K) / signal")
        ax.set_title(ttl, fontsize=9); ax.legend(fontsize=6.3, loc="lower right"); ax.set_ylim(1e-6, 3)
    fig.suptitle("Synthetic observed vs reconstructed trajectories: good in-window fits do not fix extrapolation", y=0.995, fontsize=10)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96)); p4_footer(fig, "sigma_ln = 0.10; dashed line: detection threshold; fits from 8 random starts.")
    fig.savefig(FIG / "fig_p4_07_observed_vs_reconstructed_trajectories.png"); plt.close(fig)
    log("done")


CHANNELS_KEYS = {"N": ["N"], "PD": ["P", "D"], "Ne": ["N", "e"], "PDe": ["P", "D", "e"]}

if __name__ == "__main__":
    main()
