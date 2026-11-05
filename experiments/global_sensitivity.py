"""Phase 3, Task 2: global sensitivity (LHS + Spearman/PRCC + bootstrap, Morris, Sobol incl. second order).

Run: .venv/bin/python experiments/global_sensitivity.py
Outputs: results/tables/p3_global_*.csv, results/figures/fig_p3_02*.png, fig_p3_03*.png
Sensitivity is relative to the HYPOTHETICAL screening ranges in configs/phase3.yaml; it is not parameter uncertainty.
"""
import math

from _p3common import *  # noqa: F401,F403
import matplotlib.pyplot as plt
import matplotlib.ticker

from dormancy.sensitivity import (bootstrap_ci, kendall_tau, morris_design, morris_effects, prcc, saltelli_unit_matrices,
                                  sobol_bootstrap, sobol_from_outputs, spearman_columns)
from dormancy.uncertainty import build_global_prior, default_options, evaluate_row, parallel_map

NAMES = ["r", "q", "a", "delta_D", "beta", "rho", "eta", "N0", "fD", "N_detect"]
LAB = {"r": "r", "q": "q", "a": "a", "delta_D": "$\\delta_D$", "beta": "$\\beta$", "rho": "$\\rho$", "eta": "$\\eta$",
       "N0": "$N_0$", "fD": "$f_D$", "N_detect": "$N_{det}$"}
OUT_LAB = {"logN240": "$\\log_{10}N(240)$", "Tr": "restricted time $\\min(T,240)$", "event": "event by 240 mo (0/1)"}
opt = default_options(CFG, thresholds=[], check_against_simulate=False, long_run=False, add_main=True)


def run(rows):
    return pd.DataFrame(parallel_map(evaluate_row, [(r, opt) for r in rows]))


def derive(df):
    out = pd.DataFrame({"logN240": np.log10(np.maximum(df.N_240.values, 1e-12)),
                        "event": df.event_main.values,
                        "Tr": np.where(df.event_main.values == 1, df.T_main.values, 240.0),
                        "lambda_TF": df.reg_lambda_tf.values, "maxN_stable": df.reg_max_stable_N.values})
    return out


def main():
    prior = build_global_prior(CFG)
    assert prior.names == NAMES
    seeds = CFG["seeds"]

    # ------------------------------------------------------------------ 1. LHS and convergence in n
    sizes = [250, 500, 1000, 2000, 4000]
    store = {}
    for n in sizes:
        rows = prior.sample(n, seeds["global_lhs"], "lhs")
        df = run(rows)
        store[n] = (pd.DataFrame(rows)[NAMES], df)
        log(f"LHS n={n}: P(event<=240)={np.nanmean(df.event_main):.3f}, failures={(~df.solver_ok).sum()}, "
            f"regimes={df.reg_regime.value_counts().to_dict()}")
    X, D = store[4000]
    Y = derive(D)
    save_table(X.assign(**D[["event_main", "T_main", "N_240", "reg_regime", "beta_over_betac", "r_eff"]].to_dict("list")),
               "p3_global_lhs_samples")
    save_table(D.reg_regime.value_counts().rename("count").reset_index(), "p3_global_design_regimes")

    # rank statistics with bootstrap CIs on the final n
    rows_out = []
    Xv = X.values
    for o in ["logN240", "Tr", "event", "lambda_TF", "maxN_stable"]:
        y = Y[o].values
        rho = spearman_columns(Xv, y)
        pr = prcc(Xv, y)
        lo, hi = bootstrap_ci(spearman_columns, Xv, y, B=400, seed=seeds["bootstrap"])
        plo, phi = bootstrap_ci(prcc, Xv, y, B=400, seed=seeds["bootstrap"])
        for i, nm in enumerate(NAMES):
            rows_out.append({"output": o, "input": nm, "spearman": rho[i], "spearman_lo": lo[i], "spearman_hi": hi[i],
                             "prcc": pr[i], "prcc_lo": plo[i], "prcc_hi": phi[i], "n": len(y)})
    # conditional event time among events only (explicitly conditional)
    ev = (D.event_main == 1).values
    for nm_i, nm in enumerate(NAMES):
        pass
    yT = D.T_main.values[ev]
    rho = spearman_columns(Xv[ev], yT)
    lo, hi = bootstrap_ci(spearman_columns, Xv[ev], yT, B=400, seed=seeds["bootstrap"])
    for i, nm in enumerate(NAMES):
        rows_out.append({"output": "T_given_event (CONDITIONAL, events only)", "input": nm, "spearman": rho[i], "spearman_lo": lo[i],
                         "spearman_hi": hi[i], "prcc": float("nan"), "prcc_lo": float("nan"), "prcc_hi": float("nan"), "n": int(ev.sum())})
    rank_df = pd.DataFrame(rows_out)
    save_table(rank_df, "p3_global_rank_correlations")
    # AUC (Mann-Whitney) of each input for event vs no-event
    from scipy.stats import rankdata
    auc = []
    for i, nm in enumerate(NAMES):
        r = rankdata(Xv[:, i]); n1, n0 = ev.sum(), (~ev).sum()
        auc.append({"input": nm, "AUC_event_vs_none": (r[ev].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)})
    save_table(pd.DataFrame(auc), "p3_global_event_auc")

    # convergence of rankings with n
    conv = []
    ref_rank = np.abs(spearman_columns(store[4000][0].values, derive(store[4000][1])["Tr"].values))
    for n in sizes:
        Xn, Dn = store[n]
        Yn = derive(Dn)
        for o in ["logN240", "Tr", "event"]:
            s = np.abs(spearman_columns(Xn.values, Yn[o].values))
            sref = np.abs(spearman_columns(store[4000][0].values, derive(store[4000][1])[o].values))
            lo, hi = bootstrap_ci(spearman_columns, Xn.values, Yn[o].values, B=200, seed=seeds["bootstrap"])
            conv.append({"n": n, "output": o, "kendall_tau_vs_n4000": kendall_tau(s, sref),
                         "max_abs_diff_spearman": float(np.max(np.abs(s - sref))),
                         "median_ci_halfwidth": float(np.median((hi - lo) / 2)), "P_event": float(np.nanmean(Dn.event_main))})
    save_table(pd.DataFrame(conv), "p3_global_rank_convergence")

    # ------------------------------------------------------------------ 2. Morris screening
    d = len(NAMES)
    Xm, orders, signs = morris_design(d, 40, 4, seeds["morris"])
    delta = 4 / (2 * 3)
    flat = prior.from_unit(Xm.reshape(-1, d))
    Dm = run(flat)
    Ym = derive(Dm)
    mor = []
    for o in ["logN240", "Tr", "event"]:
        res = morris_effects(Ym[o].values.reshape(40, d + 1), orders, signs, delta)
        for i, nm in enumerate(NAMES):
            mor.append({"output": o, "input": nm, "mu": res["mu"][i], "mu_star": res["mu_star"][i], "sigma": res["sigma"][i]})
    mor = pd.DataFrame(mor)
    save_table(mor, "p3_global_morris")
    log("Morris done")

    # ------------------------------------------------------------------ 3. Sobol (Saltelli) with second-order indices
    Nb = 2048
    A, B, AB, BA = saltelli_unit_matrices(d, Nb, seeds["sobol"], second_order=True)
    mats = [A, B] + AB + BA
    cache = TAB / "p3_global_sobol_cache.npz"
    key = f"{CFG_HASH}|{Nb}|{seeds['sobol']}|{d}"
    Ys = None
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        if str(z["key"]) == key:
            Ys = pd.DataFrame({k: z[k] for k in ("logN240", "event", "Tr", "lambda_TF", "maxN_stable")})
            n_fail = int(z["n_fail"])
            log("Sobol design loaded from cache (same config hash, seed, N)")
    if Ys is None:
        allrows = prior.from_unit(np.vstack(mats))
        Ds = run(allrows)
        Ys = derive(Ds)
        n_fail = int((~Ds.solver_ok).sum())
        np.savez(cache, key=key, n_fail=n_fail, **{c: Ys[c].values for c in Ys.columns})
    log(f"Sobol design evaluated: {len(mats) * Nb} runs; failures={n_fail}")
    sob_rows, s2_store, conv_s = [], {}, []
    for o in ["logN240", "Tr", "event"]:
        f = Ys[o].values.reshape(len(mats), Nb)
        fA, fB, fAB, fBA = f[0], f[1], f[2:2 + d], f[2 + d:]
        res = sobol_from_outputs(fA, fB, fAB, fBA)
        ci = sobol_bootstrap(fA, fB, fAB, fBA, B=200, seed=seeds["bootstrap"])
        for i, nm in enumerate(NAMES):
            sob_rows.append({"output": o, "input": nm, "S1": res["S1"][i], "S1_lo": ci["S1"][0][i], "S1_hi": ci["S1"][1][i],
                             "ST": res["ST"][i], "ST_lo": ci["ST"][0][i], "ST_hi": ci["ST"][1][i],
                             "interaction_ST_minus_S1": res["ST"][i] - res["S1"][i], "variance": res["V"], "N_base": Nb})
        s2_store[o] = (res["S2"], ci["S2"])
        S2df = pd.DataFrame(res["S2"], index=NAMES, columns=NAMES)
        save_table(S2df, f"p3_global_sobol_S2_{o}", index=True)
        save_table(pd.DataFrame(ci["S2"][0], index=NAMES, columns=NAMES), f"p3_global_sobol_S2_{o}_ci_lo", index=True)
        save_table(pd.DataFrame(ci["S2"][1], index=NAMES, columns=NAMES), f"p3_global_sobol_S2_{o}_ci_hi", index=True)
        for n in [128, 256, 512, 1024, 2048]:      # nested scrambled-Sobol subsets (powers of two)
            r2 = sobol_from_outputs(fA[:n], fB[:n], fAB[:, :n], None)
            for i, nm in enumerate(NAMES):
                conv_s.append({"output": o, "N_base": n, "input": nm, "S1": r2["S1"][i], "ST": r2["ST"][i]})
    sob = pd.DataFrame(sob_rows)
    save_table(sob, "p3_global_sobol")
    save_table(pd.DataFrame(conv_s), "p3_global_sobol_convergence")
    log("Sobol done")
    print(sob[sob.output == "Tr"].round(3).to_string())
    print(sob[sob.output == "event"].round(3).to_string())

    # ------------------------------------------------------------------ figures
    # Fig 2: global sensitivity ranking (rank correlations + Sobol)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for c, o in enumerate(["logN240", "Tr", "event"]):
        ax = axes[0, c]
        sub = rank_df[rank_df.output == o].set_index("input").reindex(NAMES)
        order = sub.spearman.abs().sort_values().index
        yy = np.arange(len(order))
        sp, pc = sub.loc[order], sub.loc[order]
        ax.barh(yy - 0.18, sp.spearman, height=0.34, color=BLUE, label="Spearman")
        ax.errorbar(sp.spearman, yy - 0.18, xerr=[sp.spearman - sp.spearman_lo, sp.spearman_hi - sp.spearman], fmt="none", ecolor=INK2, lw=0.8)
        ax.barh(yy + 0.18, pc.prcc, height=0.34, color=ORANGE, label="PRCC")
        ax.errorbar(pc.prcc, yy + 0.18, xerr=[pc.prcc - pc.prcc_lo, pc.prcc_hi - pc.prcc], fmt="none", ecolor=INK2, lw=0.8)
        ax.set_yticks(yy); ax.set_yticklabels([LAB[n] for n in order]); ax.axvline(0, color=INK, lw=0.8)
        ax.set_title(f"{OUT_LAB[o]}\nrank correlation (n=4000 LHS, 95% bootstrap CI)", fontsize=9)
        ax.set_xlim(-1, 1)
        if c == 0:
            ax.legend(loc="lower right")
        ax = axes[1, c]
        sb = sob[sob.output == o].set_index("input").reindex(NAMES).loc[order]
        ax.barh(yy - 0.18, sb.S1, height=0.34, color=AQUA, label="first-order $S_1$")
        ax.errorbar(sb.S1, yy - 0.18, xerr=[(sb.S1 - sb.S1_lo).clip(lower=0), (sb.S1_hi - sb.S1).clip(lower=0)], fmt="none", ecolor=INK2, lw=0.8)
        ax.barh(yy + 0.18, sb.ST, height=0.34, color=VIOLET, label="total $S_T$")
        ax.errorbar(sb.ST, yy + 0.18, xerr=[(sb.ST - sb.ST_lo).clip(lower=0), (sb.ST_hi - sb.ST).clip(lower=0)], fmt="none", ecolor=INK2, lw=0.8)
        ax.set_yticks(yy); ax.set_yticklabels([LAB[n] for n in order]); ax.axvline(0, color=INK, lw=0.8)
        ax.set_title(f"Sobol indices (N={Nb}, {len(mats)*Nb} runs)", fontsize=9)
        ax.set_xlim(min(-0.05, sb.S1_lo.min()), 1.0)
        if c == 0:
            ax.legend(loc="lower right")
    fig.suptitle("Global sensitivity of simulated relapse outputs to ten hypothetical inputs (log/logit-uniform ranges)", y=0.998, fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 0.965)); p3_footer(fig, "Sensitivity over hypothetical ranges, not parameter uncertainty. PRCC assumes monotone responses.")
    fig.savefig(FIG / "fig_p3_02_global_sensitivity_ranking.png"); plt.close(fig)

    # Fig 3: interaction heatmaps
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    for ax, o in zip(axes, ["logN240", "Tr", "event"]):
        S2 = s2_store[o][0].copy()
        lo2, hi2 = s2_store[o][1]
        S2[(lo2 <= 0) & (hi2 >= 0)] = np.nan          # keep only second-order indices whose 95% bootstrap CI excludes zero
        np.fill_diagonal(S2, np.nan)
        # draw ST - S1 on the diagonal for orientation
        sb = sob[sob.output == o].set_index("input").reindex(NAMES)
        for i, nm in enumerate(NAMES):
            S2[i, i] = sb.loc[nm, "interaction_ST_minus_S1"]
        vmax = max(0.05, np.nanmax(np.abs(S2)))
        im = ax.imshow(S2, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(d)); ax.set_xticklabels([LAB[n] for n in NAMES], rotation=60, fontsize=8)
        ax.set_yticks(range(d)); ax.set_yticklabels([LAB[n] for n in NAMES], fontsize=8)
        for i in range(d):
            for j in range(d):
                if np.isfinite(S2[i, j]) and abs(S2[i, j]) >= 0.02:
                    ax.text(j, i, f"{S2[i, j]:.2f}", ha="center", va="center", fontsize=6.5)
        ax.set_title(f"{OUT_LAB[o]}\noff-diagonal: $S_{{ij}}$ with CI excluding 0 (blank otherwise); diagonal: $S_T-S_1$", fontsize=8.2)
        fig.colorbar(im, ax=ax, fraction=0.045); ax.grid(False)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p3_footer(fig, "Saltelli (2002) second-order estimator, N=2048 base rows; cells whose 95% bootstrap interval contains 0 are blank (estimator noise).")
    fig.savefig(FIG / "fig_p3_03_parameter_interaction_heatmap.png"); plt.close(fig)

    # Fig 2b: binned scatter for top inputs (checks monotonicity before trusting PRCC)
    top = list(rank_df[rank_df.output == "Tr"].set_index("input").spearman.abs().sort_values(ascending=False).index[:4])
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.2))
    for r_, (o, ylab) in enumerate([("logN240", OUT_LAB["logN240"]), ("Tr", OUT_LAB["Tr"])]):
        for c, nm in enumerate(top):
            ax = axes[r_, c]
            x = X[nm].values; y = Y[o].values
            ax.scatter(x, y, s=2, alpha=0.25, color=GREY, rasterized=True)
            q = np.quantile(x, np.linspace(0, 1, 21))
            mid, med, q10, q90 = [], [], [], []
            for a_, b_ in zip(q[:-1], q[1:]):
                m = (x >= a_) & (x <= b_)
                mid.append(np.median(x[m])); med.append(np.median(y[m])); q10.append(np.quantile(y[m], 0.1)); q90.append(np.quantile(y[m], 0.9))
            ax.fill_between(mid, q10, q90, color=BLUE, alpha=0.2); ax.plot(mid, med, color=BLUE)
            ax.set_xscale("log" if nm != "fD" else "linear"); ax.set_xlabel(LAB[nm]); ax.set_ylabel(ylab, fontsize=8)
            ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    fig.suptitle("Scatter of outputs against the four most influential inputs (binned median and 10-90% band)", y=0.998, fontsize=10)
    fig.tight_layout(rect=(0, 0.05, 1, 0.965)); p3_footer(fig, "n=4000 LHS; saturation at 240 months (flat Tr) is the right-censoring at the horizon.")
    fig.savefig(FIG / "fig_p3_02b_global_scatter_top_inputs.png"); plt.close(fig)

    # Fig 2c: convergence
    cv = pd.DataFrame(conv_s)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    for ax, o in zip(axes[:2], ["Tr", "event"]):
        for i, nm in enumerate(NAMES):
            sub = cv[(cv.output == o) & (cv.input == nm)]
            ax.plot(sub.N_base, sub.ST, marker="o", ms=3, label=LAB[nm])
        ax.set_xscale("log", base=2); ax.set_xlabel("Sobol base sample size N"); ax.set_ylabel("total-order index $S_T$")
        ax.set_title(f"Sobol convergence: {OUT_LAB[o]}", fontsize=9)
    axes[0].legend(fontsize=6, ncol=2)
    cc = pd.DataFrame(conv)
    for o, col in zip(["logN240", "Tr", "event"], [BLUE, ORANGE, VIOLET]):
        s_ = cc[cc.output == o]
        axes[2].plot(s_.n, s_.kendall_tau_vs_n4000, marker="o", color=col, label=OUT_LAB[o])
    axes[2].set_xscale("log"); axes[2].set_ylim(0.2, 1.02); axes[2].set_xlabel("LHS sample size n")
    axes[2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    axes[2].set_ylabel("Kendall $\\tau$ of $|$Spearman$|$ ranking vs n=4000"); axes[2].legend(fontsize=7)
    axes[2].set_title("Stability of the importance ranking", fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p3_footer(fig, "Nested scrambled-Sobol subsets (N=128..2048); independent LHS at each n.")
    fig.savefig(FIG / "fig_p3_02c_global_sensitivity_convergence.png"); plt.close(fig)

    # Fig 2d: Morris
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, o in zip(axes, ["logN240", "Tr", "event"]):
        s_ = mor[mor.output == o]
        ax.scatter(s_.mu_star, s_.sigma, color=BLUE)
        for _, r_ in s_.iterrows():
            ax.annotate(LAB[r_.input], (r_.mu_star, r_.sigma), fontsize=8, xytext=(3, 3), textcoords="offset points")
        mx = max(s_.mu_star.max(), s_.sigma.max())
        ax.plot([0, mx], [0, mx], color=GREY, ls=":", lw=1)
        ax.set_xlabel("$\\mu^*$ (mean |elementary effect|)"); ax.set_ylabel("$\\sigma$ (interaction / nonlinearity)")
        ax.set_title(f"Morris: {OUT_LAB[o]}", fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p3_footer(fig, "40 trajectories, 4 levels (440 runs). Points above the dotted line: effect dominated by interactions/nonlinearity.")
    fig.savefig(FIG / "fig_p3_02d_morris_screening.png"); plt.close(fig)
    log("done")


if __name__ == "__main__":
    main()
