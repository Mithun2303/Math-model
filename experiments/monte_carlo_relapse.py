"""Phase 3, Tasks 4-6 and 8: Monte Carlo relapse simulations under explicit (hypothetical) uncertainty distributions.

Run: .venv/bin/python experiments/monte_carlo_relapse.py
Outputs: results/tables/p3_mc_*.csv, results/figures/fig_p3_06..12*.png
Simulated event probabilities are model-defined quantities under hypothetical priors, NOT patient risk.
"""
import numpy as np

from _p3common import *  # noqa: F401,F403
import matplotlib.pyplot as plt

from dormancy.sensitivity import saltelli_unit_matrices, sobol_bootstrap, sobol_from_outputs
from dormancy.uncertainty import (build_prior, clopper_pearson, default_options, event_summary, evaluate_row, parallel_map,
                                  upper_bound_rare, wilson_interval)

HZ = CFG["report_horizons_months"]
THS = CFG["thresholds"]
PRIORS = list(CFG["priors"].keys())
PCOL = {"P1_base": INK, "P2_wide": BLUE, "P3_loguniform": ORANGE, "P4_corr_pos": AQUA, "P5_corr_neg": VIOLET, "P6_r_varied": MAGENTA,
        "P7_fixed_initial": GREY}
REG_ORDER = ["TF only", "TF + NZ", "NZ only", "NZ + NZ", "TF + NZ + NZ", "no stable equilibrium"]


def evaluate(rows, **over):
    opt = default_options(CFG, **over)
    return pd.DataFrame(parallel_map(evaluate_row, [(r, opt) for r in rows]))


def summarise(df, name, theta=0.1):
    es = event_summary(df.T_main.values, df.event_main.values, HZ, CFG["horizon_months"])
    es.update(prior=name, n_failures=int((~df.solver_ok).sum()), n_start_above=int(df.start_above_main.sum()))
    es["frac_TF_unstable"] = float(np.mean(~df.reg_tf_stable)) if "reg_tf_stable" in df else float("nan")
    es["frac_bistable_any"] = float(np.mean(df.reg_regime.isin(["TF + NZ", "NZ + NZ", "TF + NZ + NZ"])))
    es["frac_events_transient"] = float(np.nanmean(df.loc[df.event_main == 1, "transient_main"])) if (df.event_main == 1).any() else float("nan")
    es["frac_events_persistent_at_240"] = float(np.nanmean(df.loc[df.event_main == 1, "persistent_main"])) if (df.event_main == 1).any() else float("nan")
    return es


def add_censoring(df, horizon=240.0):
    """Explicit per-simulation event record: indicator, observed time if event, censoring time otherwise."""
    df = df.copy()
    ev = df.event_main == 1
    df["event_indicator"] = df.event_main
    df["time_or_censor"] = np.where(ev, df.T_main, horizon)
    df["censored_at_horizon"] = (df.event_main == 0)
    df["horizon_months"] = horizon
    return df


def main():
    seeds = CFG["seeds"]
    # ------------------------------------------------------------------ 1. base prior: 5000 samples (first 1000 = initial run)
    prior = build_prior(CFG, "P1_base")
    rows = prior.sample(5000, seeds["mc_base"], "mc")
    D = evaluate(rows, long_run=True, check_against_simulate=True)
    log(f"P1: {len(D)} samples; solver failures={int((~D.solver_ok).sum())}; simulate-agreement={D.simulate_agrees.mean():.4f}; "
        f"max|dT|={np.nanmax(D.simulate_dT):.2e}; start_above={int(D.start_above_main.sum())}")
    D = add_censoring(D)
    save_table(D, "p3_mc_P1_samples")
    s1000 = summarise(D.iloc[:1000], "P1_base (n=1000, initial run)")
    s5000 = summarise(D, "P1_base (n=5000)")
    save_table(pd.DataFrame([s1000, s5000]), "p3_mc_P1_summary")
    for s in (s1000, s5000):
        log(f"   {s['prior']}: P60={s['P_by_60']:.4f} [{s['P_by_60_wilson_lo']:.4f},{s['P_by_60_wilson_hi']:.4f}] "
            f"P120={s['P_by_120']:.4f} [{s['P_by_120_wilson_lo']:.4f},{s['P_by_120_wilson_hi']:.4f}] "
            f"P240={s['P_by_240']:.4f} [{s['P_by_240_wilson_lo']:.4f},{s['P_by_240_wilson_hi']:.4f}]; "
            f"cond.median T={s['conditional_median_T']:.1f}; RMST240={s['RMST_240']:.1f}+-{s['RMST_240_se']:.1f}; "
            f"uncond.median={s['unconditional_median_T']}")

    # regime composition and stratified probabilities
    strat = []
    for reg_, g in D.groupby("reg_regime"):
        e = event_summary(g.T_main.values, g.event_main.values, HZ, 240.0)
        strat.append({"regime": reg_, "n": len(g), "fraction": len(g) / len(D), **{k: e[k] for k in
                      ("P_by_60", "P_by_120", "P_by_240", "P_by_240_wilson_lo", "P_by_240_wilson_hi", "conditional_median_T", "RMST_240")},
                      "frac_attractor_above_threshold": float(np.mean(g.N_inf >= 0.1)) if g.N_inf.notna().any() else float("nan"),
                      "event_long_frac": float(np.nanmean(g.event_long))})
    strat = pd.DataFrame(strat)
    save_table(strat, "p3_mc_P1_by_regime")
    by_att = []
    for a_, g in D.groupby("attractor"):
        e = event_summary(g.T_main.values, g.event_main.values, HZ, 240.0)
        by_att.append({"attractor": a_, "n": len(g), "P_by_240": e["P_by_240"], "event_long_frac": float(np.nanmean(g.event_long)),
                       "transient_frac_of_events": float(np.nanmean(g.loc[g.event_main == 1, "transient_main"])) if (g.event_main == 1).any() else float("nan")})
    save_table(pd.DataFrame(by_att), "p3_mc_P1_by_attractor")

    # ------------------------------------------------------------------ 2. equilibrium predictions vs events (Task 6)
    theta = 0.1
    obs = {"event_by_240": D.event_main == 1, "event_long_horizon": D.event_long == 1}
    preds = {"TF unstable (beta < beta_c; k_c-type)": ~D.reg_tf_stable,
             "highest stable equilibrium N* >= theta": D.reg_max_stable_N >= theta,
             "lowest stable equilibrium N* >= theta": D.reg_min_stable_N >= theta,
             "attractor reached has N_inf >= theta": D.N_inf >= theta}
    conf = []
    for pn, pv in preds.items():
        for on, ov in obs.items():
            tp = int((pv & ov).sum()); fp = int((pv & ~ov).sum()); fn = int((~pv & ov).sum()); tn = int((~pv & ~ov).sum())
            conf.append({"predictor": pn, "observed": on, "TP": tp, "FP": fp, "FN": fn, "TN": tn, "accuracy": (tp + tn) / len(D),
                         "sensitivity": tp / max(tp + fn, 1), "specificity": tn / max(tn + fp, 1)})
    conf = pd.DataFrame(conf)
    save_table(conf, "p3_mc_equilibrium_predictions_vs_events")
    print(conf.round(3).to_string())
    # classes of disagreement between the k_c-type classifier / N* classifier and the long-horizon event
    cls = {
        "TF unstable but never reaches theta (long run): attractor below theta (dormancy/low-burden regime)":
            (~D.reg_tf_stable) & (D.event_long == 0),
        "TF stable but event occurs (basin / transient)": D.reg_tf_stable & (D.event_long == 1),
        "attractor N_inf >= theta but no event by 240 (slow: finite-horizon censoring)": (D.N_inf >= theta) & (D.event_main == 0),
        "event by 240 but attractor N_inf < theta (transient exceedance)": (D.event_main == 1) & (D.N_inf < theta),
        "unresolved attractor classification": D.attractor == "unresolved",
        "event_long == 1 but N_max_long < theta (inconsistency check; must be 0)": (D.event_long == 1) & (D.N_max_long < theta - 1e-6)}
    save_table(pd.DataFrame([{"class": k, "count": int(v.sum()), "fraction": float(v.mean())} for k, v in cls.items()]), "p3_mc_discrepancy_classes")
    for k, v in cls.items():
        log(f"   discrepancy: {int(v.sum()):5d}  {k}")

    # ------------------------------------------------------------------ 3. convergence in n and across independent batches
    ns = list(range(100, 5001, 100))
    conv = []
    for n in ns:
        e = event_summary(D.T_main.values[:n], D.event_main.values[:n], HZ, 240.0)
        for h in HZ:
            conv.append({"n": n, "horizon": h, "P": e[f"P_by_{h}"], "lo": e[f"P_by_{h}_wilson_lo"], "hi": e[f"P_by_{h}_wilson_hi"],
                         "binomial_se": float(np.sqrt(e[f"P_by_{h}"] * (1 - e[f"P_by_{h}"]) / n)), "count": e[f"count_by_{h}"],
                         "upper95_onesided": e[f"P_by_{h}_upper95_onesided"]})
        conv[-1]["RMST_240"] = e["RMST_240"]
    conv = pd.DataFrame(conv)
    save_table(conv, "p3_mc_convergence_running")
    batches = []
    batch_frames = []
    for b in range(5):
        rows_b = prior.sample(1000, seeds["mc_base"] + 1000 + b, "mc")
        Db = evaluate(rows_b, long_run=False, check_against_simulate=False)
        batch_frames.append(Db)
        eb = event_summary(Db.T_main.values, Db.event_main.values, HZ, 240.0)
        batches.append({"batch": b, "seed": seeds["mc_base"] + 1000 + b, **{k: eb[k] for k in ("P_by_60", "P_by_120", "P_by_240", "RMST_240", "conditional_median_T", "n_events")}})
    batches = pd.DataFrame(batches)
    save_table(batches, "p3_mc_independent_batches")
    log("batches P_by_240: " + str(batches.P_by_240.round(4).tolist()) + f" (pooled-n=5000 binomial SE={np.sqrt(s5000['P_by_240']*(1-s5000['P_by_240'])/5000):.4f})")

    # ------------------------------------------------------------------ 4. prior/assumption scenarios (model uncertainty), common random numbers
    sc = []
    scen_samples = {}
    for pn in PRIORS:
        pr = build_prior(CFG, pn)
        Dp = evaluate(pr.sample(5000, seeds["mc_assumptions"], "mc"), long_run=False, check_against_simulate=False)
        Dp = add_censoring(Dp)
        scen_samples[pn] = Dp
        save_table(Dp, "p3_mc_P1_base_run_B_samples" if pn == "P1_base" else f"p3_mc_{pn}_samples")
        s = summarise(Dp, pn)
        s["description"] = CFG["priors"][pn].get("description", "")
        sc.append(s)
        log(f"   {pn}: P60={s['P_by_60']:.4f} P120={s['P_by_120']:.4f} P240={s['P_by_240']:.4f}; TF-unstable frac={s['frac_TF_unstable']:.3f}, "
            f"multistable frac={s['frac_bistable_any']:.3f}, failures={s['n_failures']}")
    sc = pd.DataFrame(sc)
    save_table(sc, "p3_mc_prior_scenarios")
    spread = {h: (sc[f"P_by_{h}"].max() - sc[f"P_by_{h}"].min()) for h in HZ}
    mcw = {h: float(np.median(sc[f"P_by_{h}_wilson_hi"] - sc[f"P_by_{h}_wilson_lo"])) for h in HZ}
    save_table(pd.DataFrame([{"horizon": h, "assumption_spread_max_minus_min": spread[h], "median_MC_95CI_width": mcw[h],
                              "ratio_spread_to_MC_width": spread[h] / mcw[h]} for h in HZ]), "p3_mc_assumption_vs_sampling_uncertainty")

    # ------------------------------------------------------------------ 4b. pooled P1 headline from three independent runs; heterogeneity
    runs = {"run A: seed mc_base, n=5000": D, "run B: seed mc_assumptions, n=5000": scen_samples["P1_base"],
            "run C: 5 batches x 1000 (seeds mc_base+1000..1004)": pd.concat(batch_frames, ignore_index=True)}
    pool_rows = []
    for nm, g in list(runs.items()) + [("POOLED (A+B+C), n=15000", pd.concat(list(runs.values()), ignore_index=True))]:
        e = event_summary(g.T_main.values, g.event_main.values, HZ, 240.0)
        row = {"run": nm, "n": e["n"], "min_first_crossing_T": float(np.nanmin(g.T_main.values)), "n_events_by_240": e["n_events"]}
        for h in HZ:
            row.update({f"count_by_{h}": e[f"count_by_{h}"], f"P_by_{h}": e[f"P_by_{h}"],
                        f"P_by_{h}_cp_lo": e[f"P_by_{h}_cp_lo"], f"P_by_{h}_cp_hi": e[f"P_by_{h}_cp_hi"],
                        f"P_by_{h}_wilson_lo": e[f"P_by_{h}_wilson_lo"], f"P_by_{h}_wilson_hi": e[f"P_by_{h}_wilson_hi"],
                        f"P_by_{h}_upper95_onesided": e[f"P_by_{h}_upper95_onesided"], f"RMST_{h}": e[f"RMST_{h}"]})
        row.update(conditional_median_T=e["conditional_median_T"], conditional_q25_T=e["conditional_q25_T"], conditional_q75_T=e["conditional_q75_T"])
        pool_rows.append(row)
    from scipy.stats import chi2_contingency
    for h in HZ:
        tab = np.array([[r_[f"count_by_{h}"], r_["n"] - r_[f"count_by_{h}"]] for r_ in pool_rows[:3]])
        chi2, pv, _, _ = chi2_contingency(tab) if (tab > 0).all() or h > 60 else (np.nan, np.nan, 0, 0)
        pool_rows[-1][f"homogeneity_chi2_p_across_3_runs_{h}"] = float(pv) if np.isfinite(pv) else float("nan")
    save_table(add_censoring(runs["run C: 5 batches x 1000 (seeds mc_base+1000..1004)"]), "p3_mc_P1_batches_run_C_samples")
    pooled = pd.DataFrame(pool_rows)
    save_table(pooled, "p3_mc_P1_pooled")
    pl = pooled.iloc[-1]
    log(f"POOLED P1 (n=15000): P60={pl['P_by_60']:.5f} ({int(pl['count_by_60'])}/15000, CP [{pl['P_by_60_cp_lo']:.5f},{pl['P_by_60_cp_hi']:.5f}]); "
        f"P120={pl['P_by_120']:.4f} [{pl['P_by_120_cp_lo']:.4f},{pl['P_by_120_cp_hi']:.4f}]; P240={pl['P_by_240']:.4f} [{pl['P_by_240_cp_lo']:.4f},{pl['P_by_240_cp_hi']:.4f}]; "
        f"min T={pl['min_first_crossing_T']:.1f}; runs' P60 counts: {[int(r_['count_by_60']) for r_ in pool_rows[:3]]}")

    # persistent-exceedance definition (sensitivity): an event counts only if N never crosses back below the threshold within 240 mo
    prow = []
    for pn, Dp in scen_samples.items():
        ev_first = (Dp.event_main == 1)
        ev_pers = ev_first & (Dp.persistent_main == 1)
        row = {"prior": pn, "n": len(Dp)}
        for h in HZ:
            a = int((ev_first & (Dp.T_main <= h)).sum()); b = int((ev_pers & (Dp.T_main <= h)).sum())
            row.update({f"P_by_{h}_first_crossing": a / len(Dp), f"P_by_{h}_persistent_only": b / len(Dp),
                        f"count_by_{h}_first_crossing": a, f"count_by_{h}_persistent_only": b,
                        f"P_by_{h}_difference": (a - b) / len(Dp)})
        row["fraction_of_first_crossings_transient"] = float(1 - ev_pers.sum() / max(ev_first.sum(), 1))
        prow.append(row)
    save_table(pd.DataFrame(prow), "p3_mc_persistent_definition")
    print(pd.DataFrame(prow)[["prior", "P_by_120_first_crossing", "P_by_120_persistent_only", "P_by_240_first_crossing", "P_by_240_persistent_only"]].round(4).to_string())

    # ------------------------------------------------------------------ 5. solver tolerance, reproducibility, numerical failures
    sub_rows = rows[:600]
    ref = evaluate(sub_rows, long_run=False, check_against_simulate=False)
    tol_rows = []
    for nm, kw in {"rtol1e-6/atol1e-9": dict(rtol=1e-6, atol=1e-9), "rtol1e-9/atol1e-12 (default)": dict(),
                   "rtol1e-12/atol1e-15": dict(rtol=1e-12, atol=1e-15), "Radau rtol1e-9": dict(method="Radau"),
                   "LSODA rtol1e-9": dict(method="LSODA")}.items():
        Dt = evaluate(sub_rows, long_run=False, check_against_simulate=False, **kw)
        both = (ref.event_main == 1) & (Dt.event_main == 1)
        tol_rows.append({"setting": nm, "event_flips_vs_default_tight_ref": int((ref.event_main != Dt.event_main).sum()),
                         "max_abs_dT": float(np.nanmax(np.abs(ref.T_main[both] - Dt.T_main[both]))) if both.any() else float("nan"),
                         "failures": int((~Dt.solver_ok).sum()), "P_by_240": float(np.nanmean(Dt.event_main))})
    save_table(pd.DataFrame(tol_rows), "p3_mc_solver_tolerance")
    print(pd.DataFrame(tol_rows).to_string())
    again = evaluate(prior.sample(5000, seeds["mc_base"], "mc")[:600], long_run=False, check_against_simulate=False)
    serial = pd.DataFrame(parallel_map(evaluate_row, [(r, default_options(CFG, long_run=False, check_against_simulate=False)) for r in rows[:150]], workers=1))
    repro = {"same_seed_rerun_identical_events": bool(np.array_equal(again.event_main.values, ref.event_main.values)),
             "same_seed_rerun_max_abs_dT": float(np.nanmax(np.abs(again.T_main.values - ref.T_main.values))),
             "parallel_vs_serial_identical_events": bool(np.array_equal(serial.event_main.values, ref.event_main.values[:150])),
             "parallel_vs_serial_max_abs_dT": float(np.nanmax(np.abs(serial.T_main.values - ref.T_main.values[:150])))}
    save_table(pd.DataFrame([repro]), "p3_mc_reproducibility")
    log("reproducibility: " + str(repro))

    # ------------------------------------------------------------------ 6. uncertainty propagation: Sobol under the P1 prior
    d = len(prior.names)
    Nb = 1024
    A, B, AB, BA = saltelli_unit_matrices(d, Nb, seeds["sobol"], second_order=False)
    mats = [A, B] + AB
    Ds = evaluate(prior.from_unit(np.vstack(mats)), long_run=False, check_against_simulate=False)
    Tr = np.where(Ds.event_main == 1, Ds.T_main, 240.0).reshape(len(mats), Nb)
    Ev = Ds.event_main.values.reshape(len(mats), Nb)
    logN = np.log10(np.maximum(Ds.N_240.values, 1e-12)).reshape(len(mats), Nb)
    sob = []
    for oname, f in [("restricted time min(T,240)", Tr), ("event by 240", Ev), ("log10 N(240)", logN)]:
        res = sobol_from_outputs(f[0], f[1], f[2:], None)
        ci = sobol_bootstrap(f[0], f[1], f[2:], None, B=200, seed=seeds["bootstrap"])
        for i, nm in enumerate(prior.names):
            sob.append({"prior": "P1_base", "output": oname, "input": nm, "S1": res["S1"][i], "S1_lo": ci["S1"][0][i], "S1_hi": ci["S1"][1][i],
                        "ST": res["ST"][i], "ST_lo": ci["ST"][0][i], "ST_hi": ci["ST"][1][i], "variance": res["V"]})
    sob = pd.DataFrame(sob)
    save_table(sob, "p3_mc_sobol_under_P1")

    # ------------------------------------------------------------------ figures
    T = D.T_main.values; ev = D.event_main.values == 1
    tt = np.linspace(0, 240, 481)
    F = np.array([np.mean(ev & (T <= t)) for t in tt])

    # Fig 6: distribution of simulated relapse times
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    ax = axes[0]
    ax.hist(T[ev], bins=np.arange(0, 241, 8), color=BLUE, edgecolor="white")
    ax.set_xlabel("simulated first-crossing time [months]"); ax.set_ylabel(f"count (of {int(ev.sum())} events among 5000)")
    ax.set_title("Observed crossing times (events by 240 months)", fontsize=9)
    ax.axvline(np.median(T[ev]), color=ORANGE, ls="--"); ax.text(np.median(T[ev]) + 3, ax.get_ylim()[1] * 0.9, f"conditional median\n{np.median(T[ev]):.0f} mo", fontsize=8, color=INK2)
    ax = axes[1]
    ax.hist(T[ev], bins=np.arange(0, 241, 8), color=BLUE, edgecolor="white", weights=np.full(ev.sum(), 1 / 5000))
    ax.bar([252], [1 - ev.mean()], width=10, color=GREY)
    ax.set_xlim(0, 262); ax.set_xticks([0, 60, 120, 180, 252]); ax.set_xticklabels(["0", "60", "120", "180", "no event\n(censored\nat 240)"])
    ax.set_ylabel("probability mass per 8-month bin / censored"); ax.set_title("All 5000 simulations incl. right-censored at 240", fontsize=9)
    ax = axes[2]
    for lab_, gm in D.groupby("reg_regime"):
        e_ = gm.event_main.values == 1
        if len(gm) >= 30 and e_.sum() >= 5:
            ax.hist(gm.T_main.values[e_], bins=np.arange(0, 241, 16), histtype="step", lw=1.6, density=True, label=f"{lab_} (n={len(gm)}, ev={e_.sum()})")
    ax.set_xlabel("first-crossing time [months]"); ax.set_ylabel("density (events only)"); ax.legend(fontsize=7)
    ax.set_title("Crossing-time distribution by equilibrium regime", fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p3_footer(fig, "Prior P1_base, n=5000, N_detect=0.1 (model-defined), E0=tumour-free level.")
    fig.savefig(FIG / "fig_p3_06_relapse_time_distribution.png"); plt.close(fig)

    # Fig 7: cumulative relapse probability
    fig, ax = plt.subplots(figsize=(8.2, 5))
    for pn in PRIORS:
        Dp = scen_samples[pn]
        Fp = np.array([np.mean((Dp.event_main.values == 1) & (Dp.T_main.values <= t)) for t in tt])
        ax.plot(tt, Fp, color=PCOL[pn], lw=1.0, alpha=0.8, label=pn)
    lo = np.array([wilson_interval(int(np.sum(ev & (T <= t))), 5000)[0] for t in tt]); hi = np.array([wilson_interval(int(np.sum(ev & (T <= t))), 5000)[1] for t in tt])
    ax.fill_between(tt, lo, hi, color=INK, alpha=0.15); ax.plot(tt, F, color=INK, lw=2.2, label="P1_base (95% Wilson band)")
    for h in HZ:
        ax.axvline(h, color=GREY, ls=":", lw=0.8)
    ax.set_xlabel("time since start of simulation [months]"); ax.set_ylabel("model-defined cumulative relapse probability $\\hat F(t)$")
    ax.set_title("Simulated cumulative probability of first crossing $N_{detect}=0.1$ (hypothetical priors)", fontsize=9.5)
    ax.legend(fontsize=7, loc="upper left"); ax.set_xlim(0, 240)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p3_footer(fig, "Dotted lines: 5, 10, 20 years. Band: Monte Carlo sampling error only (not assumption uncertainty).")
    fig.savefig(FIG / "fig_p3_07_cumulative_relapse_probability.png"); plt.close(fig)

    # Fig 8: 5/10/20-year probabilities with intervals, across priors
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    for ax, h in zip(axes, HZ):
        xs = np.arange(len(PRIORS))
        p_ = sc[f"P_by_{h}"].values
        lo = sc[f"P_by_{h}_wilson_lo"].values; hi = sc[f"P_by_{h}_wilson_hi"].values
        ax.bar(xs, p_, color=[PCOL[k] for k in PRIORS], alpha=0.8)
        ax.errorbar(xs, p_, yerr=[p_ - lo, hi - p_], fmt="none", ecolor=INK, lw=1.2, capsize=3)
        for x_, v in zip(xs, p_):
            ax.text(x_, v + (hi - p_).max() + 0.01, f"{v:.3f}", ha="center", fontsize=7)
        ax.set_xticks(xs); ax.set_xticklabels([k.split("_")[0] for k in PRIORS], fontsize=8)
        ax.set_ylabel("simulated P(first crossing by horizon)"); ax.set_title(f"{h} months ({h//12} years)", fontsize=9)
        ax.set_ylim(0, max(hi.max() * 1.25, 0.05))
    fig.suptitle("Simulated crossing probabilities at 5, 10, 20 years under seven hypothetical uncertainty scenarios (n=5000 each; bars: 95% Wilson MC interval)", y=0.999, fontsize=9.5)
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    p3_footer(fig, "P1 base | P2 wide | P3 log-uniform | P4/P5 copula +/-0.5 | P6 r varied | P7 fixed initial condition. Spread across priors = assumption uncertainty.")
    fig.savefig(FIG / "fig_p3_08_event_probabilities_5_10_20y.png"); plt.close(fig)

    # Fig 10: convergence
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    for ax, h in zip(axes, HZ):
        c = conv[conv.horizon == h]
        ax.fill_between(c.n, c.lo, c.hi, color=BLUE, alpha=0.2); ax.plot(c.n, c.P, color=BLUE, lw=1.8)
        for b_, r_ in batches.iterrows():
            ax.plot([1000], [r_[f"P_by_{h}"]], "o", color=ORANGE, ms=4, label="independent batches (n=1000)" if b_ == 0 else None)
        ax.set_xlabel("Monte Carlo sample size n"); ax.set_ylabel(f"P(crossing by {h} mo)")
        ax.set_title(f"{h} months: running estimate with 95% Wilson band", fontsize=9)
        ax.set_ylim(0, max(c.hi.max() * 1.2, 0.05))
        if h == HZ[0]:
            ax.legend(fontsize=7)
            ax.text(0.97, 0.5, f"count at n=5000: {int(c['count'].iloc[-1])}\n95% one-sided upper bound\n= {c['upper95_onesided'].iloc[-1]:.4f}", transform=ax.transAxes, ha="right", fontsize=7, color=INK2)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p3_footer(fig, "Prior P1_base; prefix sub-samples of one seeded run. Rare events (5 years): report count and one-sided bound.")
    fig.savefig(FIG / "fig_p3_10_monte_carlo_convergence.png"); plt.close(fig)

    # Fig 11: parameter uncertainty vs event-time uncertainty
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    ax = axes[0]
    bins = [0, 0.15, 0.25, 0.35, 0.5, 0.75, 1.0, 1.5, 5]
    D["bb"] = pd.cut(D.beta_over_betac.replace(np.inf, 99), bins)
    data, labs = [], []
    for b_, g in D.groupby("bb", observed=True):
        e_ = g.event_main.values == 1
        if e_.sum() >= 5:
            data.append(g.T_main.values[e_]); labs.append(f"{b_.left:g}-{b_.right:g}\n(n={len(g)}, ev={e_.sum()})")
    ax.boxplot(data, tick_labels=labs, showfliers=False, patch_artist=True, boxprops=dict(facecolor=BLUE, alpha=0.4))
    ax.set_xlabel("$\\beta/\\beta_c$ (killing relative to the tumour-free stability threshold)"); ax.set_ylabel("first-crossing time of events [months]")
    ax.tick_params(axis="x", labelsize=6.5); ax.set_title("Event-time spread within parameter strata", fontsize=9)
    ax = axes[1]
    sb = sob[sob.output == "restricted time min(T,240)"].set_index("input").sort_values("ST")
    yy = np.arange(len(sb))
    ax.barh(yy - 0.18, sb.S1, height=0.34, color=AQUA, label="$S_1$"); ax.barh(yy + 0.18, sb.ST, height=0.34, color=VIOLET, label="$S_T$")
    ax.errorbar(sb.ST, yy + 0.18, xerr=[(sb.ST - sb.ST_lo).clip(lower=0), (sb.ST_hi - sb.ST).clip(lower=0)], fmt="none", ecolor=INK2, lw=0.8)
    ax.set_yticks(yy); ax.set_yticklabels(sb.index); ax.legend(fontsize=7, loc="lower right")
    ax.set_xlabel("share of Var[min(T,240)] under the P1 prior"); ax.set_title("Which parameter uncertainty drives event-time uncertainty\n(Sobol under the P1 prior, N=1024)", fontsize=9)
    ax = axes[2]
    sb2 = sob[sob.output == "event by 240"].set_index("input").reindex(sb.index)
    ax.barh(yy - 0.18, sb2.S1, height=0.34, color=AQUA, label="$S_1$"); ax.barh(yy + 0.18, sb2.ST, height=0.34, color=VIOLET, label="$S_T$")
    ax.errorbar(sb2.ST, yy + 0.18, xerr=[(sb2.ST - sb2.ST_lo).clip(lower=0), (sb2.ST_hi - sb2.ST).clip(lower=0)], fmt="none", ecolor=INK2, lw=0.8)
    ax.set_yticks(yy); ax.set_yticklabels(sb2.index); ax.legend(fontsize=7, loc="lower right")
    ax.set_xlabel("share of Var[event indicator]"); ax.set_title("Same, for the event indicator", fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p3_footer(fig, "Parameter uncertainty (hypothetical lognormal/beta priors) propagated to simulated event time; not a clinical uncertainty statement.")
    fig.savefig(FIG / "fig_p3_11_parameter_vs_event_time_uncertainty.png"); plt.close(fig)

    # Fig 12: event probability stratified by regime / attractor
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    ax = axes[0]
    st_ = strat.set_index("regime").reindex([r for r in REG_ORDER if r in set(strat.regime)])
    xs = np.arange(len(st_))
    ax.bar(xs - 0.2, st_.P_by_240, width=0.38, color=BLUE, label="event by 240 mo")
    ax.bar(xs + 0.2, st_.event_long_frac, width=0.38, color=ORANGE, label="event by long horizon (>=1e5 mo)")
    for x_, (n_, f_) in enumerate(zip(st_.n, st_.fraction)):
        ax.text(x_, 1.06, f"n={n_}\n({100*f_:.1f}%)", ha="center", fontsize=7)
    for x_, (a_, b_) in enumerate(zip(st_.P_by_240, st_.event_long_frac)):
        ax.text(x_ - 0.2, a_ + 0.01, f"{a_:.3f}", ha="center", fontsize=6.5); ax.text(x_ + 0.2, b_ + 0.01, f"{b_:.3f}", ha="center", fontsize=6.5)
    ax.set_xticks(xs); ax.set_xticklabels(st_.index, rotation=15, fontsize=8); ax.set_ylim(0, 1.38)
    ax.set_ylabel("simulated event probability ($N_{detect}$=0.1)"); ax.legend(fontsize=7, loc="upper left", ncol=2)
    ax.set_title("By equilibrium regime of the sampled parameter vector", fontsize=9)
    ax = axes[1]
    ba = pd.DataFrame(by_att).set_index("attractor")
    xs = np.arange(len(ba))
    ax.bar(xs - 0.2, ba.P_by_240, width=0.38, color=BLUE, label="event by 240 mo")
    ax.bar(xs + 0.2, ba.event_long_frac, width=0.38, color=ORANGE, label="event by long horizon")
    for x_, n_ in enumerate(ba.n):
        ax.text(x_, 1.06, f"n={n_}", ha="center", fontsize=7)
    for x_, (a_, b_) in enumerate(zip(ba.P_by_240, ba.event_long_frac)):
        ax.text(x_ - 0.2, a_ + 0.01, f"{a_:.3f}", ha="center", fontsize=6.5); ax.text(x_ + 0.2, b_ + 0.01, f"{b_:.3f}", ha="center", fontsize=6.5)
    ax.set_xticks(xs); ax.set_xticklabels(ba.index); ax.set_ylim(0, 1.38); ax.legend(fontsize=7, loc="upper left", ncol=2)
    ax.set_title("By attractor actually reached from the sampled initial state", fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p3_footer(fig, "Prior P1_base, n=5000. Regime = stable equilibria of the parameter vector; stability does not replace event simulation.")
    fig.savefig(FIG / "fig_p3_12_event_probability_by_regime.png"); plt.close(fig)
    log("done")


if __name__ == "__main__":
    main()
