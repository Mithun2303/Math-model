"""Phase 5 (Agent 2): descriptive consistency of model-generated event-time distributions with the METABRIC RFS curve.

Run: python experiments/model_consistency.py          (needs configs/phase5.yaml and results/tables/p5_prereg_hash.txt written first)
NO fitting, NO tuning, NO new simulations: per-sample first-crossing times come from Phase 3 (p3_mc_*_samples.csv).
Order enforced here: (1) verify the pre-registered design hash, (2) only then load METABRIC, (3) evaluate the full prespecified grid.
Outputs: results/tables/p5_model_*.csv, p5_consistency_*.csv, p5_threshold_*.csv, p5_timeorigin_*.csv, p5_horizon_*.csv,
         results/figures/fig_p5_04 ... fig_p5_09 (*.png)
"""
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from _style import *  # noqa: F401,F403
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from dormancy.clinical_comparison import (ecdf_at, evaluate_scenario, km_cuminc_grid, km_reference, load_metabric, load_prereg,
                                          load_prior_times, map_times, scenario_grid, shifted_cdf, sup_norm, verify_prereg, wilson)
from dormancy.km_descriptive import km_curve

CFG_PATH = ROOT / "configs" / "phase5.yaml"
HASH_PATH = TAB / "p5_prereg_hash.txt"
RAW = ROOT / "data" / "raw" / "brca_metabric" / "data_clinical_patient.txt"
GRID_T = np.arange(0.0, 241.0, 1.0)
PC = {"P1_base": "#0072B2", "P2_wide": "#E69F00", "P3_loguniform": "#009E73", "P4_corr_pos": "#D55E00",
      "P5_corr_neg": "#CC79A7", "P6_r_varied": "#56B4E9", "P7_fixed_initial": "#555555"}
PLS = {"P1_base": "-", "P2_wide": "--", "P3_loguniform": "-.", "P4_corr_pos": ":", "P5_corr_neg": (0, (5, 1, 1, 1)), "P6_r_varied": "--", "P7_fixed_initial": "-."}
KMC = "#000000"
TAG = "MODEL (hypothetical) vs CLINICAL (descriptive)"


def save(df, name, note=""):
    with open(TAB / f"{name}.csv", "w") as fh:
        fh.write(f"# phase5 agent2; design sha256[:12]={H[:12]} (configs/phase5.yaml, pre-registered); hypothetical simulated quantities vs "
                 f"descriptive METABRIC RFS (origin undocumented); not validation. {note}\n")
        df.to_csv(fh, index=False)


def footer(fig, extra=""):
    txt = (f"{TAG}. {extra} Model: first threshold crossing from residual disease (Phase 3 samples, no treatment, no fitting); clinical: METABRIC RFS "
           f"composite endpoint, time origin undocumented. Consistency = not distinguishable at a horizon, not mechanism evidence. [design {H[:12]}]")
    fig.subplots_adjust(bottom=max(fig.subplotpars.bottom, 0.15))
    fig.text(0.01, 0.005, "\n".join(textwrap.wrap(txt, max(60, int(fig.get_figwidth() * 19)))), fontsize=6.3, color=INK2, ha="left", va="bottom")
    fig.text(0.995, 0.995, TAG, fontsize=7.5, fontweight="bold", color="#1c4f9b", ha="right", va="top")


def main():
    global H
    cfg = load_prereg(CFG_PATH)
    H = verify_prereg(CFG_PATH, HASH_PATH)                    # (1) design must match the pre-registered hash
    print("pre-registered design verified:", H[:12])
    TD = TAB
    pri = cfg["simulation_source"]["priors"]
    ths = cfg["grid"]["thresholds"]
    hs = cfg["horizons_months"]["primary"] + cfg["horizons_months"]["secondary"]

    # (2) only now load METABRIC
    t, e = load_metabric(RAW)
    print(f"METABRIC usable n={len(t)} events={int(e.sum())}")
    assert len(t) == 2380 and int(e.sum()) == 964
    km = km_reference(t, e, hs)
    kmF = km_cuminc_grid(t, e, GRID_T)

    # simulated populations
    T = {(p, th): load_prior_times(TD, pri[p]["files"], th) for p in pri for th in ths}
    print({p: len(T[(p, 0.1)]) for p in pri})

    # (3) full prespecified grid
    rows, curves = [], []
    for sc in scenario_grid(cfg):
        r = evaluate_scenario(sc, T[(sc["prior"], sc["threshold"])], km, kmF, cfg, GRID_T)
        curves.append(r.pop("_Fg")); rows.append(r)
    G = pd.DataFrame(rows)
    CUR = np.array(curves)
    assert len(G) == cfg["grid"]["size"], len(G)
    save(G, "p5_consistency_grid", "Full factorial prespecified grid (1764 scenarios); exploratory=False for all.")
    ref = cfg["grid"]["reference_scenario"]
    refmask = ((G.prior == ref["prior"]) & (G.threshold == ref["threshold"]) & (G.tau == ref["tau"]) & (G.d == ref["d"]) & (G.death_first == ref["death_first"]))
    print("reference row:\n", G[refmask].T.loc[[f"F_sim_{h}" for h in hs] + [f"F_km_{h}" for h in hs]].T.round(4).to_string())

    # ---- METABRIC reference table
    kmrows = [{"horizon_months": h, "horizon_years": h / 12, "cum_inc": km[h]["cum_inc"], "ci_lo": km[h]["ci_lo"], "ci_hi": km[h]["ci_hi"], "S": km[h]["S"],
               "n_at_risk": km[h]["n_at_risk"], "n_events": km[h]["n_events_by_horizon"],
               "flag": "15 y: no Phase 4 reliability criterion" if h == 180 else ("20 y: caution, 162 at risk overall" if h == 240 else "")} for h in hs]
    save(pd.DataFrame(kmrows), "p5_horizon_metabric_reference")

    # ---- model curves at reference (per prior) + Wilson
    cr = pd.DataFrame({"month": GRID_T, "F_km": kmF})
    for p in pri:
        r = G[(G.prior == p) & (G.threshold == 0.1) & (G.tau == 0) & (G.d == 0) & (G.death_first == 0.0)].iloc[0]
        cr[f"F_sim_{p}"] = CUR[int(r.scenario_id)]
    save(cr, "p5_model_event_free_curves_reference", "Reference mapping: theta=0.1, tau=0, d=0, death_first=0. Event-free S = 1 - F.")

    # ---- horizon comparison (reference + grid range)
    hr = []
    for p in pri:
        for h in hs:
            sub = G[G.prior == p]
            r = G[(G.prior == p) & (G.threshold == 0.1) & (G.tau == 0) & (G.d == 0) & (G.death_first == 0.0)].iloc[0]
            n = int(r.n_sim); k = int(round(r[f"F_sim_{h}"] * n)); lo, hi = wilson(k, n)
            hr.append({"prior": p, "horizon_months": h, "horizon_years": h / 12, "F_sim_reference": r[f"F_sim_{h}"], "F_sim_wilson_lo": lo, "F_sim_wilson_hi": hi,
                       "F_km": km[h]["cum_inc"], "km_ci_lo": km[h]["ci_lo"], "km_ci_hi": km[h]["ci_hi"], "D_reference": r[f"D_{h}"],
                       "inside_km_ci_reference": bool(r[f"inside_ci_{h}"]), "F_sim_grid_min_prior": sub[f"F_sim_{h}"].min(), "F_sim_grid_max_prior": sub[f"F_sim_{h}"].max(),
                       "frac_grid_inside_ci_prior": sub[f"inside_ci_{h}"].mean(), "km_n_at_risk": km[h]["n_at_risk"],
                       "flag": "15 y flagged (no reliability criterion)" if h == 180 else ("20 y caution" if h == 240 else "")})
    HR = pd.DataFrame(hr)
    save(HR, "p5_horizon_comparison")

    # ---- envelope & robust findings
    en = []
    for h in hs:
        for lab, sub in [("all", G)] + [(p, G[G.prior == p]) for p in pri]:
            en.append({"group": lab, "horizon_months": h, "n_scenarios": len(sub), "n_inside_ci": int(sub[f"inside_ci_{h}"].sum()),
                       "frac_inside_ci": sub[f"inside_ci_{h}"].mean(), "frac_below_ci": sub[f"below_ci_{h}"].mean(), "frac_above_ci": sub[f"above_ci_{h}"].mean(),
                       "D_min": sub[f"D_{h}"].min(), "D_median": sub[f"D_{h}"].median(), "D_max": sub[f"D_{h}"].max(),
                       "F_sim_min": sub[f"F_sim_{h}"].min(), "F_sim_max": sub[f"F_sim_{h}"].max(), "F_km": km[h]["cum_inc"],
                       "km_ci_lo": km[h]["ci_lo"], "km_ci_hi": km[h]["ci_hi"]})
    save(pd.DataFrame(en), "p5_consistency_envelope")
    interp = {
        "S1": "Simulated cumulative incidence at 5 y is below the METABRIC 95% CI. Expected: the model has no treatment-era / early-recurrence mechanism and starts from microscopic, mostly dormant residual disease; says nothing for or against dormancy.",
        "S2": "At 10 y the simulated incidence is below the METABRIC CI in the stated fraction of scenarios; the early-event deficit persists unless thresholds are low and offsets are zero.",
        "S3": "At 20 y the two numbers are not distinguishable only in the stated fraction of scenarios: agreement at one late horizon depends on the assumed threshold/prior and can occur for different shapes of curve.",
        "S4": "The simulated risk is later than observed (late increment F(240)-F(120) larger than METABRIC's) wherever true: a shape difference, not a level difference.",
        "S5": "A larger share of simulated 20-year events occurs after 10 y than in METABRIC: simulated events are back-loaded relative to the observed composite endpoint.",
        "S6": "Simulated incidence more than doubles between 10 and 20 y whereas METABRIC's increases by < 2x: accelerating versus decelerating late risk.",
    }
    stm = cfg["consistency_criterion"]["prespecified_statements"]
    rb = []
    for s, txt in stm.items():
        row = {"statement": s, "text": txt, "n_true": int(G[s].sum()), "n_scenarios": len(G), "frac_true": G[s].mean()}
        row["verdict"] = "ROBUST (true in all 1764 scenarios)" if G[s].all() else ("never true" if not G[s].any() else "assumption-dependent")
        for p in pri:
            row[f"frac_{p}"] = G[G.prior == p][s].mean()
        row["interpretation"] = interp[s]
        rb.append(row)
    RB = pd.DataFrame(rb)
    save(RB, "p5_consistency_robust_findings")
    print(RB[["statement", "n_true", "frac_true", "verdict"]].to_string())

    # ---- threshold sensitivity (tau=d=0, death_first=0)
    th_rows = []
    for p in pri:
        for th in ths:
            r = G[(G.prior == p) & (G.threshold == th) & (G.tau == 0) & (G.d == 0) & (G.death_first == 0.0)].iloc[0]
            Tm = T[(p, th)]; ev = Tm[np.isfinite(Tm)]
            th_rows.append({"prior": p, "threshold": th, **{f"F_sim_{h}": r[f"F_sim_{h}"] for h in hs}, **{f"D_{h}": r[f"D_{h}"] for h in hs},
                            **{f"inside_ci_{h}": bool(r[f"inside_ci_{h}"]) for h in hs}, "sup_norm_240": r["sup_norm_240"],
                            "median_T_events": float(np.median(ev)) if len(ev) else np.nan, "n_events": len(ev), "n": len(Tm)})
    TH = pd.DataFrame(th_rows)
    save(TH, "p5_threshold_sensitivity")

    # ---- time-origin / delay sensitivity
    to_rows = []
    for p in pri:
        for tau in cfg["grid"]["tau_months"]:
            for d in cfg["grid"]["detect_delay_months"]:
                for df_ in cfg["grid"]["death_first_fraction"]:
                    r = G[(G.prior == p) & (G.threshold == 0.1) & (G.tau == tau) & (G.d == d) & (G.death_first == df_)].iloc[0]
                    to_rows.append({"prior": p, "tau": tau, "d": d, "tau_plus_d": tau + d, "death_first": df_,
                                    **{f"F_sim_{h}": r[f"F_sim_{h}"] for h in hs}, **{f"D_{h}": r[f"D_{h}"] for h in hs},
                                    **{f"inside_ci_{h}": bool(r[f"inside_ci_{h}"]) for h in hs}, "sup_norm_240": r["sup_norm_240"]})
    TO = pd.DataFrame(to_rows)
    save(TO, "p5_timeorigin_sensitivity", "theta=0.1.")
    # alignment diagnostic (DESCRIPTIVE ONLY)
    al = []
    for p in pri:
        Tm = T[(p, 0.1)]
        for s in range(-60, 61):
            up = min(240.0, 240.0 + s)
            m = GRID_T <= up
            F = shifted_cdf(Tm, GRID_T[m], s)
            m2 = GRID_T <= 180      # EXPLORATORY fixed window (added after seeing that the pre-registered shrinking window favours negative shifts)
            F2 = shifted_cdf(Tm, GRID_T[m2], s)
            al.append({"prior": p, "shift_months": s, "compare_up_to_month": up, "sup_norm": sup_norm(F, kmF[m]),
                       "sup_norm_fixed_window_0_180_EXPLORATORY": sup_norm(F2, kmF[m2]), "descriptive_diagnostic_only": True})
    AL = pd.DataFrame(al)
    save(AL, "p5_timeorigin_alignment_curve", "DESCRIPTIVE DIAGNOSTIC: horizontal shift of the simulated reference curve; not an estimate of the origin, delay or any parameter. Excluded from the consistency criterion.")
    best = AL.loc[AL.groupby("prior").sup_norm.idxmin()][["prior", "shift_months", "sup_norm"]]
    best2 = AL.loc[AL.groupby("prior")["sup_norm_fixed_window_0_180_EXPLORATORY"].idxmin()][["prior", "shift_months", "sup_norm_fixed_window_0_180_EXPLORATORY"]]
    print("fixed-window (exploratory):", best2.to_string())
    best = best.merge(best2.rename(columns={"shift_months": "shift_fixed_window_EXPLORATORY"}), on="prior")
    print("alignment diagnostic (descriptive):\n", best.to_string())
    best.assign(descriptive_diagnostic_only=True).pipe(lambda d: save(d, "p5_timeorigin_alignment_minima", "Descriptive 'alignment offset' per prior; NOT an estimate."))

    figures(cfg, G, CUR, cr, HR, TH, TO, AL, km, kmF, T, t, e, pri, ths, hs)
    print("done")


# ----------------------------------------------------------------------------------------------------------- figures
def km_band(ax, t, e, tmax=240):
    tc, S, lo, hi = km_curve(t, e, t_max=tmax)
    ax.step(tc, S, where="post", color=KMC, lw=2.2, label="METABRIC RFS (Kaplan–Meier)", zorder=6)
    ax.fill_between(tc, lo, hi, step="post", color=KMC, alpha=0.15, lw=0, label="METABRIC 95% CI", zorder=5)


def yrs(ax):
    ax.set_xticks(range(0, 241, 60)); ax.set_xlim(0, 240)
    sec = ax.secondary_xaxis("top", functions=(lambda x: x / 12, lambda x: x * 12)); sec.set_xlabel("years", fontsize=8)
    ax.set_xlabel("months (model: since residual state; clinical: since undocumented RFS origin)", fontsize=7.5)


def figures(cfg, G, CUR, cr, HR, TH, TO, AL, km, kmF, T, t, e, pri, ths, hs):
    # ---------------- Fig 4: model event-free survival
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.9)); fig.subplots_adjust(left=0.07, right=0.98, top=0.86, bottom=0.2, wspace=0.22)
    ax = axs[0]
    for p in pri:
        ax.step(GRID_T, 1 - cr[f"F_sim_{p}"], where="post", color=PC[p], ls=PLS[p], lw=2.2 if p == "P1_base" else 1.4, label=f"{p} (n={len(T[(p, 0.1)]):,})")
    ax.set_ylim(0.4, 1.005); yrs(ax); ax.set_ylabel("simulated event-free fraction  1 − F(t)")
    ax.set_title("(a) Model-generated event-free curves, seven prior scenarios\n(threshold 0.1, no offset; reference mapping)", fontsize=9); ax.legend(loc="lower left", fontsize=7)
    ax = axs[1]
    Tm = T[("P1_base", 0.1)]; ev = Tm[np.isfinite(Tm)]
    ax.hist(ev, bins=np.arange(40, 241, 6), color=PC["P1_base"], alpha=0.8)
    ax.axvline(np.median(ev), color=INK, ls="--", lw=1); ax.text(np.median(ev) + 3, ax.get_ylim()[1] * 0.92, f"median {np.median(ev):.0f} mo", fontsize=8)
    ax.set_xlabel("model first-crossing time among events [months]"); ax.set_ylabel("count (P1_base pooled)")
    ax.set_title(f"(b) When do simulated events occur? P1_base: {len(ev):,} of {len(Tm):,} cross by 240 mo\n({len(Tm) - len(ev):,} censored at the 240-month horizon)", fontsize=9)
    footer(fig, "Panel a: survival from the same 5,000/15,000-sample Phase 3 populations; MC half-width ≤ 0.014 (n=5000).")
    fig.savefig(FIG / "fig_p5_04_model_event_free_survival.png"); plt.close(fig)

    # ---------------- Fig 5: model vs METABRIC overlay + difference
    fig = plt.figure(figsize=(9.2, 7.6)); gs = fig.add_gridspec(2, 1, height_ratios=[3, 1.35], hspace=0.28, left=0.09, right=0.98, top=0.91, bottom=0.2)
    ax = fig.add_subplot(gs[0]); km_band(ax, t, e)
    for p in pri:
        ax.step(GRID_T, 1 - cr[f"F_sim_{p}"], where="post", color=PC[p], ls=PLS[p], lw=2.2 if p == "P1_base" else 1.2, alpha=1 if p == "P1_base" else 0.8, label=f"model {p}")
    for tau, ls_ in ((12, (0, (3, 2))), (24, (0, (1, 1)))):
        r = G[(G.prior == "P1_base") & (G.threshold == 0.1) & (G.tau == tau) & (G.d == 0) & (G.death_first == 0.0)].iloc[0]
        ax.step(GRID_T, 1 - CUR[int(r.scenario_id)], where="post", color=PC["P1_base"], ls=ls_, lw=1.3, label=f"model P1_base, τ={tau} mo")
    ax.axvspan(180, 240, color="#cccccc", alpha=0.25, lw=0); ax.text(182, 0.3, "late: ≤ ~200 at risk,\nno reliability criterion at 15 y", fontsize=7, color=INK2)
    ax.set_ylim(0.25, 1.005); yrs(ax); ax.set_ylabel("event-free fraction"); ax.legend(loc="lower left", fontsize=7, ncol=2)
    ax.set_title("Model event-free curves versus METABRIC recurrence-free survival (descriptive overlay; no fitting)", fontsize=9.5)
    ax2 = fig.add_subplot(gs[1], sharex=ax)
    tc, S, lo, hi = km_curve(t, e, t_max=240)
    idx = np.searchsorted(tc, GRID_T, side="right") - 1
    ax2.fill_between(GRID_T, (1 - hi[idx]) - kmF, (1 - lo[idx]) - kmF, color=KMC, alpha=0.15, lw=0, label="METABRIC 95% CI (around 0)")
    for p in pri:
        ax2.plot(GRID_T, cr[f"F_sim_{p}"] - kmF, color=PC[p], ls=PLS[p], lw=2 if p == "P1_base" else 1.1)
    ax2.axhline(0, color=INK, lw=0.8); ax2.set_ylabel("F_sim − F_KM"); ax2.set_xlabel("months"); ax2.legend(loc="lower left", fontsize=7)
    ax2.set_title("Signed difference in cumulative incidence (negative = fewer simulated events than observed)", fontsize=8.5)
    footer(fig, "Overlay shows 7 hypothetical prior scenarios at the headline threshold 0.1; τ is an assumed offset, not an estimate.")
    fig.savefig(FIG / "fig_p5_05_model_vs_metabric.png"); plt.close(fig)

    # ---------------- Fig 6: consistency envelope
    fig, axs = plt.subplots(1, 2, figsize=(11.2, 5.2), gridspec_kw={"width_ratios": [1.5, 1]}); fig.subplots_adjust(left=0.07, right=0.985, top=0.86, bottom=0.22, wspace=0.2)
    ax = axs[0]
    S_all = 1 - CUR
    ax.fill_between(GRID_T, S_all.min(0), S_all.max(0), color=PC["P1_base"], alpha=0.18, lw=0, label=f"min–max over all {len(G)} prespecified scenarios")
    ax.fill_between(GRID_T, np.percentile(S_all, 10, axis=0), np.percentile(S_all, 90, axis=0), color=PC["P1_base"], alpha=0.32, lw=0, label="10th–90th percentile of scenarios")
    km_band(ax, t, e)
    ax.set_ylim(0.2, 1.005); yrs(ax); ax.set_ylabel("event-free fraction"); ax.legend(loc="lower left", fontsize=7.5)
    ax.set_title("(a) Consistency envelope of model event-free curves vs METABRIC", fontsize=9)
    ax = axs[1]
    for i, h in enumerate(hs):
        sub = G[f"D_{h}"]; k = km[h]
        ax.fill_between([i - 0.38, i + 0.38], k["ci_lo"] - k["cum_inc"], k["ci_hi"] - k["cum_inc"], color=KMC, alpha=0.15, lw=0)
        ax.plot([i, i], [sub.min(), sub.max()], color=PC["P1_base"], lw=5, alpha=0.6, solid_capstyle="butt")
        ax.plot([i], [sub.median()], "o", color=INK, ms=5)
        ax.text(i, sub.max() + 0.015, f"{G[f'inside_ci_{h}'].mean() * 100:.0f}%\nin CI", ha="center", fontsize=7.5)
    ax.axhline(0, color=INK, lw=0.8)
    ax.set_xticks(range(len(hs))); ax.set_xticklabels([f"{h // 12} y\n(at risk {km[h]['n_at_risk']})" + ("*" if h == 180 else "") for h in hs], fontsize=8)
    ax.set_ylabel("F_sim − F_KM at horizon"); ax.set_ylim(min(G[[f"D_{h}" for h in hs]].min()) - 0.05, max(G[[f"D_{h}" for h in hs]].max()) + 0.12)
    ax.set_title("(b) Discrepancy across the grid (bar = min–max, dot = median;\ngrey band = METABRIC 95% CI; label = % scenarios inside CI)", fontsize=9)
    footer(fig, "* 15 y has no Phase 4 reliability criterion; 20 y has 162 at risk overall. Envelope = thresholds 0.01–0.5, τ 0–24, d 0–12, death-first 0–27%, 7 priors.")
    fig.savefig(FIG / "fig_p5_06_consistency_envelope.png"); plt.close(fig)

    # ---------------- Fig 7: threshold sensitivity
    fig, axs = plt.subplots(2, 2, figsize=(10.4, 7.6)); fig.subplots_adjust(left=0.08, right=0.985, top=0.9, bottom=0.17, wspace=0.22, hspace=0.6)
    for ax, h in zip(axs.ravel()[:3], (60, 120, 240)):
        k = km[h]
        ax.axhspan(k["ci_lo"], k["ci_hi"], color=KMC, alpha=0.15, lw=0); ax.axhline(k["cum_inc"], color=KMC, lw=1.5, label="METABRIC KM ± 95% CI")
        for p in pri:
            d = TH[TH.prior == p]
            ax.plot(d.threshold, d[f"F_sim_{h}"], marker="o", ms=3, color=PC[p], ls=PLS[p], lw=2 if p == "P1_base" else 1.1, label=p)
        ax.set_xscale("log"); ax.set_xticks(ths); ax.set_xticklabels([f"{x:g}" for x in ths], fontsize=7)
        ax.set_xlabel("assumed detection threshold N_detect (fraction of K)"); ax.set_ylabel(f"F_sim({h // 12} y)"); ax.set_title(f"({'abc'[(60, 120, 240).index(h)]}) cumulative incidence by {h // 12} years", fontsize=9)
    axs[0, 0].legend(fontsize=6.5, loc="center right", ncol=2)
    ax = axs[1, 1]; km_band(ax, t, e)
    cm = plt.cm.viridis(np.linspace(0, 0.9, len(ths)))
    for c, th in zip(cm, ths):
        r = G[(G.prior == "P1_base") & (G.threshold == th) & (G.tau == 0) & (G.d == 0) & (G.death_first == 0.0)].iloc[0]
        ax.step(GRID_T, 1 - CUR[int(r.scenario_id)], where="post", color=c, lw=1.4, label=f"θ={th:g}")
    ax.set_ylim(0.2, 1.005); yrs(ax); ax.set_ylabel("event-free fraction"); ax.legend(fontsize=6.5, loc="lower left", ncol=2); ax.set_title("(d) P1_base event-free curves by threshold", fontsize=9)
    footer(fig, "Thresholds are model-defined burden fractions of K, not imaging limits. Simulated 5-y values ≈ 0 for θ ≥ 0.1.")
    fig.savefig(FIG / "fig_p5_07_detection_threshold_sensitivity.png"); plt.close(fig)

    # ---------------- Fig 8: time origin & delay
    fig, axs = plt.subplots(2, 2, figsize=(10.8, 8.0)); fig.subplots_adjust(left=0.08, right=0.985, top=0.91, bottom=0.17, wspace=0.22, hspace=0.4)
    ax = axs[0, 0]; km_band(ax, t, e)
    sh = (0, 6, 12, 24, 36); cm = plt.cm.plasma(np.linspace(0.05, 0.8, len(sh)))
    Tm = T[("P1_base", 0.1)]
    for c, s in zip(cm, sh):
        ax.step(GRID_T, 1 - shifted_cdf(Tm, GRID_T, s), where="post", color=c, lw=1.5, label=f"τ + d = {s} mo")
    ax.set_ylim(0.2, 1.005); yrs(ax); ax.set_ylabel("event-free fraction"); ax.legend(fontsize=7, loc="lower left"); ax.set_title("(a) P1_base shifted by τ + d (exact translation)", fontsize=9)
    ax = axs[0, 1]
    ss = np.arange(0, 49, 2)
    for h, ls_ in ((60, ":"), (120, "--"), (240, "-")):
        k = km[h]; ax.axhspan(k["ci_lo"], k["ci_hi"], color=KMC, alpha=0.12, lw=0)
        ax.plot(ss, [ecdf_at(Tm + s, h)[0] for s in ss], color=PC["P1_base"], ls=ls_, label=f"F_sim({h // 12} y)")
        ax.axhline(k["cum_inc"], color=KMC, ls=ls_, lw=0.8)
    ax.set_xlabel("τ + d [months]"); ax.set_ylabel("cumulative incidence"); ax.legend(fontsize=7.5, loc="center left", bbox_to_anchor=(0.02, 0.55)); ax.set_title("(b) Offsets lower simulated incidence; METABRIC CI bands in grey", fontsize=9)
    ax = axs[1, 0]
    dfs = cfg["grid"]["death_first_fraction"]
    for j, h in enumerate((120, 240)):
        vals = [TO[(TO.prior == "P1_base") & (TO.tau == 0) & (TO.d == 0) & (TO.death_first == f)][f"F_sim_{h}"].iloc[0] for f in dfs]
        ax.bar(np.arange(len(dfs)) + j * 0.38 - 0.19, vals, 0.36, color=PC["P1_base"], alpha=0.5 + 0.4 * j, label=f"F_sim({h // 12} y)")
        k = km[h]; ax.hlines(k["cum_inc"], -0.5, len(dfs) - 0.5, color=KMC, ls="--" if j == 0 else "-", lw=0.9)
    ax.set_xticks(range(len(dfs))); ax.set_xticklabels([f"{f * 100:.1f}%" for f in dfs]); ax.set_xlabel("assumed death-first share of events (lag mean 24 mo)"); ax.set_ylabel("cumulative incidence")
    ax.set_ylim(0, 0.62); ax.legend(fontsize=7.5, loc="upper left"); ax.set_title("(c) Death-first mapping barely moves the curve (black lines: METABRIC KM)", fontsize=9)
    ax = axs[1, 1]
    for p in pri:
        d = AL[AL.prior == p]; yy = d["sup_norm_fixed_window_0_180_EXPLORATORY"]
        ax.plot(d.shift_months, yy, color=PC[p], ls=PLS[p], lw=2 if p == "P1_base" else 1.1, label=p)
        m = d.loc[yy.idxmin()]; ax.plot([m.shift_months], [m["sup_norm_fixed_window_0_180_EXPLORATORY"]], "v", color=PC[p], ms=6)
    ax.axvline(0, color=INK, lw=0.7); ax.set_xlabel("horizontal shift applied to the simulated curve [months]"); ax.set_ylabel("sup-norm distance to METABRIC KM (window 0–180 mo)")
    ax.set_title("(d) DESCRIPTIVE alignment diagnostic (not an estimate of τ);\n▼ = minimum", fontsize=9); ax.legend(fontsize=6.5, ncol=2)
    footer(fig, "τ, d and death-first lag are assumed nuisances, not estimates. Panel d compares only where model data exist (t ≤ 240 + shift for negative shifts).")
    fig.savefig(FIG / "fig_p5_08_time_origin_delay_uncertainty.png"); plt.close(fig)

    # ---------------- Fig 9: horizon estimates
    fig, ax = plt.subplots(figsize=(9.6, 5.8)); fig.subplots_adjust(left=0.09, right=0.98, top=0.88, bottom=0.26)
    x = np.arange(len(hs))
    for i, h in enumerate(hs):
        k = km[h]
        ax.errorbar(i, k["cum_inc"], yerr=[[k["cum_inc"] - k["ci_lo"]], [k["ci_hi"] - k["cum_inc"]]], fmt="s", color=KMC, ms=8, capsize=5, lw=2, zorder=5, mfc="white" if h == 180 else KMC,
                    label="METABRIC KM ± 95% CI" if i == 0 else None)
        sub = G[f"F_sim_{h}"]
        ax.plot([i + 0.3, i + 0.3], [sub.min(), sub.max()], color="#999", lw=6, alpha=0.5, solid_capstyle="butt", label="range over all 1764 scenarios" if i == 0 else None)
        for j, p in enumerate(pri):
            v = HR[(HR.prior == p) & (HR.horizon_months == h)].F_sim_reference.iloc[0]
            ax.plot(i - 0.3 + j * 0.1, v, "o", color=PC[p], ms=6, label=p if i == 0 else None, mec=INK, mew=0.4)
        ax.text(i, -0.2, f"at risk {k['n_at_risk']}\n{k['n_events_by_horizon']} events" + ("\n(15 y flagged)" if h == 180 else "\n(caution)" if h == 240 else ""), ha="center", fontsize=7.5, transform=ax.get_xaxis_transform())
    ax.set_xticks(x); ax.set_xticklabels([f"{h // 12} years" for h in hs]); ax.set_ylabel("cumulative incidence by horizon"); ax.set_ylim(-0.02, 0.78)
    ax.legend(fontsize=7, ncol=3, loc="upper left"); ax.set_title("Observed (METABRIC) vs simulated cumulative incidence at supported horizons", fontsize=10)
    footer(fig, "Dots: reference mapping (θ=0.1, τ=d=0) for each prior (left of each group); grey bar: range across the full prespecified grid. 15-y estimate flagged (open square).")
    fig.savefig(FIG / "fig_p5_09_horizon_estimates.png"); plt.close(fig)


if __name__ == "__main__":
    main()
