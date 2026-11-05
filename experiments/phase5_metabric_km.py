"""Phase 5 (P5a): descriptive METABRIC Kaplan-Meier analysis.  NO model is fitted, overlaid or compared here.

Run: .venv/bin/python experiments/phase5_metabric_km.py
Reads the raw cBioPortal export (never modified).  Endpoint: RFS_STATUS/RFS_MONTHS as exported (composite of loco-regional or
distant relapse and disease-specific death; time origin undocumented -- see results/p5_rfs_time_origin_note.md).
Outputs: results/tables/p5_km_overall.csv (written FIRST), p5_km_subgroups.csv, p5_logrank_subgroups.csv, p5_landmark_5y.csv,
results/p5_rfs_time_origin_note.md, figures fig_p5_01_metabric_km_overall.png, fig_p5_02_metabric_km_subgroups.png.
"""
import hashlib
import sys
import textwrap
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from _style import *  # noqa: F401,F403
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from dormancy.km_descriptive import (km_point, km_curve, logrank_test, pairwise_logrank, landmark_incidence, reliable_horizon_mask)

RAW = ROOT / "data" / "raw" / "brca_metabric"
RES = ROOT / "results"
HORIZONS = (5, 10, 15, 20)               # years; months = 12 * years
RELIABLE_HORIZONS = (5, 10, 20)          # horizons for which Phase 4 defined a reliability criterion
CLIN = "CLINICAL (descriptive): cBioPortal METABRIC export, RFS as defined there (time origin undocumented); no model fitted or overlaid."
CB = {"black": "#000000", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442", "blue": "#0072B2",
      "verm": "#D55E00", "purple": "#CC79A7"}  # Okabe-Ito


def provenance():
    h = hashlib.sha256((RAW / "data_clinical_patient.txt").read_bytes()).hexdigest()[:12]
    return f"# phase5 descriptive METABRIC; raw patient file sha256[:12]={h}; RFS time origin undocumented; not model validation"


def save(df, name):
    path = TAB / f"{name}.csv"
    with open(path, "w") as fh:
        fh.write(provenance() + "\n")
        df.to_csv(fh, index=False)
    return path


def load():
    pat = pd.read_csv(RAW / "data_clinical_patient.txt", sep="\t", comment="#", dtype=str)
    for c in ("OS_MONTHS", "RFS_MONTHS"):
        pat[c] = pd.to_numeric(pat[c], errors="coerce")
    pat["event"] = pat.RFS_STATUS.str.startswith("1").astype(float).where(pat.RFS_STATUS.notna())
    both = pat.dropna(subset=["RFS_MONTHS", "event"]).copy()
    both["event"] = both["event"].astype(int)
    return both


def group_masks(d):
    return {"ER (IHC) positive": d.ER_IHC == "Positve",  # sic: spelling in the export
            "ER (IHC) negative": d.ER_IHC == "Negative",
            "PAM50: LumA": d.CLAUDIN_SUBTYPE == "LumA", "PAM50: LumB": d.CLAUDIN_SUBTYPE == "LumB",
            "PAM50: Her2": d.CLAUDIN_SUBTYPE == "Her2", "PAM50: Basal": d.CLAUDIN_SUBTYPE == "Basal",
            "PAM50: claudin-low": d.CLAUDIN_SUBTYPE == "claudin-low", "PAM50: Normal": d.CLAUDIN_SUBTYPE == "Normal"}


def main():
    d = load()
    t, e = d.RFS_MONTHS.to_numpy(), d.event.to_numpy()
    print(f"usable n={len(d)} events={int(e.sum())}")

    # ---------------------------------------------------------------- 1. overall (written first)
    rows = []
    for y in HORIZONS:
        r = km_point(t, e, 12 * y)
        rows.append({"horizon_years": y, "n_at_risk": r["n_at_risk"], "n_events_by_horizon": r["n_events_by_horizon"],
                     "cum_inc": r["cum_inc"], "ci_lo": r["ci_lo"], "ci_hi": r["ci_hi"], "S": r["S"], "S_lo": r["S_lo"], "S_hi": r["S_hi"], "se_S": r["se"]})
    overall = pd.DataFrame(rows)
    save(overall, "p5_km_overall")
    print(overall.round(4).to_string())
    ph1 = pd.read_csv(TAB / "metabric_km_overall.csv", comment="#")
    dif = np.abs(ph1.set_index("horizon_years").cumulative_incidence_1_minus_S - overall.set_index("horizon_years").cum_inc).max()
    print(f"max |difference| vs Phase 1 metabric_km_overall.csv cumulative incidence: {dif:.2e}")

    # ---------------------------------------------------------------- 2. subgroups: reliability filter
    masks = group_masks(d)
    rel = pd.read_csv(TAB / "p4_metabric_subgroup_reliability.csv", comment="#").set_index("group")
    rel_cols = {5: "reliable_5y (>=50 events, >=100 at risk at 60 mo)", 10: "reliable_10y (>=50 events, >=100 at risk at 120 mo)",
                20: "reliable_20y (>=50 events, >=30 at risk at 240 mo)"}
    groups = {"all usable": pd.Series(True, index=d.index), **masks}
    sub_rows, reliable = [], {}
    for g, m in groups.items():
        tg, eg = t[m.to_numpy()], e[m.to_numpy()]
        for y in RELIABLE_HORIZONS:
            r = km_point(tg, eg, 12 * y)
            # recompute the Phase 4 criterion here and check it against the Phase 4 table
            ok = reliable_horizon_mask(int(eg.sum()), r["n_at_risk"], y)
            assert ok == bool(rel.loc[g, rel_cols[y]]), (g, y)
            reliable[(g, y)] = ok
            reason = "" if ok else ("fewer than 50 events in group" if eg.sum() < 50 else
                                    f"at risk at {12 * y} mo = {r['n_at_risk']} < {100 if y < 20 else 30}")
            sub_rows.append({"group": g, "horizon_years": y, "n_group": int(m.sum()), "n_events_total": int(eg.sum()),
                             "n_at_risk": r["n_at_risk"], "n_events_by_horizon": r["n_events_by_horizon"], "reliable": ok,
                             "exclusion_reason": reason,
                             "cum_inc": r["cum_inc"] if ok else np.nan, "ci_lo": r["ci_lo"] if ok else np.nan, "ci_hi": r["ci_hi"] if ok else np.nan,
                             "S": r["S"] if ok else np.nan})
    sub = pd.DataFrame(sub_rows)
    save(sub, "p5_km_subgroups")
    print(sub[["group", "horizon_years", "n_at_risk", "reliable", "cum_inc", "ci_lo", "ci_hi"]].round(3).to_string())

    # ---------------------------------------------------------------- 3. log-rank within families (Holm within family x horizon)
    families = {"ER": ["ER (IHC) positive", "ER (IHC) negative"],
                "PAM50": [g for g in masks if g.startswith("PAM50")]}
    lr_rows = []
    for fam, gl in families.items():
        for y in RELIABLE_HORIZONS:
            use = [g for g in gl if reliable[(g, y)]]
            if len(use) < 2:
                lr_rows.append({"family": fam, "horizon_years": y, "test": "global", "groups": "; ".join(use), "n_groups": len(use), "chi2": np.nan,
                                "df": np.nan, "p_raw": np.nan, "p_holm": np.nan, "note": "fewer than two reliable strata"})
                continue
            lab = np.full(len(d), "", dtype=object)
            for g in use:
                lab[masks[g].to_numpy()] = g
            sel = lab != ""
            chi2, df, p, O, E = logrank_test(t[sel], e[sel], lab[sel], horizon=12 * y)
            lr_rows.append({"family": fam, "horizon_years": y, "test": "global", "groups": "; ".join(use), "n_groups": len(use), "chi2": chi2, "df": df,
                            "p_raw": p, "p_holm": np.nan, "note": "log-rank, administratively censored at horizon; strata restricted to those reliable at the horizon"})
            if len(use) > 2:
                for pr in pairwise_logrank(t[sel], e[sel], lab[sel], horizon=12 * y):
                    lr_rows.append({"family": fam, "horizon_years": y, "test": "pairwise", "groups": f"{pr['group_a']} vs {pr['group_b']}", "n_groups": 2,
                                    "chi2": pr["chi2"], "df": 1, "p_raw": pr["p_raw"], "p_holm": pr["p_holm"],
                                    "note": f"Holm over {len(use) * (len(use) - 1) // 2} pairwise tests in this family/horizon"})
    lr = pd.DataFrame(lr_rows)
    save(lr, "p5_logrank_subgroups")
    print(lr[["family", "horizon_years", "test", "groups", "chi2", "p_raw", "p_holm"]].round(4).to_string())

    # ---------------------------------------------------------------- 4. 5-year landmark
    lm_rows = []
    for g, m in groups.items():
        tg, eg = t[m.to_numpy()], e[m.to_numpy()]
        for ahead_y in (5, 10, 15):
            abs_y = 5 + ahead_y
            if abs_y in RELIABLE_HORIZONS:
                ok = reliable[(g, abs_y)]; why = "" if ok else f"group not reliable at {abs_y} y (Phase 4 criterion)"
            else:  # 15 y: no Phase 4 criterion; overall cohort only
                ok = (g == "all usable"); why = "" if ok else "no Phase 4 reliability criterion at 15 y; reported for the whole cohort only"
            r = landmark_incidence(tg, eg, 60.0, 12 * ahead_y)
            lm_rows.append({"group": g, "landmark_years": 5, "ahead_years": ahead_y, "absolute_horizon_years": abs_y,
                            "n_event_free_at_landmark": r["n_landmark"], "n_at_risk_at_horizon": r["n_at_risk"], "n_events_in_window": r["n_events_by_horizon"],
                            "reported": ok, "exclusion_reason": why,
                            "cond_cum_inc": r["cum_inc"] if ok else np.nan, "ci_lo": r["ci_lo"] if ok else np.nan, "ci_hi": r["ci_hi"] if ok else np.nan})
    lm = pd.DataFrame(lm_rows)
    save(lm, "p5_landmark_5y")
    print(lm[lm.reported][["group", "ahead_years", "n_event_free_at_landmark", "n_events_in_window", "cond_cum_inc", "ci_lo", "ci_hi"]].round(3).to_string())

    write_note(d, overall, dif)
    fig_overall(t, e, overall)
    fig_subgroups(d, t, e, masks, reliable)


# ---------------------------------------------------------------------------------------------------------------- note
def write_note(d, overall, dif):
    txt = f"""# RFS time-origin limitation (Phase 5 descriptive METABRIC analysis)

**Verified (documentation / exact computation, Phase 4):**
- RFS_STATUS is a composite endpoint (loco-regional relapse, distant relapse or disease-specific death); the cBioPortal export does not
  provide the time of each component, nor which event came first.
- OS_MONTHS is documented as measured from initial diagnosis. **No documentation of the RFS_MONTHS origin was found** (export metadata,
  open text of the METABRIC papers, datahub history). The export has no treatment or surgery dates.

**Inference (not proof of the origin):** RFS_MONTHS equals OS_MONTHS exactly for all 745 patients not recurred and living and all 432 who
died of other causes (0 violations of RFS <= OS in 1,980 patients with both times), which suggests a common clock for the two
endpoints. This does not establish whether the clock starts at diagnosis, surgery or another date.

**Consequence for this analysis:** every Kaplan-Meier time in Phase 5 is "months since the (unknown) RFS origin". The offset between
diagnosis and the RFS origin, the diagnosis-to-treatment interval and any detection delay are nuisance quantities (Phase 4, section 5-6)
and are NOT estimated here. Results are descriptive of the exported variable; they are not comparable to model times without a
stated time-origin assumption, and they are not model validation.

Analysis sample: {len(d)} patients with usable RFS time and status ({int(d.event.sum())} events). Cumulative incidence in this file agrees with the
Phase 1 table results/tables/metabric_km_overall.csv to within {dif:.1e}.
"""
    (RES / "p5_rfs_time_origin_note.md").write_text(txt)


# --------------------------------------------------------------------------------------------------------------- figures
def style_fig_footer(fig, kind_txt):
    lines = textwrap.wrap(kind_txt, max(60, int(fig.get_figwidth() * 19)))
    fig.text(0.01, 0.005, "\n".join(lines), fontsize=6.5, color="#1c4f9b", ha="left", va="bottom")
    fig.text(0.995, 0.995, "CLINICAL (descriptive)", fontsize=8, fontweight="bold", color="#1c4f9b", ha="right", va="top")


def fig_overall(t, e, overall):
    fig = plt.figure(figsize=(7.6, 6.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[4, 1], hspace=0.08, left=0.10, right=0.97, top=0.92, bottom=0.17)
    ax = fig.add_subplot(gs[0]); axt = fig.add_subplot(gs[1], sharex=ax)
    tc, S, lo, hi = km_curve(t, e, t_max=300)
    ax.step(tc / 12, 1 - S, where="post", color=CB["blue"], lw=2, label="Kaplan–Meier cumulative incidence (1 − S)")
    ax.fill_between(tc / 12, 1 - hi, 1 - lo, step="post", color=CB["blue"], alpha=0.20, lw=0, label="95% pointwise CI (log–log Greenwood)")
    for _, r in overall.iterrows():
        ax.plot([r.horizon_years], [r.cum_inc], "o", color=CB["verm"], ms=6, zorder=5)
        ax.annotate(f"{r.cum_inc:.3f}", (r.horizon_years, r.cum_inc), textcoords="offset points", xytext=(6, -12), fontsize=8, color=INK)
    ax.set_ylabel("Cumulative incidence of RFS event"); ax.set_ylim(0, 0.75)
    ax.set_title("METABRIC (cBioPortal): recurrence-free survival event, all usable patients (n = 2,380; 964 events)")
    ax.legend(loc="upper left"); plt.setp(ax.get_xticklabels(), visible=False)
    axt.set_ylim(0, 1); axt.set_yticks([]); axt.grid(False)
    for s in ("left",):
        axt.spines[s].set_visible(False)
    for y in range(0, 26, 5):
        n = int(np.sum(t >= 12 * y))
        axt.text(y, 0.5, f"{n:,}", ha="center", va="center", fontsize=8.5)
    axt.text(-0.9, 0.5, "At risk", ha="right", va="center", fontsize=8.5)
    axt.set_xlabel("Years since the RFS time origin (undocumented in the export)")
    axt.set_xlim(-1.0, 25)
    style_fig_footer(fig, "Markers at 5, 10, 15, 20 y. Greenwood variance, log-log CI. " + CLIN)
    fig.savefig(FIG / "fig_p5_01_metabric_km_overall.png"); plt.close(fig)


def fig_subgroups(d, t, e, masks, reliable):
    fams = [("ER (IHC) status", ["ER (IHC) positive", "ER (IHC) negative"], [CB["blue"], CB["verm"]]),
            ("PAM50 / claudin-low subtype", [g for g in masks if g.startswith("PAM50")],
             [CB["blue"], CB["orange"], CB["green"], CB["black"], CB["purple"], CB["sky"]])]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.2), sharey=True)
    fig.subplots_adjust(left=0.07, right=0.985, top=0.88, bottom=0.30, wspace=0.08)
    ls = ["-", "--", "-.", ":", "-", "--"]
    for ax, (title, gl, cols) in zip(axes, fams):
        for i, g in enumerate(gl):
            maxy = max([y for y in RELIABLE_HORIZONS if reliable[(g, y)]], default=0)
            if maxy == 0:
                continue
            m = masks[g].to_numpy()
            tc, S, lo, hi = km_curve(t[m], e[m], t_max=12 * maxy)
            ax.step(tc / 12, 1 - S, where="post", color=cols[i], ls=ls[i], lw=1.8, label=f"{g.replace('PAM50: ', '').replace('ER (IHC) ', 'ER ')} (reliable to {maxy} y)")
            ax.fill_between(tc / 12, 1 - hi, 1 - lo, step="post", color=cols[i], alpha=0.10, lw=0)
            ax.plot([tc[-1] / 12], [1 - S[-1]], "o", color=cols[i], ms=4)
        ax.set_title(title); ax.set_xlabel("Years since the RFS time origin (undocumented)"); ax.set_xlim(0, 20.3)
        ax.legend(loc="upper left", fontsize=7.5)
    axes[0].set_ylabel("Cumulative incidence of RFS event"); axes[0].set_ylim(0, 0.9)
    # at-risk tables below
    for ax, (title, gl, cols) in zip(axes, fams):
        txt = ["At risk (5 / 10 / 15 / 20 y):"]
        for g in gl:
            m = masks[g].to_numpy()
            nr = []
            for y in (5, 10, 15, 20):
                ok = y == 15 and any(reliable[(g, yy)] for yy in (20,)) or (y != 15 and reliable[(g, y)])
                nr.append(f"{int(np.sum(t[m] >= 12 * y))}" if ok else "–")
            txt.append(f"{g.replace('PAM50: ', '').replace('ER (IHC) ', 'ER ')}: " + " / ".join(nr))
        ax.text(0.0, -0.20, "\n".join(txt), transform=ax.transAxes, fontsize=7.2, va="top", family="monospace")
    style_fig_footer(fig, "Curves are drawn only up to the longest horizon (5, 10 or 20 y) at which the Phase 4 reliability criterion holds "
                          "(≥50 events in the group; ≥100 at risk at 5 and 10 y, ≥30 at 20 y); '–' = criterion not met. Strata overlap between panels; "
                          "differences are associations, not mechanisms. " + CLIN)
    fig.savefig(FIG / "fig_p5_02_metabric_km_subgroups.png"); plt.close(fig)


if __name__ == "__main__":
    main()
