"""Phase 5 (P5a): descriptive METABRIC Kaplan-Meier analysis.  NO model is fitted, overlaid or compared here.

Run: .venv/bin/python experiments/phase5_metabric_km.py
Reads the raw cBioPortal export (never modified).  Endpoint: RFS_STATUS/RFS_MONTHS as exported (composite of loco-regional or
distant relapse and disease-specific death; time origin undocumented -- see results/p5_rfs_time_origin_note.md).
Outputs: results/tables/p5_km_overall.csv (written FIRST), p5_km_subgroups.csv, p5_logrank_subgroups.csv, p5_landmark_5y.csv,
results/p5_rfs_time_origin_note.md, p5_clinical_coding_checks.csv, p5_late_recurrence_hazard.csv, figures fig_p5_01_metabric_km_overall.png, fig_p5_02_er_positive_vs_negative.png, fig_p5_03_molecular_subgroups.png.
"""
import hashlib
from scipy import stats
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
    coding_checks(d)

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

    late_hazard(d, masks)
    write_note(d, overall, dif)
    fig_overall(t, e, overall)
    fig_er(d, t, e, masks, lr)
    fig_pam50(t, e, masks, reliable)


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


def _short(g):
    return g.replace("PAM50: ", "").replace("ER (IHC) ", "ER ")


def _at_risk_rows(ax, t, strata, years, colors, maxy=None):
    """Number-at-risk rows under an axis: strata = {label: boolean mask}; maxy = {label: last horizon shown}."""
    ax.set_ylim(-0.5, len(strata) - 0.5); ax.set_yticks([]); ax.grid(False)
    ax.spines["left"].set_visible(False)
    for r, (lab, m) in enumerate(strata.items()):
        row = len(strata) - 1 - r
        ax.text(-0.015, row, lab, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8, color=colors[r], fontweight="bold")
        for y in years:
            if maxy is None or y <= maxy[lab]:
                ax.text(y, row, f"{int(np.sum(t[m] >= 12 * y)):,}", ha="center", va="center", fontsize=8.5)
            else:
                ax.text(y, row, "–", ha="center", va="center", fontsize=8.5, color=INK2)


def fig_er(d, t, e, masks, lr):
    pos, neg = masks["ER (IHC) positive"].to_numpy(), masks["ER (IHC) negative"].to_numpy()
    fig = plt.figure(figsize=(7.8, 6.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[4.2, 1], hspace=0.12, left=0.17, right=0.97, top=0.93, bottom=0.17)
    ax = fig.add_subplot(gs[0]); axt = fig.add_subplot(gs[1], sharex=ax)
    for m, lab, c, ls in ((pos, "ER-positive (IHC)", CB["blue"], "-"), (neg, "ER-negative (IHC)", CB["verm"], "--")):
        tc, S, lo, hi = km_curve(t[m], e[m], t_max=240)
        ax.step(tc / 12, 1 - S, where="post", color=c, ls=ls, lw=2, label=f"{lab} (n = {int(m.sum()):,}; {int(e[m].sum())} events)")
        ax.fill_between(tc / 12, 1 - hi, 1 - lo, step="post", color=c, alpha=0.15, lw=0)
    sel = lr[(lr.family == "ER") & (lr.test == "global")].set_index("horizon_years")
    parts = []
    for y in (5, 10, 20):
        if y in sel.index and np.isfinite(sel.loc[y, "p_raw"]):
            p = sel.loc[y, "p_raw"]
            parts.append(f"{y} y: p {'< 0.001' if p < 1e-3 else '= %.3f' % p}")
    ax.text(0.98, 0.04, "Log-rank (censored at horizon) " + "; ".join(parts), transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color=INK2)
    ax.set_ylabel("Cumulative incidence of RFS event (1 − S)"); ax.set_ylim(0, 0.85); ax.legend(loc="upper left")
    ax.set_title("METABRIC RFS event by ER (IHC) status, 0–20 years (95% log-log CI)")
    plt.setp(ax.get_xticklabels(), visible=False)
    _at_risk_rows(axt, t, {"ER+": pos, "ER−": neg}, range(0, 21, 5), [CB["blue"], CB["verm"]])
    axt.set_xlabel("Years since the RFS time origin (undocumented in the export)"); axt.set_xlim(-0.5, 20.5)
    style_fig_footer(fig, "At-risk counts at 0, 5, 10, 15, 20 y. Association only (treatment, stage, grade not adjusted). Composite endpoint (relapse or "
                          "disease-specific death). " + CLIN)
    fig.savefig(FIG / "fig_p5_02_er_positive_vs_negative.png"); plt.close(fig)


def fig_pam50(t, e, masks, reliable):
    gl = [g for g in masks if g.startswith("PAM50")]
    cols = [CB["blue"], CB["orange"], CB["green"], CB["black"], CB["purple"], CB["sky"]]
    ls = ["-", "--", "-.", ":", "-", "--"]
    fig = plt.figure(figsize=(8.4, 6.8))
    gs = fig.add_gridspec(2, 1, height_ratios=[4, 1.5], hspace=0.12, left=0.20, right=0.97, top=0.93, bottom=0.17)
    ax = fig.add_subplot(gs[0]); axt = fig.add_subplot(gs[1], sharex=ax)
    strata, colors, maxy = {}, [], {}
    for i, g in enumerate(gl):
        my = max([y for y in RELIABLE_HORIZONS if reliable[(g, y)]], default=0)
        if my == 0:
            continue
        m = masks[g].to_numpy()
        tc, S, lo, hi = km_curve(t[m], e[m], t_max=12 * my)
        ax.step(tc / 12, 1 - S, where="post", color=cols[i], ls=ls[i], lw=1.9, label=f"{_short(g)} (n = {int(m.sum())}; shown to {my} y)")
        ax.fill_between(tc / 12, 1 - hi, 1 - lo, step="post", color=cols[i], alpha=0.10, lw=0)
        strata[_short(g)] = m; colors.append(cols[i]); maxy[_short(g)] = my
    ax.set_ylabel("Cumulative incidence of RFS event (1 − S)"); ax.set_ylim(0, 0.9); ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("METABRIC RFS event by PAM50 / claudin-low subtype (reliable strata and horizons only)")
    plt.setp(ax.get_xticklabels(), visible=False)
    _at_risk_rows(axt, t, strata, range(0, 21, 5), colors, maxy)
    axt.set_xlabel("Years since the RFS time origin (undocumented in the export)"); axt.set_xlim(-0.5, 20.5)
    style_fig_footer(fig, "Curves end at the longest horizon (5, 10 or 20 y) at which the Phase 4 reliability criterion holds (>=50 events in group; >=100 at risk at "
                          "5 and 10 y; >=30 at 20 y); '–' = not shown. NC (n = 6) and missing subtype omitted. Unadjusted associations. " + CLIN)
    fig.savefig(FIG / "fig_p5_03_molecular_subgroups.png"); plt.close(fig)


def coding_checks(d):
    raw = pd.read_csv(RAW / "data_clinical_patient.txt", sep="\t", comment="#", dtype=str)
    st = raw.RFS_STATUS.value_counts(dropna=False)
    rfs = pd.to_numeric(raw.RFS_MONTHS, errors="coerce")
    os_ = pd.to_numeric(raw.OS_MONTHS, errors="coerce")
    both = rfs.notna() & os_.notna()
    ph1 = pd.read_csv(TAB / "metabric_km_overall.csv", comment="#")
    rows = []

    def add(name, value, expected, ok, note=""):
        rows.append({"check": name, "value": value, "expected": expected, "pass": bool(ok), "note": note})
    add("patients in export", len(raw), 2509, len(raw) == 2509)
    add("RFS_STATUS levels", "; ".join(f"{k}={v}" for k, v in st.items()), "0:Not Recurred, 1:Recurred, missing",
        set(raw.RFS_STATUS.dropna().unique()) == {"0:Not Recurred", "1:Recurred"}, "event indicator = leading digit")
    add("events among usable (status startswith '1')", int(d.event.sum()), 964, int(d.event.sum()) == 964)
    add("usable (status and time non-missing)", len(d), 2380, len(d) == 2380, "Phase 1/4 cohort")
    add("RFS_STATUS missing", int(raw.RFS_STATUS.isna().sum()), 21, int(raw.RFS_STATUS.isna().sum()) == 21)
    add("RFS_MONTHS missing", int(rfs.isna().sum()), "", True, "excluded from the usable cohort")
    add("negative RFS_MONTHS among usable", int((d.RFS_MONTHS < 0).sum()), 0, (d.RFS_MONTHS < 0).sum() == 0)
    add("RFS_MONTHS == 0 among usable (kept, as in Phase 1)", int((d.RFS_MONTHS == 0).sum()), "", True, "zero-time records never enter a risk set after t = 0; kept so the cohort equals the Phase 1/4 cohort")
    n_viol = int((rfs[both] > os_[both] + 1e-6).sum())
    add("RFS > OS among patients with both times", n_viol, 0, n_viol == 0, "consistent with a common clock; does not identify the RFS origin")
    add("max RFS_MONTHS among usable", round(float(d.RFS_MONTHS.max()), 2), "", True, "months")
    add("events with RFS_MONTHS > 240", int(((d.RFS_MONTHS > 240) & (d.event == 1)).sum()), "", True, "horizons censor administratively")
    tt, ee = d.RFS_MONTHS.to_numpy(), d.event.to_numpy()
    for y in HORIZONS:
        r = ph1[ph1.horizon_years == y].iloc[0]
        ci = km_point(tt, ee, 12 * y)["cum_inc"]
        add(f"cumulative incidence at {y} y vs Phase 1 table", round(ci, 4), round(r.cumulative_incidence_1_minus_S, 4),
            abs(ci - r.cumulative_incidence_1_minus_S) < 5e-4, f"|diff| = {abs(ci - r.cumulative_incidence_1_minus_S):.1e}")
    df = pd.DataFrame(rows)
    save(df, "p5_clinical_coding_checks")
    print(df[["check", "value", "expected", "pass"]].to_string())
    assert df["pass"].all(), "coding check failed"


def poisson_ci(k, T, conf=0.95):
    a = 1 - conf
    lo = 0.0 if k == 0 else stats.chi2.ppf(a / 2, 2 * k) / 2
    hi = stats.chi2.ppf(1 - a / 2, 2 * k + 2) / 2
    return lo / T, hi / T


def late_hazard(d, masks):
    t = d.RFS_MONTHS.to_numpy(); e = d.event.to_numpy()
    wins = [(0, 5), (5, 10), (10, 15), (15, 20), (5, 20)]
    groups = {"all usable": np.ones(len(d), bool), "ER positive": masks["ER (IHC) positive"].to_numpy(), "ER negative": masks["ER (IHC) negative"].to_numpy()}
    rows = []
    for g, m in groups.items():
        for a, b in wins:
            lo_m, hi_m = 12 * a, 12 * b
            py = float(np.sum(np.clip(np.minimum(t[m], hi_m) - lo_m, 0, None))) / 12
            ev = int(np.sum(m & (e == 1) & (t > lo_m) & (t <= hi_m)))
            n0 = int(np.sum(t[m] > lo_m))
            lo, hi = poisson_ci(ev, py) if py > 0 else (np.nan, np.nan)
            rows.append({"group": g, "window_years": f"{a}-{b}",
                         "note": "5-y landmark window: conditional on being event-free and observed at 5 y" if (a, b) == (5, 20) else "",
                         "n_under_observation_at_window_start": n0, "events": ev, "person_years": py,
                         "hazard_per_100py": 100 * ev / py if py > 0 else np.nan, "ci_lo_per_100py": 100 * lo, "ci_hi_per_100py": 100 * hi})
    out = pd.DataFrame(rows)
    save(out, "p5_late_recurrence_hazard")
    print(out.round(3).to_string())
    w = pd.read_csv(TAB / "metabric_events_by_window.csv", comment="#").iloc[:4]
    a = out[out.group == "all usable"].set_index("window_years")
    keys = ("0-5", "5-10", "10-15", "15-20")
    cmp = pd.DataFrame({"window_years": keys, "phase1_events": w.events.to_numpy(), "phase5_events": [a.loc[x, "events"] for x in keys],
                        "phase1_records": w.n_records.to_numpy(), "phase5_under_observation_at_start": [a.loc[x, "n_under_observation_at_window_start"] for x in keys]})
    save(cmp, "p5_late_vs_phase1_windows")
    print(cmp.to_string())
    # death-vs-recurrence bounds (reused Phase 4 table; no competing-risk analysis possible)
    b = pd.read_csv(TAB / "p4_metabric_death_vs_recurrence_bounds.csv", comment="#")
    save(b, "p5_death_vs_recurrence_bounds")


if __name__ == "__main__":
    main()
