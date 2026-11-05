"""Phase 4, Task 5: what the METABRIC export supports about the RFS endpoint (descriptive; NO model is fitted or overlaid).

Run: .venv/bin/python experiments/metabric_endpoint_investigation.py
Reads the raw cBioPortal export only (data/raw/brca_metabric/). Outputs results/tables/p4_metabric_*.csv.
Every statement is classified VERIFIED (stated in documentation or computed exactly from the files) or INFERENCE (a data-based
argument that does not prove the claim). The time origin of RFS_MONTHS is NOT inferred from RFS_MONTHS <= OS_MONTHS.
"""
from _p4common import *  # noqa: F401,F403

RAW = ROOT / "data" / "raw" / "brca_metabric"


def load():
    pat = pd.read_csv(RAW / "data_clinical_patient.txt", sep="\t", comment="#", dtype=str)
    smp = pd.read_csv(RAW / "data_clinical_sample.txt", sep="\t", comment="#", dtype=str)
    for c in ("OS_MONTHS", "RFS_MONTHS", "AGE_AT_DIAGNOSIS"):
        pat[c] = pd.to_numeric(pat[c], errors="coerce")
    pat["rfs_event"] = pat.RFS_STATUS.str.startswith("1").astype(float).where(pat.RFS_STATUS.notna())
    return pat, smp


def main():
    pat, smp = load()
    both = pat.dropna(subset=["RFS_MONTHS", "rfs_event"]).copy()
    log(f"usable RFS records: {len(both)} (events {int(both.rfs_event.sum())})")
    tol = 1e-6

    # ------------------------------------------------------------------ 1. common-origin tests (stronger than RFS <= OS)
    rows = []

    def eq_stats(label, sub, rationale):
        s = sub.dropna(subset=["OS_MONTHS", "RFS_MONTHS"])
        d = s.OS_MONTHS - s.RFS_MONTHS
        rows.append({"group": label, "n_with_both_times": len(s), "n_RFS_equals_OS": int((d.abs() < tol).sum()),
                     "fraction_equal": float((d.abs() < tol).mean()) if len(s) else float("nan"), "n_RFS_less_than_OS": int((d > tol).sum()),
                     "n_RFS_greater_than_OS (violations)": int((d < -tol).sum()), "median_OS_minus_RFS_when_less": float(d[d > tol].median()) if (d > tol).any() else float("nan"),
                     "rationale": rationale})

    eq_stats("Not Recurred & Living", both[(both.rfs_event == 0) & (both.VITAL_STATUS == "Living")],
             "both clocks end at the last contact date; equality requires the same origin (does not say WHICH origin)")
    eq_stats("Not Recurred & Died of Other Causes", both[(both.rfs_event == 0) & (both.VITAL_STATUS == "Died of Other Causes")],
             "both clocks end at death by other causes (RFS censored there); equality requires same origin")
    eq_stats("Recurred & Died of Disease", both[(both.rfs_event == 1) & (both.VITAL_STATUS == "Died of Disease")],
             "equality = death recorded as the first event (recurrence not separately recorded); inequality = recurrence preceded death")
    eq_stats("Recurred & Living", both[(both.rfs_event == 1) & (both.VITAL_STATUS == "Living")], "recurrence while alive: RFS < OS expected")
    eq_stats("Recurred & Died of Other Causes", both[(both.rfs_event == 1) & (both.VITAL_STATUS == "Died of Other Causes")], "recurrence then death by other cause")
    eq_stats("All with both times", both, "(for reference)")
    org = pd.DataFrame(rows)
    save_table(org, "p4_metabric_origin_tests")
    print(org[["group", "n_with_both_times", "n_RFS_equals_OS", "fraction_equal", "n_RFS_greater_than_OS (violations)"]].round(4).to_string())

    # ------------------------------------------------------------------ 2. composition of the composite event
    ev = both[both.rfs_event == 1].dropna(subset=["OS_MONTHS"])
    ev_all = both[both.rfs_event == 1]
    comp = []
    n_ev = len(ev_all)
    cat = {
        "Recurred & Living (a recurrence; no death)": (ev_all.VITAL_STATUS == "Living"),
        "Recurred & Died of Other Causes (recurrence recorded; death by other cause)": (ev_all.VITAL_STATUS == "Died of Other Causes"),
        "Died of Disease & RFS < OS (recurrence preceded disease-specific death)": (ev_all.VITAL_STATUS == "Died of Disease") & ((ev_all.OS_MONTHS - ev_all.RFS_MONTHS) > tol),
        "Died of Disease & RFS = OS (death recorded as first event; recurrence may be unrecorded)": (ev_all.VITAL_STATUS == "Died of Disease") & ((ev_all.OS_MONTHS - ev_all.RFS_MONTHS).abs() <= tol),
        "Vital status missing": ev_all.VITAL_STATUS.isna(),
    }
    for k, m in cat.items():
        comp.append({"category": k, "n": int(m.sum()), "fraction_of_964_usable_events": float(m.sum() / n_ev)})
    comp = pd.DataFrame(comp)
    comp.loc[len(comp)] = {"category": "TOTAL usable events", "n": n_ev, "fraction_of_964_usable_events": 1.0}
    save_table(comp, "p4_metabric_event_composition")
    print(comp.round(4).to_string())
    # bounds on the share of events that could be 'death first' (cannot be separated from recurrence)
    dd = ev_all[(ev_all.VITAL_STATUS == "Died of Disease")]
    n_first_possible = int(((dd.OS_MONTHS - dd.RFS_MONTHS).abs() <= tol).sum())
    n_unknown = int(ev_all.VITAL_STATUS.isna().sum())
    ub = float(n_first_possible / n_ev)
    ub_incl_unknown = float((n_first_possible + n_unknown) / n_ev)
    # post-relapse interval among those who died of disease with a recorded earlier RFS (descriptive; conditional on death)
    gap = (dd.OS_MONTHS - dd.RFS_MONTHS)[(dd.OS_MONTHS - dd.RFS_MONTHS) > tol]
    save_table(pd.DataFrame([{"quantity": "upper bound on fraction of events that are death-first, among events with KNOWN vital status (RFS == OS, died of disease)", "value": ub},
                             {"quantity": "upper bound if ALL events with missing vital status were also death-first", "value": ub_incl_unknown},
                             {"quantity": "events with documented relapse before death or no death (Recurred & Living, Recurred & other-cause death, Died of Disease with RFS < OS)", "value": float((92 + 65 + 545) / n_ev) if False else float((ev_all.VITAL_STATUS == "Living").sum() + (ev_all.VITAL_STATUS == "Died of Other Causes").sum() + ((ev_all.VITAL_STATUS == "Died of Disease") & ((ev_all.OS_MONTHS - ev_all.RFS_MONTHS) > tol)).sum()) / n_ev},
                             {"quantity": "lower bound on that fraction", "value": 0.0},
                             {"quantity": "n died of disease with RFS < OS", "value": len(gap)},
                             {"quantity": "OS - RFS among them: median / IQR (months), CONDITIONAL on death (biased; not a delay distribution)",
                              "value": f"{gap.median():.1f} / {gap.quantile(0.25):.1f}-{gap.quantile(0.75):.1f}"},
                             {"quantity": "note", "value": "OS-RFS is the time from the recorded RFS event to death only if both clocks share an origin; descriptive"}]), "p4_metabric_death_vs_recurrence_bounds")

    # ------------------------------------------------------------------ 3. early-time structure (origin plausibility; INFERENCE only)
    early = []
    for thr in (0.5, 1, 3, 6, 12):
        sub = both[both.RFS_MONTHS < thr]
        early.append({"RFS_MONTHS_below": thr, "n_records": len(sub), "n_events": int(sub.rfs_event.sum())})
    early.append({"RFS_MONTHS_below": 0.0001, "n_records": int((both.RFS_MONTHS == 0).sum()), "n_events": int(both[both.RFS_MONTHS == 0].rfs_event.sum())})
    save_table(pd.DataFrame(early), "p4_metabric_early_times")

    # ------------------------------------------------------------------ 4. treatment timing / what is recorded
    avail = []
    for col in ("CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY", "BREAST_SURGERY", "AGE_AT_DIAGNOSIS"):
        avail.append({"field": col, "n_nonmissing": int(pat[col].notna().sum()), "type": "binary flag / category (no date, no duration)" if col != "AGE_AT_DIAGNOSIS" else "age at diagnosis (years)"})
    save_table(pd.DataFrame(avail), "p4_metabric_treatment_fields")

    # ------------------------------------------------------------------ 5. subgroup reliability (events / at-risk criteria)
    km = pd.read_csv(TAB / "metabric_km_subgroups.csv")
    cols = km.columns.tolist()
    save_table(km, "p4_metabric_subgroup_table_phase1_copy")
    # re-derive at-risk by subgroup from raw data for the criterion
    both = both.merge(smp[["PATIENT_ID", "ER_STATUS", "HER2_STATUS", "GRADE", "TUMOR_STAGE", "TUMOR_SIZE"]], on="PATIENT_ID", how="left")
    groups = {"ER (IHC) positive": both.ER_IHC == "Positve", "ER (IHC) negative": both.ER_IHC == "Negative",
              "PAM50: LumA": both.CLAUDIN_SUBTYPE == "LumA", "PAM50: LumB": both.CLAUDIN_SUBTYPE == "LumB", "PAM50: Her2": both.CLAUDIN_SUBTYPE == "Her2",
              "PAM50: Basal": both.CLAUDIN_SUBTYPE == "Basal", "PAM50: claudin-low": both.CLAUDIN_SUBTYPE == "claudin-low", "PAM50: Normal": both.CLAUDIN_SUBTYPE == "Normal",
              "all usable": both.rfs_event.notna()}
    rel = []
    for g, m in groups.items():
        sub = both[m]
        row = {"group": g, "n": len(sub), "events": int(sub.rfs_event.sum())}
        for h in (60, 120, 180, 240):
            row[f"at_risk_ge_{h}mo"] = int((sub.RFS_MONTHS >= h).sum())
            row[f"events_after_{h}mo"] = int(((sub.RFS_MONTHS > h) & (sub.rfs_event == 1)).sum())
        row["reliable_5y (>=50 events, >=100 at risk at 60 mo)"] = bool(row["events"] >= 50 and row["at_risk_ge_60mo"] >= 100)
        row["reliable_10y (>=50 events, >=100 at risk at 120 mo)"] = bool(row["events"] >= 50 and row["at_risk_ge_120mo"] >= 100)
        row["reliable_20y (>=50 events, >=30 at risk at 240 mo)"] = bool(row["events"] >= 50 and row["at_risk_ge_240mo"] >= 30)
        rel.append(row)
    rel = pd.DataFrame(rel)
    save_table(rel, "p4_metabric_subgroup_reliability")
    print(rel[["group", "n", "events", "at_risk_ge_60mo", "at_risk_ge_120mo", "at_risk_ge_240mo", "reliable_5y (>=50 events, >=100 at risk at 60 mo)", "reliable_10y (>=50 events, >=100 at risk at 120 mo)", "reliable_20y (>=50 events, >=30 at risk at 240 mo)"]].to_string())
    missing = {c: int(pat[c].isna().sum()) for c in ("ER_IHC", "CLAUDIN_SUBTYPE", "INTCLUST", "COHORT", "CHEMOTHERAPY")}
    save_table(pd.DataFrame([{"variable": k, "n_missing_of_2509": v} for k, v in missing.items()]), "p4_metabric_subgroup_variable_missingness")

    # ------------------------------------------------------------------ 6. documentation facts (verified vs not)
    facts = [
        ("RFS_STATUS definition", "VERIFIED (export metadata, verbatim)", "Relapse Free Status (Recurred = cases where there was loco-regional relapse, distant relapse or death-specific death)", "cBioPortal meta_clinical_patient / data_clinical_patient header"),
        ("RFS_MONTHS definition", "VERIFIED (export metadata, verbatim)", "Time until last follow-up or relapse (in months); NO time origin stated", "same"),
        ("OS_MONTHS definition", "VERIFIED (export metadata, verbatim)", "Overall survival in months since initial diagonosis", "same"),
        ("Composite event includes disease-specific death", "VERIFIED (export metadata)", "yes; all 646 'Died of Disease' are Recurred", "export + Phase 1 audit"),
        ("Loco-regional vs distant relapse distinguishable", "VERIFIED NOT POSSIBLE", "not in export", "export"),
        ("Disease-specific death distinguishable from recurrence", "PARTLY (bounds only)", "see p4_metabric_event_composition: death-first cannot be separated when RFS == OS", "this analysis"),
        ("RFS origin = diagnosis", "NOT VERIFIED", "no documentation in the export, the datahub commit history (Dec 2020 'Added RFS data', 'Update RFS STATUS description'; Apr 2024 'update survival data'; no source named) or the three METABRIC papers' open text", "documentation search 2026-11-04 (Europe PMC PMC3440846, PMC4866047, PMC6647838; GitHub API)"),
        ("Rueda 2019 uses surgery as the start of the disease-free state", "VERIFIED (paper text)", "'transition probabilities from surgery ... disease-free/after surgery'; says nothing about the cBioPortal RFS_MONTHS field", "PMC6647838"),
        ("Common origin of OS and RFS clocks", "INFERENCE (strong, not proof)", "see p4_metabric_origin_tests: equality for censored-alive patients", "this analysis"),
        ("Treatment dates / durations available", "VERIFIED NOT AVAILABLE", "only binary flags; no dates", "export"),
        ("Diagnosis-to-treatment interval", "NOT AVAILABLE", "unknown; remains a nuisance parameter tau", "-"),
    ]
    save_table(pd.DataFrame(facts, columns=["item", "status", "content", "source"]), "p4_metabric_documentation_facts")
    log("done")


if __name__ == "__main__":
    main()
