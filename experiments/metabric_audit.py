"""Phase 1 data audit of the cBioPortal METABRIC clinical export (brca_metabric).

Reads data/raw/brca_metabric/data_clinical_{patient,sample}.txt, prints/saves every statistic quoted in
data/dataset_feasibility.md and writes data/clinical_data_dictionary.csv. Nothing here is imputed.
"""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "brca_metabric"
OUT = ROOT / "results" / "tables"
OUT.mkdir(parents=True, exist_ok=True)


def read_cbio(path: Path):
    """cBioPortal clinical file: 4 '#' metadata rows (label, description, type, priority) then header."""
    meta = [l.rstrip("\n").split("\t") for _, l in zip(range(4), open(path))]
    df = pd.read_csv(path, sep="\t", skiprows=4, dtype=str, keep_default_na=False)
    df = df.replace({"": np.nan, "NA": np.nan, "NaN": np.nan})
    info = {c: dict(label=meta[0][i].lstrip("#"), description=meta[1][i].lstrip("#"), declared_type=meta[2][i].lstrip("#"))
            for i, c in enumerate(df.columns)}
    return df, info


pat, pinfo = read_cbio(RAW / "data_clinical_patient.txt")
smp, sinfo = read_cbio(RAW / "data_clinical_sample.txt")
print(f"patient rows={len(pat)} unique ids={pat.PATIENT_ID.nunique()}; sample rows={len(smp)} "
      f"unique patients={smp.PATIENT_ID.nunique()} unique samples={smp.SAMPLE_ID.nunique()}")
df = pat.merge(smp, on="PATIENT_ID", how="outer", indicator=True)
print("merge:", df["_merge"].value_counts().to_dict())
for c in ["RFS_MONTHS", "OS_MONTHS", "AGE_AT_DIAGNOSIS", "NPI", "TUMOR_SIZE", "TMB_NONSYNONYMOUS",
          "LYMPH_NODES_EXAMINED_POSITIVE"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# ---------- value inventories of categorical fields ----------
print("\n=== categorical value counts (incl. NaN) ===")
for c in ["RFS_STATUS", "OS_STATUS", "VITAL_STATUS", "COHORT", "ER_IHC", "ER_STATUS", "HER2_SNP6", "HER2_STATUS",
          "PR_STATUS", "CLAUDIN_SUBTYPE", "INTCLUST", "THREEGENE", "CHEMOTHERAPY", "HORMONE_THERAPY",
          "RADIO_THERAPY", "BREAST_SURGERY", "SAMPLE_TYPE", "GRADE", "TUMOR_STAGE", "SEX",
          "INFERRED_MENOPAUSAL_STATE", "CANCER_TYPE_DETAILED"]:
    print(f"\n{c}\n{df[c].value_counts(dropna=False).to_string()}")

# ---------- endpoint semantics ----------
print("\n=== endpoint semantics ===")
df["rfs_event"] = df.RFS_STATUS.str.split(":").str[0].map({"0": 0, "1": 1})
df["os_event"] = df.OS_STATUS.str.split(":").str[0].map({"0": 0, "1": 1})
print("RFS status x vital status:\n", pd.crosstab(df.RFS_STATUS, df.VITAL_STATUS, dropna=False).to_string())
print("\nRFS status x OS status:\n", pd.crosstab(df.RFS_STATUS, df.OS_STATUS, dropna=False).to_string())
ok = df.RFS_MONTHS.notna() & df.OS_MONTHS.notna()
print("\nRFS_MONTHS > OS_MONTHS (tol 1e-6):", int((df.RFS_MONTHS[ok] > df.OS_MONTHS[ok] + 1e-6).sum()),
      "| RFS == OS:", int((np.abs(df.RFS_MONTHS[ok] - df.OS_MONTHS[ok]) < 1e-6).sum()), "of", int(ok.sum()))
for status in df.VITAL_STATUS.dropna().unique():
    m = (df.VITAL_STATUS == status) & ok
    same = (np.abs(df.RFS_MONTHS[m] - df.OS_MONTHS[m]) < 1e-6).mean() if m.sum() else np.nan
    ev = df.rfs_event[m].mean() if m.sum() else np.nan
    print(f"  {status:28s} n={int(m.sum()):5d}  frac RFS_time==OS_time={same:.3f}  frac RFS event={ev:.3f}")

# ---------- time variable ----------
d = df.dropna(subset=["RFS_MONTHS", "rfs_event"]).copy()
print(f"\nusable RFS records: {len(d)} of {len(df)}  (missing time={int(df.RFS_MONTHS.isna().sum())}, "
      f"missing status={int(df.rfs_event.isna().sum())})")
print("RFS_MONTHS describe:\n", d.RFS_MONTHS.describe(percentiles=[.05, .25, .5, .75, .95]).round(2).to_string())
print("zero or negative times:", int((d.RFS_MONTHS <= 0).sum()))
print("events:", int(d.rfs_event.sum()), "censored:", int((d.rfs_event == 0).sum()))

# reverse Kaplan-Meier median follow-up
def km(time, event):
    order = np.argsort(time)
    t, e = np.asarray(time)[order], np.asarray(event)[order]
    ut = np.unique(t[e == 1])
    S, var_sum, out = 1.0, 0.0, []
    for u in ut:
        n = (t >= u).sum()
        dd = ((t == u) & (e == 1)).sum()
        S *= 1 - dd / n
        var_sum += dd / (n * (n - dd)) if n > dd else 0.0
        out.append((u, S, S * np.sqrt(var_sum), n))
    return pd.DataFrame(out, columns=["t", "S", "se", "n_at_risk"])

def km_at(curve, tq):
    c = curve[curve.t <= tq]
    return (1.0, 0.0) if c.empty else (c.S.iloc[-1], c.se.iloc[-1])

rev = km(d.RFS_MONTHS.values, 1 - d.rfs_event.values)
med_fu = rev.t[rev.S <= 0.5].min()
print(f"median potential follow-up (reverse KM) = {med_fu:.1f} months")

print("\nevents by window (months since time origin):")
bins = [0, 60, 120, 180, 240, 1e9]
labels = ["0-60", "60-120", "120-180", "180-240", ">240"]
d["win"] = pd.cut(d.RFS_MONTHS, bins=bins, labels=labels, right=True)
win = d.groupby("win", observed=False).agg(n_records=("rfs_event", "size"), events=("rfs_event", "sum"),
                                           censored=("rfs_event", lambda x: int((x == 0).sum())))
print(win.to_string())
win.to_csv(OUT / "metabric_events_by_window.csv")
print("\nnumber still at risk (RFS_MONTHS >= t):",
      {t: int((d.RFS_MONTHS >= t).sum()) for t in (60, 120, 180, 240)})

curve = km(d.RFS_MONTHS.values, d.rfs_event.values)
rows = []
for yrs in (5, 10, 15, 20):
    s, se = km_at(curve, 12 * yrs)
    nrisk = int((d.RFS_MONTHS >= 12 * yrs).sum())
    rows.append(dict(horizon_years=yrs, n_at_risk=nrisk, km_event_free=round(s, 4), se=round(se, 4),
                     ci95_low=round(max(0, s - 1.96 * se), 4), ci95_high=round(min(1, s + 1.96 * se), 4),
                     cumulative_incidence_1_minus_S=round(1 - s, 4)))
kmtab = pd.DataFrame(rows)
print("\nOverall KM (RFS, as defined by export):\n", kmtab.to_string(index=False))
kmtab.to_csv(OUT / "metabric_km_overall.csv", index=False)

# ---------- subgroup KM ----------
def subgroup(col, minimum_events=20):
    out = []
    for g, sub in d.groupby(col, dropna=False):
        ev, n = int(sub.rfs_event.sum()), len(sub)
        row = dict(variable=col, group=("(missing)" if pd.isna(g) else g), n=n, events=ev,
                   usable=ev >= minimum_events)
        c = km(sub.RFS_MONTHS.values, sub.rfs_event.values) if ev else None
        for yrs in (5, 10, 20):
            if c is None:
                row[f"S{yrs}y"] = np.nan
            else:
                row[f"S{yrs}y"] = round(km_at(c, 12 * yrs)[0], 3)
            row[f"atrisk{yrs}y"] = int((sub.RFS_MONTHS >= 12 * yrs).sum())
        out.append(row)
    return pd.DataFrame(out)

sg = pd.concat([subgroup(c) for c in ["ER_IHC", "CLAUDIN_SUBTYPE", "INTCLUST", "THREEGENE", "HER2_SNP6", "COHORT",
                                      "GRADE", "CHEMOTHERAPY", "HORMONE_THERAPY"]])
sg.to_csv(OUT / "metabric_km_subgroups.csv", index=False)
print("\nSubgroup KM (event-free fraction at 5/10/20 y; atrisk = n with follow-up >= horizon):")
print(sg.to_string(index=False))

# ---------- late recurrence feasibility: patients recurrence-free beyond 5 y ----------
lm = d[d.RFS_MONTHS > 60]
print(f"\nLandmark at 60 mo: {len(lm)} patients event-free & under follow-up at 5 y; "
      f"subsequent events={int(lm.rfs_event.sum())}")
er = lm[(lm.ER_IHC == "Positve")]
print(f"  ER-positive (IHC) landmark cohort: n={len(er)}, subsequent events={int(er.rfs_event.sum())}")

# ---------- missingness / data dictionary ----------
all_info = {**pinfo, **sinfo}
n = len(df)
dict_rows = []
for c in list(pat.columns) + [c for c in smp.columns if c != "PATIENT_ID"]:
    s = df[c]
    numeric = c in ("RFS_MONTHS", "OS_MONTHS", "AGE_AT_DIAGNOSIS", "NPI", "TUMOR_SIZE", "TMB_NONSYNONYMOUS",
                    "LYMPH_NODES_EXAMINED_POSITIVE")
    dict_rows.append(dict(
        field=c, source_file=("data_clinical_patient.txt" if c in pat.columns else "data_clinical_sample.txt"),
        cbioportal_label=all_info[c]["label"], declared_type=all_info[c]["declared_type"],
        cbioportal_description=all_info[c]["description"],
        n_non_missing=int(s.notna().sum()), n_missing=int(s.isna().sum()), pct_missing=round(100 * s.isna().mean(), 2),
        n_unique=int(s.nunique(dropna=True)),
        observed_values=("numeric: min=%.4g max=%.4g" % (s.min(), s.max())) if numeric and s.notna().any()
        else "identifier (not listed)" if c in ("PATIENT_ID", "SAMPLE_ID") else "; ".join(f"{k} ({v})" for k, v in s.value_counts().head(12).items()),
        units=("months" if c in ("RFS_MONTHS", "OS_MONTHS") else "years" if c == "AGE_AT_DIAGNOSIS" else
               "mm (not verified)" if c == "TUMOR_SIZE" else "count" if c == "LYMPH_NODES_EXAMINED_POSITIVE" else
               "not applicable/not verified"),
    ))
dd = pd.DataFrame(dict_rows)
semantic = {
    "RFS_STATUS": "EVENT INDICATOR. '1:Recurred' = loco-regional relapse, distant relapse OR disease-specific death (export description); "
                  "'0:Not Recurred' = censored. Single composite event: loco-regional vs distant NOT separable. No date of relapse.",
    "RFS_MONTHS": "Time to first RFS event or last follow-up (months). Time origin not stated for RFS in export (OS is stated as 'since initial diagnosis'); "
                  "inferred as diagnosis - NOT treatment completion. See feasibility doc.",
    "OS_STATUS": "Overall survival event (death from any cause); 1:DECEASED.",
    "OS_MONTHS": "Months since initial diagnosis (stated in export) to death or last follow-up.",
    "VITAL_STATUS": "Cause-specific vital status (Living / Died of Disease / Died of Other Causes). LEAKS outcome; never use as a predictor.",
    "SAMPLE_TYPE": "Specimen type of the sequenced tumour (primary); not a recurrence indicator.",
}
dd["event_semantics_or_caveat"] = dd.field.map(semantic).fillna("")
leak = {"RFS_STATUS", "RFS_MONTHS", "OS_STATUS", "OS_MONTHS", "VITAL_STATUS"}
post = {"CHEMOTHERAPY", "HORMONE_THERAPY", "RADIO_THERAPY", "BREAST_SURGERY"}
dd["leakage_risk"] = dd.field.map(lambda f: "OUTCOME - exclude from predictors" if f in leak else
                                  "treatment (timing/assignment unverified; confounded by indication)" if f in post else
                                  "baseline at diagnosis (presumed; not date-verified)")
dd.to_csv(ROOT / "data" / "clinical_data_dictionary.csv", index=False)
print("\nwrote clinical_data_dictionary.csv with", len(dd), "fields")
print(dd[["field", "n_non_missing", "n_missing", "pct_missing"]].to_string(index=False))
d.to_csv(OUT / "metabric_rfs_analysis_cohort_summary_NOT_PATIENT_LEVEL.csv", index=False) if False else None
