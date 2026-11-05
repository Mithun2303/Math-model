"""Phase 1 alternative-dataset audit: TCGA-BRCA PanCancer Atlas clinical export (DFS/PFS/OS/DSS fields)."""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
path = ROOT / "data" / "raw" / "brca_tcga_pan_can_atlas_2018" / "data_clinical_patient.txt"
meta = [l.rstrip("\n").split("\t") for _, l in zip(range(4), open(path))]
df = pd.read_csv(path, sep="\t", skiprows=4, dtype=str, keep_default_na=False).replace({"": np.nan, "[Not Available]": np.nan})
print("patients:", len(df), "unique:", df.PATIENT_ID.nunique())
for c in ["DFS_STATUS", "PFS_STATUS", "OS_STATUS", "DSS_STATUS", "SUBTYPE", "NEOADJUVANT" if False else "HISTORY_NEOADJUVANT_TRTYN",
          "NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT", "PERSON_NEOPLASM_CANCER_STATUS"]:
    print(f"\n{c}\n{df[c].value_counts(dropna=False).to_string()}")
for c in ["DFS", "PFS", "OS", "DSS"]:
    i = list(df.columns).index(f"{c}_STATUS")
    print(f"\n{c}_STATUS description: {meta[1][i]!r} | {c}_MONTHS description: {meta[1][i+1]!r}")
    t = pd.to_numeric(df[f"{c}_MONTHS"], errors="coerce")
    e = df[f"{c}_STATUS"].str.split(":").str[0].map({"0": 0, "1": 1})
    ok = t.notna() & e.notna()
    print(f"  usable n={int(ok.sum())} events={int(e[ok].sum())} median time={t[ok].median():.1f} max={t[ok].max():.1f} "
          f"at-risk>=60mo={int((t[ok] >= 60).sum())} >=120mo={int((t[ok] >= 120).sum())} >=240mo={int((t[ok] >= 240).sum())} "
          f"events after 60mo={int(((t > 60) & (e == 1)).sum())} after 120mo={int(((t > 120) & (e == 1)).sum())}")
