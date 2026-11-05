# Dataset Feasibility: METABRIC (cBioPortal `brca_metabric`) for Recurrence-Time Validation

*Phase 1 deliverable. Every number below was computed from the downloaded files by `experiments/metabric_audit.py` and `experiments/tcga_audit.py` (logs: `results/tables/metabric_audit_log.txt`, `results/tables/tcga_audit_log.txt`). Nothing is imputed.*

## 1. Provenance

| Item | Value |
|---|---|
| Portal | https://www.cbioportal.org/study?id=brca_metabric (portal page itself not scraped; files taken from the repository) |
| Repository | https://github.com/cBioPortal/datahub/tree/master/public/brca_metabric |
| Download date | 2026-11-03 (UTC) |
| Latest repo commit touching the folder at that time | `86690e1ed9`, 2026-01-30 ("Gene Table v7: migrate outdated gene symbols") |
| Files used | `data_clinical_patient.txt` (SHA-256 `2a6bbadc…6e4b2b07c…d6454995`, see `report/references.bib`), `data_clinical_sample.txt`, `meta_*.txt` |
| Storage | Large files are Git LFS objects; the plain `raw.githubusercontent.com` URL returns a 131-byte pointer. Real content was obtained from `media.githubusercontent.com`. |
| Licence | Open Database License (ODbL 1.0) per repository `LICENSE`. Share-alike and attribution apply to derived datasets. |
| Study metadata | "Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)"; "Targeted sequencing of 2509 primary breast tumors with 548 matched normals"; citations: Pereira 2016, Rueda 2019, Curtis 2012. |
| Raw files kept in project | `data/raw/brca_metabric/` and `data/raw/brca_tcga_pan_can_atlas_2018/` (not committed anywhere; there is no git repository) |

## 2. Files actually available

The repository folder contains clinical patient and sample tables plus molecular matrices (CNA, expression z-scores, mutations, methylation, gene-panel matrix; not downloaded in Phase 1) and metadata. Only the two clinical files are needed for the endpoint audit.

- `data_clinical_patient.txt`: 2,509 patients × 24 columns.
- `data_clinical_sample.txt`: 2,509 samples × 13 columns (one primary sample per patient; `SAMPLE_TYPE` = Primary for all 2,509).
- Patient and sample tables join 1:1 on `PATIENT_ID` (2,509 matched, 0 unmatched).

Field-by-field documentation, units, missingness and leakage classification: `data/clinical_data_dictionary.csv` (36 fields).

## 3. Verified endpoint fields

**The export contains `RFS_STATUS` and `RFS_MONTHS`. It does not contain `DFS_*`, `PFS_*` or any recurrence date.**

| Field | Export description (verbatim) | Observed values |
|---|---|---|
| `RFS_STATUS` | "Relapse Free Status (Recurred = cases where there was loco-regional relapse, distant relapse or death-specific death)" | `0:Not Recurred` 1,486; `1:Recurred` 1,002; missing 21 |
| `RFS_MONTHS` | "Time until last follow-up or relapse (in months)" | 2,388 non-missing; 121 missing |
| `OS_STATUS` / `OS_MONTHS` | "Overall patient survival status." / "Overall survival in months since initial diagonosis." | 1,981 non-missing; 528 missing |
| `VITAL_STATUS` | "The survival state of the person." | Living 837; Died of Disease 646; Died of Other Causes 497; missing 529 |

### Interpretation (what is and is not established)

1. **Event semantics.** `Recurred` is a **composite**: loco-regional relapse, distant relapse, or disease-specific death (the export's text says "death-specific death"; we read this as disease/breast-cancer-specific death). Loco-regional and distant events **cannot be separated** in this export. This differs from the four-stage multistate framing in `Rueda2019`, so results here must not be described as "distant recurrence".
2. **Competing risks are censored, not counted.** Of 497 patients who died of other causes, 432 are `Not Recurred` (censored) and 65 `Recurred`. Cross-tabulation: all 646 "Died of Disease" patients are `Recurred`. So this endpoint is *recurrence-free survival with non-disease deaths censored*. It is **not** "disease-free survival including death from any cause" and **not** overall survival. This is consistent with the concern raised by `Tolaney2021` (ABS) that endpoint definitions differ in whether non-breast deaths and second primaries are included. Censoring at other-cause death is informative for cumulative incidence: 1 − KM overstates the incidence of recurrence relative to a competing-risk (Aalen–Johansen) estimate.
3. **Time origin.** The export states OS time is "since initial diagnosis" but gives **no origin for RFS**. Evidence for the same origin: RFS_MONTHS ≤ OS_MONTHS in all 1,981 patients with both (0 violations); RFS equals OS exactly for 1,282 of them. This is consistent with a diagnosis-based origin but does **not prove it**. The time origin is certainly *not* the end of treatment, and **no treatment dates are provided**, so the model's t = 0 (post-treatment residual disease) cannot be aligned with the data's t = 0 beyond an unknown offset.
4. **Units.** Months, non-integer (e.g. 84.6333…). Max RFS_MONTHS = 389.3 (about 32 years).
5. **Censoring.** Administrative/loss-to-follow-up censoring is not distinguished from death censoring. A patient with `Not Recurred` and `Living` is the only unambiguous "alive and recurrence-free at last contact" case (745 patients).
6. **Data quirks.** 4 patients have RFS_MONTHS = 0 (3 of them events). 13 patients are missing both time and status, 108 only time, 8 only status.

## 4. Verified usable cohort and counts

| Quantity | Value |
|---|---|
| Patients in export | 2,509 (female 2,509) |
| Both RFS time and status present (**usable for survival analysis**) | **2,380** |
| Excluded for missing RFS information | 129 (13 missing both + 108 missing time only + 8 missing status only) |
| Events (`Recurred`) | **964** |
| Censored | 1,416 |
| Median potential follow-up (reverse Kaplan–Meier) | 157.7 months |
| RFS_MONTHS median / IQR / 95th percentile | 100.1 / 41.0–169.4 / 253.0 months |
| Usable patients who also have the extended clinical/molecular fields | 1,979 (401 usable patients lack them) |

The 1,980 patients with `INTCLUST`, PAM50, therapy fields etc. coincides with the "1,980 for whom molecular data are available" in the Rueda et al. 2019 abstract. This is a numerical agreement only, not proof that the files are the same dataset.

### Events and follow-up by window (months)

| Window | Records with time in window | Events | Censored |
|---|---|---|---|
| 0–60 | 795 | 601 | 194 |
| 60–120 | 597 | 215 | 382 |
| 120–180 | 464 | 104 | 360 |
| 180–240 | 358 | 32 | 326 |
| >240 | 162 | 9 | 153 |

(Three events at time 0 are not in these windows, since the first bin excludes 0.) Number still under follow-up at ≥ 60/120/180/240 months: 1,581 / 984 / 521 / 162.

### Overall Kaplan–Meier (RFS as defined; own implementation, Greenwood standard error)

| Horizon | At risk | Event-free | 95% CI | 1 − S |
|---|---|---|---|---|
| 5 y | 1,581 | 0.737 | 0.719–0.755 | 0.263 |
| 10 y | 984 | 0.624 | 0.604–0.645 | 0.376 |
| 15 y | 521 | 0.544 | 0.521–0.567 | 0.456 |
| 20 y | 162 | 0.495 | 0.468–0.522 | 0.505 |

The 20-year estimate rests on 162 patients and 9 events beyond 240 months, and the follow-up in that tail comes mainly from two sub-cohorts (cohort 2: 68, cohort 3: 60 patients at risk at 20 y). Treat it as low-precision. A CI that ignores between-cohort heterogeneity understates uncertainty.

### Subgroups (event-free fraction; counts of patients still under follow-up at each horizon in `results/tables/metabric_km_subgroups.csv`)

| Group | n | Events | S(5 y) | S(10 y) | S(20 y) | At risk at 20 y |
|---|---|---|---|---|---|---|
| ER (IHC) positive | 1,735 | 689 | 0.775 | 0.642 | 0.485 | 110 |
| ER (IHC) negative | 573 | 258 | 0.610 | 0.558 | 0.503 | 49 |
| PAM50+claudin-low: LumA | 699 | 231 | 0.844 | 0.733 | 0.539 | 61 |
| LumB | 475 | 227 | 0.697 | 0.539 | 0.421 | 27 |
| Her2 | 224 | 116 | 0.609 | 0.490 | 0.426 | 18 |
| Basal | 209 | 93 | 0.591 | 0.563 | 0.530 | 13 |
| IntClust 3 / 7 / 8 / 10 (shown only to illustrate size and precision) | 290 / 189 / 299 / 226 | 87 / 56 / 115 / 86 | 0.873 / 0.862 / 0.782 / 0.652 | 0.788 / 0.739 / 0.665 / 0.641 | 0.532 / 0.620 / 0.530 / 0.598 | 22 / 22 / 23 / 19 |

(The IntClust rows are arbitrary examples. They are **not** claimed to be the four late-recurring clusters of `Rueda2019`; the abstract does not name them, so identifying those clusters requires the full text.) Any subgroup curve at 20 years has ≤ 110 patients at risk; for several (e.g. Basal, 13) it is not interpretable.

**Landmark cohort for late relapse.** Patients with RFS_MONTHS > 60: 1,581, with 360 subsequent events. ER-positive (IHC) among these: 1,225 with 312 subsequent events. This is the most direct slice for studying relapse after the typical treatment period.

### Missingness of key baseline variables (of 2,509)

| Variable | Missing | Variable | Missing |
|---|---|---|---|
| ER_IHC | 83 (3.3%) | GRADE | 121 (4.8%) |
| AGE_AT_DIAGNOSIS | 11 (0.4%) | TUMOR_SIZE | 149 (5.9%) |
| LYMPH_NODES_EXAMINED_POSITIVE | 266 (10.6%) | TUMOR_STAGE | 721 (28.7%) |
| NPI | 222 (8.8%) | INTCLUST / CLAUDIN_SUBTYPE / CHEMOTHERAPY / HORMONE_THERAPY / RADIO_THERAPY | 529 (21.1%) each |
| THREEGENE | 745 (29.7%) | BREAST_SURGERY | 554 (22.1%) |

## 5. Potential data leakage and bias

| Issue | Detail | Handling |
|---|---|---|
| Outcome-derived columns | `VITAL_STATUS` (includes "Died of Disease"), `OS_*`, `RFS_*` | Never use as predictors |
| Treatment variables | `CHEMOTHERAPY`, `HORMONE_THERAPY`, `RADIO_THERAPY`, `BREAST_SURGERY`: assignment depends on prognosis (confounding by indication), timing unverified | Exclude from prediction baselines or treat as strata with explicit caveat |
| Derived baseline index | `NPI` is derived from size, nodes and grade | Do not combine with its components in one model without noting redundancy |
| Missing-not-at-random | 529 patients lack extended fields; 21 lack RFS status | Report complete-case size; compare RFS of those with vs without extended fields (their KM at 5/10/20 y is in `metabric_km_subgroups.csv` as "(missing)") |
| Cohort heterogeneity | `COHORT` 1–9 have very different long-term follow-up (at risk at 20 y: cohort 1: 16, cohort 2: 68, cohort 3: 60) | Stratify or cluster by cohort |
| Selection / era | Primary tumours treated in earlier eras; treatment practice differs from today | State limitations |
| Competing risk | Non-disease deaths are censored | Use Aalen–Johansen where cause is available (1,980 patients) |

## 6. Suitability for validating simulated relapse-time distributions

| Requirement | Met? | Comment |
|---|---|---|
| Event indicator documented | **Yes** | `RFS_STATUS`, composite definition stated |
| Event/censoring time | **Yes (2,380 usable)** | months |
| Defensible time origin | **Partly** | Diagnosis, inferred; unrelated to treatment end; no dates |
| Enough events and long follow-up | **Yes for ≤ 10 y; weak at 15–20 y** | 964 events; 162 at risk beyond 20 y |
| Event type matches model (threshold-crossing of tumour burden) | **Partly** | Composite incl. loco-regional relapse and cancer death |
| Access permission | **Yes** | ODbL, public |
| Measures P, D or E | **No** | Never; only the event time |

## 7. Alternative datasets considered

**TCGA-BRCA PanCancer Atlas (`brca_tcga_pan_can_atlas_2018`), downloaded and audited.**

| Item | Value |
|---|---|
| Patients | 1,084 |
| `DFS_STATUS`: DiseaseFree / Recurred-Progressed / missing | 858 / 84 / 142 |
| DFS definition (export) | "Disease free status since initial treatment." (origin = initial treatment) |
| Usable DFS records | 941; 84 events; median time 25.0 months; max 281.3 |
| Under follow-up ≥ 60 / 120 / 240 months | 190 / 29 / 4 |
| DFS events after 60 / 120 months | 10 / 1 |
| `PFS` events (usable 1,082) | 145; 22 after 60 months |

**Verdict:** TCGA offers a better-defined origin (initial treatment) but its follow-up is too short and events too few for long-latency relapse. It is **not suitable for 10- or 20-year validation**; it could serve only to *cross-check endpoint definitions*. TCGA-CDR (curated endpoints) and original METABRIC data from the European Genome-phenome Archive were **not examined** and remain candidates. Datasets are **not combined**: METABRIC (origin = diagnosis, RFS composite) and TCGA (origin = initial treatment, DFS) have non-reconcilable endpoints and origins.

## 8. Decision

### **CONDITIONAL GO: group-level, exploratory, ≤ 10-year (cautiously 15-year) comparison only. NO-GO for patient-level prediction and for any claim about dormant-cell parameters.**

**Justification.** METABRIC provides a documented relapse-type event (`RFS_STATUS`), 2,380 usable records, 964 events and enough follow-up to constrain 5- and 10-year event-free probabilities and the post-5-year hazard. It cannot validate individual-level or mechanistic claims, and the endpoint differs from the model's definition in three ways (composite event, diagnosis origin, competing-risk censoring).

**Conditions for the clinical track to proceed:**

1. Always label the target "RFS as defined in the cBioPortal METABRIC export", never "recurrence", "DFS" or "OS".
2. Introduce an explicit offset τ between diagnosis and the start of the model's post-treatment state (not recorded; to be varied as a nuisance, or fixed with justification).
3. Use 5- and 10-year horizons as primary targets. Report 15 and 20 years with their at-risk counts and as exploratory.
4. Stratify or cluster by `COHORT` when estimating uncertainty; report complete-case sizes.
5. Treat any fit between simulated and empirical curves as **consistency**, not validation. Parameters remain weakly identified.
6. If conditions cannot be met (for example, if the offset τ dominates the shape), retain the project as an in-silico study and record the clinical comparison as a documented limitation.

## 9. Analyses feasible with the verified data

| Analysis | Feasible? | Notes |
|---|---|---|
| Kaplan–Meier overall and by ER/PAM50/IntClust/grade/cohort | **Yes** | Done in Phase 1 (tables in `results/tables/`) |
| Right-censoring-aware 5/10/20-year probabilities with CIs | **Yes** | 20-year low precision |
| Landmark analysis after 5 years (ER+) | **Yes** | 1,225 patients, 312 subsequent events |
| Competing-risk cumulative incidence (Aalen–Johansen) | **Partly** | Needs cause: `VITAL_STATUS`, 1,980 patients |
| Cox proportional-hazards baseline | **Yes** | Use baseline covariates only: age, NPI/size/nodes/grade, ER, IntClust/PAM50. Not tuned on treatment or outcome columns. Package install needed. |
| Concordance index, time-dependent AUC, IPCW Brier score, calibration | **Yes** | Patient-level train/test split with cohort stratification |
| Separate loco-regional vs distant relapse | **No** | Not in export |
| Absolute time since treatment end | **No** | No treatment dates |
| Estimate dormant-cell counts or transition rates from patients | **No** | Never measured |
| Direct comparison of model event time with an individual's relapse | **No** | Model is deterministic per parameter set; data are group-level |

## 10. Distinction to preserve in all later reports

A **simulation-based event probability** F̂(t) is the fraction of simulated scenarios, under chosen parameter distributions, whose modelled burden crosses N_detect by time t. A **clinically validated patient risk** is a calibrated, externally validated prediction for real patients. This project produces only the former. Kaplan–Meier curves from METABRIC describe cohorts, not mechanisms.
