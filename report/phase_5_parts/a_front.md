## 1. Introduction

Late recurrence of breast cancer, years or decades after apparently successful primary treatment, is well documented clinically (Pan et al. 2017; Rueda et al. 2019). One mechanistic hypothesis is that some disseminated cells enter a non-proliferating, dormant state and are later reactivated, or escape from immune control (Aguirre-Ghiso 2007; Wilkie and Hahnfeldt 2013). Phases 1–4 of this project built and analysed a hypothetical three-compartment ordinary differential equation (ODE) model of that idea. The state variables are a proliferating population $P$, a dormant population $D$ and an effective immune activity $E$:

$$
\begin{aligned}
\dot P &= rP\left(1-\frac{P+D}{K}\right) - qP + aD - kEP,\\
\dot D &= qP - (a+\delta_D)D,\\
\dot E &= s + \rho\,\frac{P}{h+P}E - \mu E - \eta (P+D)E .
\end{aligned}
\tag{1.1}
$$

The previous phases established the following.

- **Phase 2 (mathematics).** The model has an exactly characterised equilibrium structure and a verified transcritical bifurcation at $k=k_c$, and the relapse-time estimates have documented validity limits.
- **Phase 3 (uncertainty).** Simulated relapse probabilities depend strongly on the assumed priors and on the detection threshold. Net proliferation $r$ and immune killing $\beta$ dominate the sensitivity.
- **Phase 4 (identifiability).** From burden or event-time data most mechanistic parameters are not practically identifiable. Clinical *calibration* of the mechanism was therefore judged unsuitable, and only a descriptive comparison (scope A) was recommended.

Phase 5 carries out that descriptive comparison. It has two parts. The first is a descriptive Kaplan–Meier analysis of the public METABRIC recurrence-free survival (RFS) data. The second is a transparent, pre-specified comparison of those observed event-time patterns with model-generated event-free survival curves.

**What Phase 5 is not.** It is not a clinical validation, not a prediction of any patient's risk, and not a calibration of biological parameters. All model parameters remain hypothetical. "Consistency" here means only that an observed curve does or does not lie within the range of curves the hypothetical model family can produce under declared assumptions.

## 2. Research Objectives

| # | Objective | Success criterion |
|---|---|---|
| P5-O1 | Describe observed METABRIC RFS at 5, 10, 15 and 20 years | Kaplan–Meier estimates with 95% CIs, events and numbers at risk, from the verified cohort |
| P5-O2 | Describe ER-positive versus ER-negative RFS and molecular subgroups | Only strata meeting the pre-existing reliability criterion (Section 6.4) are interpreted |
| P5-O3 | Describe late-recurrence patterns | Interval-specific event rates, landmark analysis; no causal reading |
| P5-O4 | Generate model-defined event-free survival curves from the existing Phase 3 and 4 simulations | Reuse of the validated ODE and the prior-P1 population; no new fitting |
| P5-O5 | Compare simulated and observed curves descriptively | A scenario grid and comparison criteria recorded *before* the comparison is examined |
| P5-O6 | Quantify the effect of the unknown time origin, the detection threshold and observation delay | Sensitivity plots; findings that survive all assumptions identified |
| P5-O7 | State limits of interpretation | Explicit separation of verified facts, inferences and hypotheses |

These objectives address Phase 1 research question RQ5, "can the model reproduce group-level recurrence patterns?", and its hypothesis H5. H5 was only stated, never tested, before this phase. It is read here as consistency, not as evidence of mechanism, exactly as formulated in `research_gap_and_objectives.md`.

## 3. METABRIC Dataset Description

**Source.** The data are the cBioPortal export of the METABRIC study (study identifier `brca_metabric`; "Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016)"), downloaded on 2026-11-03 from the cBioPortal datahub (cBioPortal datahub 2026). The underlying cohort is described by Curtis et al. (2012) and Pereira et al. (2016); late-recurrence dynamics in the same cohort are described by Rueda et al. (2019). The files are kept unmodified in `data/raw/brca_metabric/`.

**Licence.** The repository licence is the Open Database License (ODbL 1.0), which requires attribution and share-alike for derived databases.

**Files and size.** The patient table has 2,509 patients and 24 columns. The sample table has 2,509 samples (all primary tumours) and 13 columns. They join one-to-one on `PATIENT_ID`. All 2,509 patients are female. A 36-field dictionary with units, missingness and event semantics is in `data/clinical_data_dictionary.csv`.

**Fields used in Phase 5.**

| Field | Meaning in the export | Use |
|---|---|---|
| `RFS_STATUS` | "Recurred" = loco-regional relapse, distant relapse or disease-specific death; "Not Recurred" otherwise | event indicator |
| `RFS_MONTHS` | "Time until last follow-up or relapse (in months)"; no origin stated | time |
| `OS_STATUS`, `OS_MONTHS`, `VITAL_STATUS` | overall survival, documented "since initial diagnosis" | consistency checks and death-versus-recurrence bounds only |
| `ER_IHC` | oestrogen-receptor status by immunohistochemistry | subgroup (ER-positive, ER-negative) |
| `CLAUDIN_SUBTYPE` | PAM50 and claudin-low subtype | subgroup |

**Why TCGA was excluded.** The TCGA-BRCA PanCancer Atlas clinical export (1,084 patients; TCGA cBioPortal datahub 2026) was audited as an alternative in Phase 1. It has only 10 disease-free-survival events after 60 months, so it cannot inform late-recurrence behaviour. The two datasets also use different endpoints, so they are not combined.

## 4. Endpoint Definitions

### 4.1 Clinical endpoint: composite RFS

The exported RFS event combines loco-regional relapse, distant relapse and disease-specific death. Non-disease deaths are censored (432 of 497 other-cause deaths are "Not Recurred"; the other 65 are recorded as recurred). The endpoint is therefore recurrence-free survival with other-cause deaths censored. It is not overall survival and not disease-free survival including all deaths. The wording of endpoint definitions differs between trials (Tolaney et al. 2021), so the target of every comparison is named "RFS as defined in the cBioPortal METABRIC export".

Three facts are verified (Phase 4, `p4_metabric_*.csv`):

- The time origin of RFS is **not documented** in the export, in the open text of the three METABRIC papers, or in the datahub history.
- There are no treatment or surgery dates.
- Loco-regional and distant relapse cannot be separated, and recurrence cannot be separated from disease-specific death. Of the 964 events, at most 10.5% (known vital status) to 27.2% (if all 161 events with missing vital status are counted) could be death-first.

One inference is supported: RFS equals OS exactly for 745 of 745 "Not Recurred & Living" and 432 of 432 "Not Recurred & Died of Other Causes" patients, with no case of RFS > OS in 1,980 patients with both times. This indicates a *shared clock* for the two endpoints. It is still consistent with a diagnosis origin or a surgery origin, and does not show that the origin is the end of treatment.

### 4.2 Simulated relapse event

In the model, relapse is the **first upward crossing** of a model-defined burden threshold: $N(t)=P(t)+D(t)$ reaches $N_{\mathrm{detect}}$ (baseline 0.1 of carrying capacity $K$), starting below it, and the time is right-censored at the horizon (240 months by default). The threshold is a mathematical design choice and is not an imaging limit. The model time $t_{\mathrm{model}}$ starts at the post-treatment residual-disease state $(P_0,D_0,e_0)$.

### 4.3 Mapping from model time to clinical time

Following Phase 4 (`endpoint_mapping.py`), the clinical event time is represented by

$$
t_{\mathrm{clin}} = \tau + t_{\mathrm{model}} + d ,
\tag{4.1}
$$

where $\tau\ge 0$ is the unknown offset between the RFS time origin and the start of the model clock, and $d\ge 0$ is a detection or observation delay. Both are **nuisance assumptions, not estimates**. Only the sum $\tau+d$ enters, and it shifts the whole curve exactly.

### 4.4 Why the two events are not equivalent

1. **Composition.** Clinical RFS includes loco-regional relapse and disease-specific death; the model event is a burden crossing of a single lumped population.
2. **Time origin.** The clinical origin is unknown (diagnosis or surgery are the candidates supported by documentation); the model origin is a residual-disease state after treatment.
3. **Threshold.** Clinical detection depends on imaging, symptoms and follow-up schedule; the model threshold is arbitrary. Varying it from 0.01 to 0.5 moves the simulated 20-year event probability from 0.689 to 0.214 (Phase 3).
4. **Mechanisms absent from the model.** Treatment, early metastatic recurrence, and patient and tumour heterogeneity beyond the sampled priors are not represented.
5. **Variability.** Simulated spread is parameter and initial-condition uncertainty only, not stochastic extinction or demographic noise.

For these reasons the comparison is a statement about the *range of shapes* the hypothetical model family can produce, not a test of the mechanism.

## 5. Clinical Preprocessing

The preprocessing reuses the documented Phase 1 and Phase 4 rules. No imputation is applied.

1. Read `data_clinical_patient.txt` (tab-separated, comment lines skipped). Convert `RFS_MONTHS` and `OS_MONTHS` to numeric.
2. Code the event as 1 if `RFS_STATUS` starts with "1" ("1:Recurred") and 0 if it starts with "0" ("0:Not Recurred"); missing status stays missing.
3. Keep patients with both a usable `RFS_MONTHS` and a non-missing event indicator: **2,380 of 2,509** (964 events, 1,416 censored). The 129 excluded patients are 13 missing both fields, 108 missing time only and 8 missing status only. These counts are checked again by the Phase 5 code against the Phase 1 audit.
4. Retain records with RFS time zero (4 patients, 3 of them events) as they are.
5. Define subgroups from the same table: ER-positive and ER-negative by `ER_IHC`; PAM50 and claudin-low strata by `CLAUDIN_SUBTYPE`. About 21% of usable patients lack the extended subgroup fields and about 3.3% lack ER status; they drop out of subgroup analyses only. The groups overlap across the two families (ER and PAM50), so they are never pooled in one test.
6. The raw files are never modified. Every output table carries a provenance line with the SHA-256 prefix of the raw patient file.

## 6. Kaplan–Meier Methodology

### 6.1 Estimator

For event times $t_1<t_2<\dots$ with $d_i$ events and $n_i$ subjects at risk,

$$
\hat S(t)=\prod_{t_i\le t}\left(1-\frac{d_i}{n_i}\right),
\qquad
\widehat{\mathrm{Var}}\,\hat S(t)=\hat S(t)^2\sum_{t_i\le t}\frac{d_i}{n_i(n_i-d_i)} \quad\text{(Greenwood)} .
\tag{6.1}
$$

The reported quantity is the cumulative incidence of the RFS event, $1-\hat S(t)$. Because other-cause deaths are censored, this is not a competing-risk cumulative incidence and tends to overstate recurrence relative to an Aalen–Johansen estimate. A conventional competing-risk analysis is not performed, since distinct event causes and their times are not available.

### 6.2 Confidence intervals

Pointwise 95% intervals use the log–log transformation, $\hat S^{\exp(\pm z\,\hat\sigma_{\mathrm{ll}})}$ with $\hat\sigma_{\mathrm{ll}}=\sqrt{\widehat{\mathrm{Var}}\,\hat S}\,/\,(\hat S\,|\ln\hat S|)$ and $z=1.96$. This keeps the interval inside $[0,1]$. The own implementation (`km_descriptive.py`) was checked for internal consistency against the Phase 1 hand-written estimator; no external survival library is installed.

### 6.3 Horizons and censoring

Estimates are reported at 5, 10, 15 and 20 years (60, 120, 180 and 240 months). The number at risk at horizon $t$ is the number of patients with observed time $\ge t$. "Administrative censoring at the horizon" means that events after the horizon are censored at the horizon; it is used so that log-rank tests compare groups over the same window.

### 6.4 Reliability criterion for subgroups

The criterion was fixed in Phase 4, before any Phase 5 curve was drawn. A group is interpreted at a horizon only if it has at least **50 events** in total and at least **100 patients at risk** at 5 and 10 years, or at least **30 at risk** at 20 years. A group that fails is listed with the reason and is not interpreted. No criterion was defined for 15 years, so 15-year results are reported for the whole cohort only.

### 6.5 Group comparisons

Differences between strata use the Mantel–Cox log-rank test, restricted to strata reliable at the horizon and administratively censored there. For more than two strata, pairwise tests are adjusted with the Holm step-down procedure within each family and horizon. The ER family and the PAM50 family are tested separately. The p-values are descriptive: the strata overlap across families and no adjustment is made across families or horizons.

### 6.6 Landmark analysis

To describe late recurrence among patients who have already been event-free for five years, the analysis conditions on being event-free and still under observation at 60 months and computes the Kaplan–Meier cumulative incidence over the following 5, 10 and 15 years. The landmark removes early-event patients from the risk set, so its estimate is conditional and is not comparable to unconditional curves.

### 6.7 Interval-specific event rates

For descriptive late-recurrence patterns, the events in successive intervals (for example 0–5, 5–10, 10–15 and 15–20 years) are divided by the person-time at risk in the same interval. The interval hazard is $\hat\lambda=d/T$ for $d$ events and $T$ person-years. An exact (Garwood) 95% Poisson interval is

$$
\left[\frac{\chi^2_{0.025,\,2d}}{2T},\ \frac{\chi^2_{0.975,\,2d+2}}{2T}\right].
\tag{6.2}
$$

These rates are reported for description only. They mix recurrence and disease-specific death and depend on the unknown time origin, so no biological process is inferred from them.
