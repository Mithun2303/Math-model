# Phase 4: Parameter Identifiability, Synthetic Validation and Clinical Endpoint Mapping

*All results in sections 1–5, 7–8 are SYNTHETIC: known ground truth, a hypothetical model, an assumed observation model. Section 6 is a descriptive audit of the cBioPortal METABRIC export; no model is fitted to or overlaid on it. Nothing here is clinical evidence, a validated measurement model, or a patient risk. No clinical calibration was performed.*

Reproduce: `.venv/bin/python experiments/{identifiability_analysis,synthetic_recovery,model_reduction,model_reduction_combinations,time_origin_sensitivity,survival_validation,population_identifiability_nuisance,metabric_endpoint_investigation}.py`. Design and **pre-registered decision criteria**: `configs/phase4.yaml` (sha256[:12] `cfb17cec5c80`, stamped in every Phase 4 table and figure). Tables `results/tables/p4_*.csv` (read with `comment="#"`), figures `results/figures/fig_p4_*.png`, each marked SYNTHETIC or CLINICAL (descriptive).

---

## 0. Guardrails and what was decided in advance

* The normalised P–D–e model (Phase 2; s = μ, tumour-free e = 1, K = 1 so burdens are fractions of K) is the reference. Reduced models are special cases in a separate experiment; no source model was modified (new code is additive: `identifiability.py`, `observation.py`, `synthetic_validation.py`, `endpoint_mapping.py`).
* **Proof standard.** *Full rank* of an exact series identifiability matrix at random rational points proves **generic local structural identifiability** (finite-to-one), not global identifiability. *Rank deficiency* is only evidence until certified by an explicit null vector. A numerically full-rank sensitivity matrix proves nothing structural. Practical identifiability is a statement about noisy, sparse, truncated data (profile likelihood, Raue et al. 2009).
* **Decision criteria were written before any result** (`configs/phase4.yaml`, `decision_criteria`): A (descriptive comparison), B (restricted exploratory population calibration), C (mechanistic calibration unsuitable). They are not mutually exclusive; each is answered in §9.

---

## 1. Observation model (Task 4; `src/dormancy/observation.py`)

Five levels are kept separate: (1) hidden states x(t) = (P, D, e); (2) tumour burden N = P + D (a mathematical quantity); (3) hypothetical **measurement signals** y_N = c·N·exp(σε), likewise y_P, y_D, y_e = c_e·e·exp(σε) (multiplicative log-normal noise; unknown or known scale c), values below a limit of detection (LOD) **left-censored with a Tobit likelihood, not dropped**; (4) the **threshold-defined recurrence event** (first upward crossing of N_detect), seen only through visits as an **interval-censored** record (L < T ≤ R); (5) an observed **clinical outcome** (time origin, detection delay, competing death; `endpoint_mapping.py`).

**Observation stops at detection.** Once N reaches N_detect the relapse is diagnosed and treated, so the untreated course is observed only up to T_detect (`truncate_at_detection`). "Long follow-up" therefore means observing up to detection (216 of 240 months for the G1 baseline), not N(240). Assumptions (all hypothetical): noise sd 0.05–0.5 (natural-log), LOD 10⁻⁶ (10⁻⁷ in the identifiability runs; 3×10⁻⁴ in one stress case), visits every 1–24 months, baseline visit at t = 0, signal scale c known (1) or unknown (2 in the stress case). Synthetic data are generated with a **different solver and tolerance** (Radau, rtol 10⁻¹²) than the fitter (DOP853, rtol 10⁻⁹, analytic forward sensitivities), avoiding an inverse crime.

---

## 2. Structural identifiability (Task 1; exact)

### Method
`src/dormancy/identifiability.py` computes the Taylor coefficients y_k of the output at t = 0 (series recurrence, including the division for P/(h+P)) with forward-mode automatic differentiation, in GF(p) (p = 2⁶¹−1) or exact rationals, at 3 random rational points with positive states and generic parameters (Taylor series verified against a 10⁻¹³ ODE solution to 10⁻¹⁶; AD gradients verified against finite differences). Rank is evaluated against Taylor order (flat over the last 6 of 30 orders at every point). Null vectors are exact rationals.

### Results (`p4_structural_identifiability.csv`; Fig 11 left)
Unknown sets A–F (7 parameters; + P₀, D₀; + h, μ; + e₀; + signal scale c; + K) × observation sets N, (P, D), (N, e), (P, D, e) (with c_e unknown when e is observed on an unknown scale):

* **Sets A–E (up to 13 unknowns, K = 1 fixed, c possibly unknown): full rank for all four observation sets**: every parameter, initial state and signal scale is **generically locally structurally identifiable, even from N(t) alone**.
* **Set F (K and c both unknown): deficiency exactly 1 for every observation set.** The exact rational null vector (`p4_structural_scale_gauge_null_vector.csv`) has relative components (−0.868, +0.868, +0.868, +0.868, −0.868, +0.868) on (η, P₀, D₀, h, c, K): the **burden-scale gauge** x → λx, K → λK, h → λh, η → η/λ, c → c/λ. This is the burden analogue of the Phase 2 E-gauge; fixing K = 1 (burden in units of K) removes it.
* **Low-burden sub-models are NOT identifiable** (`p4_structural_submodels.csv`, `p4_structural_linear_regime.csv`). For the linear block x′ = Mx with M = [[m, a],[q, −(a+δ_D)]], m = r − q − β e₀ and N-only observation, the output has four invariants (λ₊, λ₋, c₊, c₋) for six unknowns (m, a, q, s, P₀, D₀): symbolic rank 4 of 6, null dimension 2. In the model with r, β separate: **rank 4 of 7** for N alone (null dimension 3), **rank 6 of 7 for (P, D)** where the single deficiency is exactly **r, β appearing only as r − β**. The logistic term r·P·N/K and the immune equation are what separate r from β, and they are second-order in the burden.

### The decisive distinction
*Generic* structural identifiability (full model, analytic germ at a generic point) coexists with *non*-identifiability of the same parameters in the linear regime. The information separating r from β scales with the (tiny) nonlinear terms; whether it is usable is a **practical** question (§3). Global identifiability was not proven (a differential-algebra input–output analysis was not carried out).

### Event-only observation (scenario 5)
* **One scalar per individual.** The first-crossing time T(θ) is a smooth scalar function. Where it exists with a transversal crossing, ∇T ≠ 0 (G1: dT/dln r = −563, dT/dln β = +218, dT/dln q = +146, dT/dln f_D = +186, dT/dln N_detect = +51 months), so by the implicit-function theorem **every level set {T = const} is a manifold of dimension p − 1 = 9** in the ten inputs: no parameter is identifiable from one event time. Numerically: changing β by ±10–25% and re-solving r to keep T = 219.9 months exactly leaves N(60) within 2–5% and N(120) within 5–15% of the original curve (`p4_event_only_iso_T_family.csv`).
* **Linear-regime confounding (exact).** With T = (1/λ₊)·ln(N_detect/(n₁u₀)) in the single-mode limit, scaling N₀ by s, scaling the threshold, and shifting the time origin by τ all enter as additive log/linear terms and are **exactly confounded at a fixed λ₊**. Numerically (G1, λ₊ = 0.044): N₀ × 0.1, 0.3, 3, 10 shifts T by +52.37, +27.40, −25.06, −52.77 months vs the prediction −ln s/λ₊ = +52.33, +27.36, −24.97, −52.33 (`p4_event_only_N0_time_origin_confounding.csv`). Population heterogeneity in λ₊ breaks the confounding only partially (§7).
* For a **population**, event data identify the push-forward law of T, i.e. a distribution of one scalar, not the parameter distribution (§8).

---

## 3. Numerical and practical identifiability (Tasks 1–2)

### Fisher spectra (`p4_fisher_spectra.csv`; Fig 4)
Noise-scaled sensitivity matrices (analytic forward sensitivities, 9 log/logit parameters r, q, a, δ_D, β, ρ, η, N₀, f_D, truncated at detection, G1, 6-month visits): the number of directions with singular value > 3 (predicted sd < 0.33 in ln) is **4 of 9 for N, 6 for (P, D), 5 for (N, e), 7 for (P, D, e)** at σ = 0.10 (3, 6, 4, 6 at σ = 0.25). Truncation at detection costs at most one direction here, because the nonlinear regime is reached only after detection. For N alone the four determined directions are dominated by ln r with q and β (the growth-rate combination), then N₀, then mixtures; the weakest five involve δ_D, f_D, η, ρ. **Information on the r–β split does not increase with the initial burden**: the estimator correlation is 0.9998 and the Wald sd of ln r is 3×10²–2.6×10³ and of ln β 5×10²–4.6×10³ for N₀ scaled 0.1–100 (`p4_fisher_rbeta_vs_initial_burden.csv`).

### Profile likelihoods (G1 baseline, σ = 0.10, one data set per observed-state set; `p4_profile_likelihood_summary.csv`; Fig 3, 11 right)
Raue et al. criterion (95% interval finite on both sides within ±ln 10; profile grid step 0.33 in ln units, so interval widths of tightly determined parameters are approximate):

| parameter | N | P, D | N, e | P, D, e |
|---|---|---|---|---|
| r | one-sided | finite (×14 wide) | finite (×11 wide) | **finite (×1.08)** |
| q | not identifiable | **finite (×1.07)** | one-sided | **finite (×1.07)** |
| a | one-sided | **finite (×1.13)** | one-sided | **finite (×1.11)** |
| δ_D | not identifiable | one-sided | not identifiable | one-sided |
| β | not identifiable | not identifiable | one-sided | **finite (×1.59)** |
| ρ | not identifiable | not identifiable | finite (×14) | **finite (×1.74)** |
| η | not identifiable | not identifiable | not identifiable | **not identifiable** |
| N₀ | finite (×1.15) | finite | finite | finite |
| f_D | not identifiable | finite | not identifiable | finite |

**From total burden alone only the initial burden N₀ is practically identifiable** (profile finite on both sides); the other eight parameters are open on at least one side. Observing P and D separately adds r, q, a, f_D; adding the immune signal adds β and ρ. **δ_D and η are not practically identifiable in any observation set** for this ground truth. (All starts but one or two of 8 reached fits within Δχ² < 3.84 of the best, i.e. the likelihood surface is a ridge, not a set of distinct basins: 7/8 for N and P, D; 6/8 for P, D, e.)

---

## 4. Synthetic parameter recovery (Task 2; `p4_recovery_*.csv`; Figs 1, 5, 6, 7)

Design: one-factor-at-a-time around a base condition (ground truth G1, observed N, σ = 0.10, 6-month visits, follow-up limit 240 months truncated at detection, initial burden and dormant fraction unknown), 23 conditions × 8 noise replicates = 184 constrained fits, each from 8 random starts spread over box bounds (truth ÷30…×30, log scale); estimates on a bound are flagged. Fitted: r, q, a, δ_D, β, ρ, η, N₀, f_D (μ, h, K fixed; e₀ known).

**Base condition (observed N).** Median |ln error| per parameter: r 0.77, q 0.64, a 0.75, **δ_D 1.80, β 1.78, ρ 1.43, η 2.19, f_D 1.62**, N₀ 0.03; only N₀ is within a factor 1.5 in all replicates (r 25%, q 25%, a 50%, the rest ≤ 12%). **62.5% of fits end on a bound.** Yet the **crossing time is predicted within 10% in 87.5% of replicates** (median error 1.1%) and the in-window RMSE of ln N is 0.038: *good fits and good event-time predictions coexist with parameters wrong by factors of 5–10* (Fig 7; example seed 20261203: median parameter error e^2.5, T error 2.4%). The Wald coverage of the well-determined *directions* of the linearisation at the truth is only 8/32, because estimates sit far from the linearisation point; profile likelihoods are the reliable tool.

| condition (observed N unless stated) | directions with s > 3 | median \|ln err\| over params | bound-hit fraction | median T error |
|---|---|---|---|---|
| base σ 0.10 | 4 | 1.43 | 0.63 | 1.1% |
| σ 0.001 / 0.01 / 0.05 / 0.25 / 0.5 | 6 / 5 / 4 / 3 / 2 | 0.35 / 0.61 / 1.12 / 2.21 / 2.11 | 0 / 0 / 0.13 / 1.0 / 0.88 | 0.0% / 0.2% / 0.6% / 1.6% / 4.6% |
| schedule 1 / 3 / 12 / 24 months | 4 | 1.28 / 1.78 / 1.46 / 1.68 | 0.50 / 0.38 / 0.75 / 0.63 | 0.6% / 0.7% / 1.2% / 1.5% |
| follow-up limit 60 / 120 months | 2 / 3 | 2.91 / 1.70 | 1.0 / 1.0 | **43.5% / 26.8%** |
| observed (P, D) | 6 | 0.36 | 0.63 | 0.4% |
| observed (N, e) | 5 | 1.69 | 1.0 | 1.1% |
| **observed (P, D, e)** | **7** | **0.07** | 0.75 | 0.3% |
| initial state known | 3 | 1.33 | 0.38 | 0.7% |
| signal scale c unknown (true 2) | 4 | 1.42 | 0.38 | **25%** (in-window RMSE of ln 1.59 for the unscaled prediction) |
| high LOD 3×10⁻⁴ (early data censored) | 3 | 1.83 | 0.50 | 1.1% |
| truth G2 fast / G3 slow / G4 controlled / G5 low-dormant | 3 / 4 / 2 / 3 | 1.17 / 1.68 / 1.87 / 1.78 | 0.88 / 0.88 / 0.88 / 0.50 | 2.2% / **27.5%** / n.d. (no crossing) / 2.0% |
| **truth G6 (nonlinear transient)** | **5** | **0.09** | 1.0 | n.d. |

Findings.
1. **Noise.** Even at σ = 0.001 (near noise-free) the median error is 0.35: six directions are determined, three are not. Increasing noise degrades the determined directions steadily (6 → 2).
2. **Sampling frequency barely matters**; more visits (1-month: 220 observations) do not repair the weak directions: the missing information is in the dynamics sampled, not the sampling density.
3. **Truncated windows hurt prediction, not parameters in a useful way:** with a 60-month follow-up limit the predicted crossing time is wrong by 43%.
4. **Regime matters.** For the controlled regime (G4: decay only) 2 directions are determined. For **G6, whose trajectory enters the immune-feedback regime, r, β, ρ and η are recovered from N alone (|ln errors| 0.07, 0.09, 0.05, 0.09) but q and δ_D are not (3.4 each, the bound)**: which parameters are identifiable depends on which regime the data visit.
5. **Observed states.** Adding P and D separately recovers r, q, a, f_D (errors 0.02–0.36); adding e as well recovers β and ρ (0.07, 0.11). δ_D (1.41) and η (2.10) stay unidentified.
6. **Prediction is not recovery.** Prediction error was reported separately throughout; the "success" of predicting T (median 1–2% at the base condition) does not support any mechanistic parameter.

---

## 5. Redundancy and model reduction (Task 3; `p4_model_reduction_*.csv`; Fig 13)

Candidate reductions were derived from the exact null space and the Fisher directions, not by hand, and tested on synthetic data generated by the **full** model (6 replicates per cell, observed N and (P, D, e), truths G1, G2, G5, G6):
* **M7** immune frozen (ρ = η = 0), **M6** M7 with δ_D = 0, **M5** M6 with f_D fixed at the equilibrium ratio q/(q+a).

**Accept/reject (three criteria).** For **burden-only data in the sub-threshold regime (G1, G2, G5)**: M7, M6, M5 are preferred by AIC (median ΔAIC −3.5 to −6.4 vs the full model, 100% of replicates not worse than full + 2), and predict the crossing time equally well (median error 0.3–3.0% vs 0.3–1.4%) → *reduction justified by criteria (i) and (iii)*. But **criterion (ii) fails**: the individual parameters of the reduced models remain badly determined (median max Wald sd of ln parameters 2–484, above 30 in 10 of 24 cells), so **the reduction removes untestable parameters, not the redundancy** between r, β, q. For **(P, D, e) data and for the nonlinear G6**, the reductions are rejected decisively (ΔAIC +29 to +30,000; ρ, η matter).

**Which combination is recovered?** From N alone the combination recovered to 3% is the **dominant eigenvalue λ₊ of the low-burden block** (median |ln error| 0.032 G1, 0.030 G2, 0.028 G5, 0.004 G6), together with N₀ (3–6%); m = r − q − β is good to 5–9% (22% G6), whereas g = r − β (38%, 30%, 42%), λ₋ (ln errors 0.42, 0.49, 3.4) and the individual r, β, q, a, δ_D are poor (`p4_model_reduction_combinations.csv`). With (P, D, e): λ₊ to 0.7–1.2%, λ₋ to 3–12%, m to 1–7%.

**Proposed reduced parameterisation (justified only for the sub-threshold, burden-only regime):** report (λ₊, N₀) (plus λ₋ only with P, D data) as the identified quantities, and treat (r, β, q, a, δ_D, ρ, η, f_D) as nuisance parameters whose individual values are not identified. Mechanistic parameters are recoverable only when P, D (and e) are observed or when the trajectory enters the nonlinear regime. The original model remains the reference implementation.

---

## 6. METABRIC RFS endpoint: what is verified (Task 5; descriptive; `p4_metabric_*.csv`)

**Documentation search (time-boxed).** Export metadata (verbatim): RFS_STATUS "Recurred = cases where there was loco-regional relapse, distant relapse or death-specific death"; RFS_MONTHS "Time until last follow-up or relapse (in months)" (**no origin stated**); OS_MONTHS "Overall survival in months **since initial diagnosis**". Europe PMC open full text of Curtis 2012 (PMC3440846), Pereira 2016 (PMC4866047) and Rueda 2019 (PMC6647838): Curtis and Pereira report disease-specific survival; Rueda defines a multistate model whose disease-free state starts **at surgery** ("transition probabilities from surgery", "relapse free five years post diagnosis" used loosely) and says nothing about the cBioPortal RFS field. cBioPortal datahub history of `data_clinical_patient.txt` shows "Added RFS data" (2020-12-01), "Update RFS STATUS description" (2020-12-02), "update survival data" (2024-04-22) with no source named (the patches are Git-LFS pointers); the Europe PMC supplementary archive for PMC6647838 contains one image. **No source gives the RFS time origin, treatment dates, or a separation of recurrence from death.**

**Common-origin test (a stronger check than RFS ≤ OS; tests a shared clock, not which origin).** If the two clocks share an origin, patients whose RFS and OS clocks both end at the same last contact must have RFS = OS exactly: **745 of 745** "Not Recurred & Living" and **432 of 432** "Not Recurred & Died of Other Causes" patients have RFS_MONTHS = OS_MONTHS to 10⁻⁶; across 1,980 patients with both times there are 0 with RFS > OS. This is strong evidence of a common clock (INFERENCE), and still consistent with diagnosis **or** surgery (days to weeks apart, unrecorded) as the origin; it does not say the origin is treatment end.

**Composition of the composite event (964 usable events).** 92 Recurred & Living (9.5%); 65 recurred with other-cause death (6.7%); **545 died of disease with RFS < OS** (recurrence preceded death; 56.5%); **101 died of disease with RFS = OS** (death possibly the first recorded event; 10.5%); **161 with missing vital status (16.7%)**. Documented relapse before death or no death: 72.8%. Upper bound on the death-first share: 10.5% (known status) to 27.2% (if all missing-status events were death-first); lower bound 0. Recurrence and disease-specific death are **not separable**, only bounded. Among the 545, OS − RFS has median 15.5 months (IQR 6.0–34.1), conditional on death (biased; not a delay distribution).

**Treatment timing.** CHEMOTHERAPY, HORMONE_THERAPY, RADIO_THERAPY, BREAST_SURGERY are binary/categorical flags (2,509 − 529…554 non-missing); **no dates or durations**: the diagnosis-to-treatment interval τ cannot be reconciled and remains a nuisance.

**Subgroups** (`p4_metabric_subgroup_reliability.csv`; criteria ≥ 50 events and ≥ 100 at risk at 60/120 months or ≥ 30 at 240 months): reliable at 5 years for all listed groups; at 10 years except Her2 (77 at risk), Basal (79) and Normal (61, 65 events); at 20 years only all-usable (162), ER+ (110), ER− (49), LumA (61); LumB (27), Her2 (18), Basal (13), claudin-low (17), Normal (10) are not interpretable at 20 years. 21% of patients lack the extended subgroup fields; 3.3% lack ER.

**Verified / inferred / unknown:** VERIFIED: OS origin = diagnosis; composite definition; no treatment dates; no loco-regional/distant split; common-clock equality (computed). INFERENCE: common origin of OS and RFS. UNKNOWN: RFS origin, diagnosis-to-treatment interval, death-first share, recurrence-to-detection delay.

---

## 7. Time-origin and detection-delay sensitivity (Task 6; `p4_time_origin_*.csv`, `p4_threshold_effects.csv`, `p4_equivalence_map.csv`; Figs 9, 10, 12)

Hypothetical mapping t_clin = τ + t_model + d (recurrence) or τ + t_model + L (death-first); **nuisance assumptions, not estimates**. Population: prior P1 run A (n = 5000, N_detect 0.1).

* **Shift is exact and additive:** τ and d only enter through τ + d: the whole curve shifts by τ + d (time to 10% cumulative probability 125 → 128/131/137/143/149 months for τ + d = 3/6/12/18/24), while F(120 mo) falls 0.085 → 0.077/0.066/0.051/0.038/0.025 and F(240 mo) falls 0.004–0.035.
* **Detection threshold 0.01 → 0.5:** F(240) 0.689 → 0.214, time to 10% 76 → 182 months (no sample starts above any threshold).
* **Initial burden ×0.1…×10:** numerical median shift matches the linear-regime law −ln s/λ₊ per sample (−18.76 vs −18.66 at ×3; −39.51 vs −39.11 at ×10; +36.68 vs +36.62 at ×0.1).
* **Death-first share** 0.105 / 0.272 (lag mean 24 months): time to 20% moves 150 → 153/155 months; F(240) 0.424 → 0.420/0.415: small.
* **Equivalence map** (sup-norm distance over [0, 240] after the best compensating change, vs DKW sampling noise ε = 0.0278 for n = 2380, 0.0192 for n = 5000): τ = 12 months is reproduced to 0.0092 by an initial-burden rescale ×0.5 (below ε for both n) and to 0.0224 by a threshold of 0.136 (below ε for n = 2380); τ = 24 by N₀ ×0.2 to 0.0250 (below ε for n = 2380). A factor-3 change in N₀ equals a shift of 19 months (N₀/3: distance 0.0134, below ε) or −28 months (N₀×3: 0.0450, above ε); a factor-3 threshold equals +49 months (0.0490, above ε). **At METABRIC sample size, an offset of up to about 12 months cannot be separated from a factor-2 change in the initial burden using event times alone.** Larger offsets (≥ 24 months, or factor-3 threshold changes) are partly separable *only if the other nuisances are known*.

---

## 8. Synthetic survival validation (Task 7; `p4_survival_*.csv`; Fig 8)

Truth: the pooled prior-P1 population (3 independent runs, 15,000 individuals; F(60/120/240) = 0.0005 / 0.0910 / 0.4363). 200 resampled cohorts of n = 2380 (and n = 1000).

* **Right censoring handled correctly.** Kaplan–Meier under administrative censoring at 60/120/240 months and under additional independent loss to follow-up is unbiased (|bias| ≤ 0.0007, ≈ 1 MC SE); 95% coverage 0.94–0.97 at 120 and 240 months. At 60 months coverage is 0.65–0.69 because the truth (0.0005) gives about one event per cohort and the interval is degenerate (not meaningful).
* **Interval censoring.** Turnbull with annual or 6-monthly visits is unbiased (|bias| ≤ 0.0007). **Turnbull interval coverage at n = 2380 was not computed** (bootstrap was run only for n = 1000, 50 cohorts).
* **Naive estimators** that treat censored records as event-free are biased low (−0.017 at 120, −0.116 at 240 months). Midpoint imputation was essentially unbiased in this setting (not a useful contrast).
* **Event definition.** Analysing persistent exceedance when the target is first crossing biases F by −0.0024 (120 mo) and −0.0042 (240 mo) (small but ≈ 10 MC SE); against the persistent target the bias is +0.0008: a pure definition effect. Early-heavy, late-heavy and mostly-no-crossing populations: |bias| ≤ 0.0013.
* **Population hyper-parameters.** For ln β ~ N(m, s²) (m = ln 0.04, s = 0.3), n = 2380, annual visits, 50 replicates, with everything else known: m bias −0.0007 (RMSE 0.0078; Wald coverage 0.94, profile 0.92), ln s bias +0.0038 (RMSE 0.022; coverage 0.96/0.96): **identifiable**. T(β) is monotone increasing, so the event-time law determines the β law.
* **Non-identifiability appears with a second varying parameter.** If ln r also varies (4 unknowns), the profile of m_β is a **ridge**: 11 of 17 grid points lie inside the 95% region (m_β spanning ≈ −3.62 … −2.62, about 1 log unit, versus ±0.014 when r is known); the best fit has m_β = −2.73 (truth −3.22) while the event-time CDF is fitted to a sup-norm error of 0.030 and moves by at most 0.009 along the ridge. Event times identify a distribution, not mechanisms.

### Pre-registered criterion B tested: population parameters with the nuisances left free (`p4_population_nuisance*.csv`)

Synthetic population ln β ~ N(m, s²) (m = ln 0.04, s = 0.3), true offset τ = 12 months and true initial-burden scale 1, n = 2380, annual visits in clinical time to 240 months (one structural source of heterogeneity, everything else known and the model structure exactly correct).
* **τ free (m, s, τ):** all three are estimable. Over 30 replicates the estimates are m̂ = −3.217 (truth −3.219; sd 0.024), ŝ = 0.300 (sd 0.027) and τ̂ = 11.6 months (truth 12; sd 4.7); the replicate correlation of m̂ and τ̂ is −0.947. An analyst who wrongly fixes τ = 0 gets biased location and scale (m̂ +0.051, ŝ −0.062; the likelihood-ratio statistic is 5.3 larger than with τ free on average). The nonlinear shape of T(β) distinguishes a shift from a change of m.
* **τ and the initial-burden scale s₀ both free (m, s, τ, s₀; 8 replicates, smoothed pushforward through a 2-D table T(ln β, ln s₀)):** m̂ = −3.230 (sd 0.055), ŝ = 0.310 (sd 0.045): **m and s still have closed profile intervals** (rep 0: m ∈ [−3.38, −3.08], i.e. median β within ±16%; ln s ∈ [−1.65, −0.96], s ∈ [0.19, 0.38]). But **the nuisances are not identified**: τ ∈ [−7, +46] months (sd across replicates 9.2), and the initial-burden scale is **open** (the whole window ×0.22–×4.5 lies inside the 95% region); τ̂ and ln ŝ₀ are correlated 0.88. This is the exact confounding of §2 and §7 showing up in the likelihood.
* **A second varying parameter destroys it** (§8 above): with ln r also varying, m_β is a ridge of width ≈ 1 log unit even when τ and s₀ are known.
Interpretation: criterion B as literally written (a finite profile interval for some population combination with τ and the threshold/N₀ free) is met **only** in the idealised setting of one known heterogeneity source and an exactly known model. Which sources of heterogeneity exist in real patients cannot be checked from event times (the data never show P, D or e), and with two sources the criterion fails. The decision below uses that judgement, which was made after seeing these results.


---

## 9. Clinical feasibility decision (Task 8; criteria pre-registered in `configs/phase4.yaml`)

**A. Descriptive population-level comparison — YES.** The endpoint is documented (composite RFS, non-disease deaths censored), 2,380 usable records and 964 events (≥ 500), reliable subgroup estimates exist at 5 and 10 years, and uncertainty can be stated. Any comparison must be labelled *consistency, not validation*, name the target "RFS as defined in the cBioPortal METABRIC export", use 5 and 10 years as primary horizons, report 15–20 years with at-risk counts, and carry τ, d, the threshold and the initial-burden scale as explicit nuisances.

**B. Restricted exploratory population calibration — **NOT SUPPORTED (stop at A).** Pre-registered criterion: at least one combination of population parameters with a finite 95% profile interval from interval-censored event times when τ and the threshold/N₀ scale are free, and not a reparametrisation of the nuisance. Evidence: (i) met in the idealised one-source test (m, s finite with τ and N₀-scale free; τ itself spans ≈ 53 months and the N₀ scale is open); (ii) the event-time curves of τ = 12 months and N₀ ×0.5 are within DKW noise at n = 2380 (§7); (iii) with a second varying parameter the location parameter is a ridge of ≈ 1 log unit (§8); (iv) METABRIC cannot verify the heterogeneity structure and has an undocumented origin. A defensible reading is that only a **phenomenological** description of the event-time distribution (location, scale, cure fraction) can be estimated, and that is a descriptive summary (A), not mechanistic calibration. B would require the project owner's separate approval and a pre-specified phenomenological scope.**

**C. Mechanistic calibration — UNSUITABLE.** Criterion: the mechanistic parameters (r, β, q, a, ρ, η, initial composition) must be practically identifiable from burden-type or event-type data under the stated observation models. They are not: from burden alone only N₀ (and the combination λ₊) is identified; r, β, q, a, δ_D, ρ, η and f_D are open (§3); an event time gives one scalar and every parameter lies on a 9-dimensional iso-T manifold (§2); the time-origin, threshold and initial-burden nuisances are exactly confounded in the linear regime and only partially separable otherwise (§7); METABRIC provides event times only (never P, D or e), with an undocumented RFS origin and no treatment dates (§6). Observing P, D and e would make r, q, a, β, ρ and f_D identifiable (§3–4) but the data cannot supply it.

**Recommendation.** Stop at descriptive comparison (A) unless B passes under its pre-registered conditions.

---

## 10. Failure cases and limitations

1. **Structural identifiability is local and generic**, tested at 3 random rational points per set with rank saturation over orders; global identifiability and special (non-generic) regimes beyond the linear and frozen-immune sub-models were not analysed with a differential-algebra method.
2. **Practical identifiability depends on the assumed observation model** (noise level, schedule, truncation at detection, unknown scale). A different, richer biomarker would change the answer; none exists in the clinical data.
3. **Fits:** 8 random starts per fit; ridge-like likelihoods make "global optimum" statements unreliable (most starts are within Δχ² < 3.84 of the best); bounds are hypothetical; 62% of base-condition fits ended on a bound. The Gauss–Newton Fisher information ignores the Tobit curvature; Wald/linearised diagnostics are poorly calibrated far from the truth (8/32 well-determined-direction intervals covered) and the profile grid (step 0.33 in ln) limits interval precision for tight parameters.
4. **Replicates are few** (8 per condition; 6 for model reduction) and the ground truths are 6 hypothetical parameter sets; results are conditional on them (e.g. G6 shows a different identifiable set).
5. **Survival validation:** truth is the pooled prior-P1 empirical law (n = 15,000, sampling error ≈ 0.004); 60-month results are degenerate (truth 0.0005); Turnbull coverage at n = 2380 not computed; informative censoring not studied.
6. **Time-origin analysis** uses a single hypothetical population (prior P1) and a single death-lag assumption; nuisance windows are not estimates. The equivalence map can only show what cannot be separated, not what τ is.
7. **METABRIC:** RFS origin and treatment timing remain unknown after a time-boxed documentation search (supplementary methods of the original papers beyond the open full text were not accessible); 161 events have missing vital status, which limits the death-versus-recurrence bounds; 21% of patients lack the extended fields.
8. **Not done:** clinical parameter calibration (out of scope by instruction), stochastic birth–death–switching model, competing-risk (Aalen–Johansen) analyses of METABRIC, loco-regional versus distant separation, any use of METABRIC genomics.
