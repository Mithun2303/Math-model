# Phase 4 Summary: Identifiability, Synthetic Validation and Clinical Endpoint Mapping

*SYNTHETIC results (known ground truth, hypothetical model, assumed observation model) plus a descriptive METABRIC audit. No clinical calibration was performed; nothing here is a patient risk. Full methods, numbers and failure cases: `report/phase_4_identifiability_and_validation.md`. Design and pre-registered decision criteria: `configs/phase4.yaml` (sha256[:12] `cfb17cec5c80`).*

## Tests
**138 passed, 0 failed** (`.venv/bin/python -m pytest`): the 111 Phase 1–3 tests unchanged + 27 new (`tests/test_phase4_identifiability.py`, `test_phase4_observation_synthetic.py`, `test_phase4_endpoint_survival.py`). They include Taylor series vs a 10⁻¹³ ODE solution, forward sensitivities vs finite differences, the linear sub-model rank 4 of 7 and r − β redundancy, the K–c scale gauge (deficiency 1), toy non-identifiable (a·b·x) and identifiable (logistic) models, GF(p) vs exact rank, synthetic-data reproducibility, noise statistics, LOD/Tobit, truncation at detection, interval/right censoring, noise-free recovery within 1%, parameter-bound enforcement, solver-tolerance convergence, time-origin transformations, Kaplan–Meier worked example, Turnbull, DKW.

## Files created
- `configs/phase4.yaml` (ground truths, observation assumptions, fitting design, seeds, **pre-registered decision criteria**)
- `src/dormancy/identifiability.py` (exact series/AD identifiability matrices over GF(p) and Q, null vectors, forward sensitivities, Fisher spectra, event-time gradient, linear-regime symbolic analysis), `observation.py`, `synthetic_validation.py` (residuals with analytic Jacobians and Tobit, multistart, Wald/Fisher, profile likelihood, recovery/prediction metrics, Kaplan–Meier, Turnbull, hyper-parameter likelihood), `endpoint_mapping.py`
- `experiments/_p4common.py`, `identifiability_analysis.py`, `synthetic_recovery.py` (+`--figures-only`), `model_reduction.py`, `model_reduction_combinations.py`, `time_origin_sensitivity.py`, `survival_validation.py`, `population_identifiability_nuisance.py`, `metabric_endpoint_investigation.py`
- the three Phase 4 test files; `report/phase_4_identifiability_and_validation.md`; this file
- `results/tables/p4_*.csv` and `results/figures/fig_p4_01 … 13*.png` (figures 1–12 as requested plus a model-reduction figure; each marked SYNTHETIC or CLINICAL)
- Files modified: `README.md` (usage). No Phase 1–3 source, table, figure or config was changed.

## Structurally identifiable (exact; generic LOCAL identifiability, not global)
- **All of r, q, a, δ_D, β, ρ, η, h, μ, P₀, D₀, e₀ and the signal scale c are generically locally structurally identifiable from total burden N(t) alone** (and for every other observation set), when K = 1 (burden in units of K). Rank is full at 3 random rational points and flat over the last 6 of 30 Taylor orders.
- **Exception: when both K and the signal scale c are unknown, deficiency exactly 1**: the burden-scale gauge x → λx with (K, h, P₀, D₀) ×λ and (η, c) ÷λ (exact rational null vector).
- **Not identifiable (low-burden/linear regime):** four invariants for six unknowns from N alone (symbolic rank 4 of 6); r and β enter only as r − β (rank 6 of 7 even from P, D).
- **Event-only:** one scalar per individual; every level set of the first-crossing time is a 9-dimensional manifold, so no parameter is identifiable from one event time. N₀-scale, threshold-scale and time-origin shift τ are exactly confounded in the linear regime (numerical shift matches −ln s/λ₊ to ≈ 1%).

## Practically identifiable (profile likelihood, G1 baseline, σ_ln 0.10, 6-month visits, truncated at detection)
| observed | practically identifiable (95% interval finite both sides) | not identifiable / open |
|---|---|---|
| N | **N₀ only** | r (one-sided), q, a (one-sided), δ_D, β, ρ, η, f_D |
| P, D | r (wide, ×14), q, a, N₀, f_D | δ_D (one-sided), β, ρ, η |
| N, e | r (wide), ρ (wide), N₀ | q, a, δ_D, β (one-sided), η, f_D |
| P, D, e | r, q, a, β, ρ, N₀, f_D (intervals ×1.07–1.74) | **δ_D, η** |
- **δ_D and η are practically non-identifiable in every observation set.**
- From burden alone the recoverable quantities are the **dominant eigenvalue λ₊ of the low-burden block (recovered to 3%)** and N₀ (3%); m = r − q − β to 7%; the guessed "r − β" and the individual parameters are not.

## Synthetic recovery
- Base condition (observed N): median |ln error| 1.4 (a factor ≈ 4); 62% of fits end on a bound; yet the **crossing time is predicted within 10% in 88%** of replicates (median 1.1%). **Good fits and good event-time predictions coexist with parameters wrong by factors of 5–10.**
- Observed (P, D, e): median |ln error| 0.07 (r, q, a, β, ρ, f_D recovered; δ_D, η not). Observed (P, D): 0.36. Observed (N, e): 1.7.
- Regime matters: for the nonlinear transient case (G6) r, β, ρ, η are recovered from N alone (|ln error| 0.05–0.09) but q and δ_D are not; for decaying (controlled) truths only 2 directions are determined.
- **Noise:** 6 → 2 well-determined directions from σ = 0.001 to 0.5; even at σ = 0.001 the median error is 0.35. **Sparse sampling** (1 → 24 months) changes little. **Short follow-up** is the main prediction hazard: with 60/120-month limits the predicted crossing time is wrong by 43%/27%. Unknown signal scale: prediction error 25%.
- Model reduction (immune frozen, δ_D = 0, f_D fixed) is preferred by AIC and predicts T equally well for burden-only sub-threshold data, but leaves individual parameters undetermined; it is rejected for (P, D, e) data and the nonlinear regime (ΔAIC up to 3×10⁴). Proposed reduced parameterisation: report (λ₊, N₀) and treat the rest as nuisance.

## Synthetic survival validation
Kaplan–Meier and Turnbull recover the true cumulative incidence without bias (|bias| ≤ 0.0007; KM coverage 0.94–0.97 at 10 and 20 years); naive "censored = event-free" is biased low (−0.017, −0.116); using the wrong event definition (persistent vs first crossing) biases by −0.002 to −0.004. A one-parameter population (ln β ~ N(m, s²)) is recovered (m bias −0.0007, ln s +0.004; coverage ≈ 0.94–0.96) while with two varying parameters (β, r) the profile of m_β is a ridge spanning ≈ 1 log unit although the event-time CDF is fitted to 0.03. Event times identify a distribution, not mechanisms. (Not computed: Turnbull coverage at n = 2380; 60-month results are degenerate.)

## METABRIC endpoint feasibility (descriptive audit)
- **Verified:** RFS is a composite (loco-regional relapse, distant relapse, disease-specific death); OS origin is diagnosis; RFS origin is **not documented** anywhere found (export metadata, the three METABRIC papers' open text, the datahub commit history); no treatment dates or durations; loco-regional and distant relapse cannot be separated.
- **Common clock (inference, not proof of which origin):** RFS = OS exactly for 745/745 alive-not-recurred and 432/432 other-cause-death patients; 0 violations in 1,980.
- **Death vs recurrence:** of 964 events, 72.8% have documented relapse before death or no death; death may be the first recorded event in at most 10.5% (27.2% counting the 161 events with missing vital status). Not separable, only bounded.
- **Subgroups:** reliable at 5 years for all listed groups, at 10 years except Her2/Basal/Normal, at 20 years only all-usable, ER+, ER−, LumA.
- **Nuisance sensitivity** (hypothetical): τ + d shifts the curve exactly; at METABRIC size an offset up to ≈ 12 months cannot be separated from a factor-2 change in initial burden by event times alone; thresholds 0.01 → 0.5 move F(20 y) from 0.69 to 0.21.

## Clinical feasibility decision (criteria pre-registered before the results)
- **A. Descriptive population-level comparison: YES** (labelled consistency, not validation; 5 and 10 years primary; τ, d, threshold, initial burden as explicit nuisances).
- **B. Restricted exploratory population calibration: NOT SUPPORTED.** The pre-registered criterion (a finite profile interval for some population combination with τ and the threshold/N₀ free) is met only in an idealised one-heterogeneity-source test: m and s stay estimable (m ± 16%) with τ and N₀-scale free, but τ itself spans ≈ 53 months and the N₀ scale is open; with a second varying parameter the location is a ridge of ≈ 1 log unit; τ = 12 months is indistinguishable from an N₀ rescale ×0.5 at METABRIC size. Only a phenomenological event-time description is estimable, which is descriptive (A).
- **C. Mechanistic calibration: UNSUITABLE** (mechanistic parameters not practically identifiable from burden or event data; event times give one scalar; nuisances confounded; METABRIC has event times only, undocumented origin, no treatment dates).
- **Recommendation: stop at descriptive comparison (A) unless the project owner separately approves a restricted B scope under the stated conditions.**

## Remaining uncertainties
- Global (not just local) structural identifiability not proven; practical results depend on the assumed observation model and on 6 hypothetical ground truths with 8 replicates each; multistart fits on ridge-like likelihoods; Wald diagnostics poorly calibrated far from the truth; profile grid step 0.33 in ln.
- METABRIC RFS origin, diagnosis-to-treatment interval, death-first share and detection delay remain unknown; 161 events lack vital status.
- Population-level analyses use one hypothetical population (prior P1) and a single death-lag assumption.

## Proposed Phase 5 implementation plan (no clinical calibration without separate approval)
1. **Descriptive METABRIC analysis (scope A):** cohort-stratified Kaplan–Meier and Aalen–Johansen (where vital status is known) with CIs at 5/10/15/20 years, ER and PAM50 strata passing the reliability criteria, landmark analysis after 5 years; report the unknown RFS origin as a stated limitation.
2. **Consistency envelope, pre-registered:** define in advance a consistency criterion (e.g. whether the METABRIC CI at 5 and 10 years lies inside the envelope of simulated cumulative-incidence curves over the nuisance grid τ ∈ [0, 24] months, d ∈ [0, 12], threshold 0.02–0.3, N₀ ×0.1–10) *before* looking at the comparison; no fitting; report which hypothetical priors are inconsistent with the descriptive curves.
3. **Only if the owner approves scope B:** fit a phenomenological event-time family (location, scale, cure fraction), not mechanistic parameters, with τ free, by interval-censored likelihood; report profile intervals and the equivalence-map limits; stop if the pre-registered B conditions are not met.
4. **Identifiability-driven data requirements** for any future mechanistic work: paired P-like and D-like (or immune) biomarker time series on an absolute scale, visits through the nonlinear regime, an anchor for the time origin; document them as a requirements list, since the current data cannot supply them.
5. Extend the Monte Carlo with a stochastic birth–death–switching model (M6) before any quantitative comparison, since all present variability is parameter variability.
6. Carry forward: truncation at detection, Tobit LOD, solver-mismatch data generation, pooled replicates for rare events, and config-hash provenance.
