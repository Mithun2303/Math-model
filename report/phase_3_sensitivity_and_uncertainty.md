# Phase 3: Global Sensitivity, Basins of Attraction and Uncertainty-Based Relapse Simulation

*Everything in this report is a statement about a **hypothetical** mathematical model with **illustrative** parameter ranges and **assumed** uncertainty distributions. "Relapse" means first upward crossing of a model-defined burden threshold. No simulated probability is a patient risk, a clinical prediction, or a calibrated estimate, and none of the thresholds is a clinical imaging limit.*

Reproduce: `.venv/bin/python experiments/{local_sensitivity,global_sensitivity,basin_analysis,monte_carlo_relapse,threshold_robustness}.py` (the last needs the Monte Carlo run first). Design in `configs/phase3.yaml` (sha256[:12] `c7e70720ebfc`): every Phase 3 figure footer and every Phase 3 table (including the per-simulation dumps) written by the experiment scripts carries this hash in a `# phase3 config sha256[:12]=…` header line or footer; the Sobol cache stores it in its key. Tables: `results/tables/p3_*`. Figures: `results/figures/fig_p3_*`.

---

## 1. Design principles

1. **Gauge.** All Phase 3 runs use the normalised model of Phase 2 (`NormalizedParams.to_params()`, s = μ), so the tumour-free immune level is E = e = 1 and killing is the gauge-invariant β = ks/μ. Initial immune level e₀ = 1 unless a basin experiment varies it. Fixed: K = 1, μ = 0.05, h = 0.02.
2. **One trajectory runner.** `uncertainty.run_trajectory` integrates once, densely and non-terminally on [0, 240] (DOP853, rtol 10⁻⁹, atol 10⁻¹²). From the dense path it derives the first upward crossing at every threshold, N(60/120/240), downward crossings after the first crossing, the dormancy-duration metric and solver diagnostics. It is **cross-checked against the verified terminal-event solver `dormancy.solver.simulate` on every Monte Carlo sample**: 5,000/5,000 agree, max |ΔT| = 8×10⁻¹¹ months. The legacy solver is reused, not replaced.
3. **Long-horizon outcome is separate.** `uncertainty.long_run` records the attractor reached and whether the trajectory ever reaches the threshold, with T_long = max(10⁵, 50 / slowest stable-equilibrium rate) months, never mixed with the 240-month outcome. The RHS is clipped to the nonnegative orthant **only in this classifier** (to keep the integrator away from the pole of P/(h+P) at P = −h; the true flow leaves the orthant invariant). A label is accepted only if the end state lies within 10⁻³ (relative) of an enumerated stable equilibrium with residual < 10⁻⁸; otherwise it is `unresolved` (count: 0 in every experiment).
4. **Samples are drawn up-front** from a seed (`numpy.random.default_rng`, scipy `qmc.LatinHypercube` / scrambled `qmc.Sobol`) and mapped over a process pool using top-level functions and plain dicts, so results are independent of scheduling (tested: parallel = serial, exactly).
5. **Start-above-threshold trap.** The legacy `simulate` returns `event = 1, time = 0` when N₀ ≥ N_detect. In Phase 3 such samples have **no defined event** (NaN, excluded and counted), and the supports are chosen so that N₀ ≤ 5×10⁻³ < every threshold ≥ 0.01. Count of excluded samples in all runs: 0 (tested).

### Event definition and exceedance types
Relapse at threshold θ: N = P + D crosses θ upward for the first time, starting below θ; right-censored at the horizon (240 months default; 60/120/240 = 5/10/20 years reported). After a first crossing, exceedance is **transient** if N crosses θ downward within the horizon, **persistent** otherwise; a trajectory that later settles at an attractor below θ is reported separately via the long-run attractor (`N_inf < θ`). Crossings are bracketed on a 0.05-month grid and refined by Brent's method on the dense interpolant (an up-and-down excursion shorter than 0.05 month would be missed).

### Operational "dormancy duration" (hypothetical definition)
Total time before the first crossing (or the horizon) during which N < θ, |d ln N/dt| < 0.01 month⁻¹ **and** D/N ≥ 0.5. At the baseline reference this is **18.1 months** and depends on the thresholds of the definition itself (`p3_dormancy_definition_sensitivity.csv`): 10.2 months for |d ln N/dt| < 0.005, 18.1 for 0.01, 32.9 (32.05 for D/N ≥ 0.7) for 0.02; the dormant-fraction cut-off hardly matters. The quantity is a descriptor of the transient, not a measurement; where the trajectory is right-censored (no exit by 240 months), the duration is not a number and derivatives are reported as undefined.

---

## 2. Local sensitivity (Task 1)

**Method** (`sensitivity.log_derivative`, `flip_margin`; `experiments/local_sensitivity.py`). Central differences in log-parameter at steps h = 5% and h/2 = 2.5%; the finer estimate is reported **only if** both agree within 20% and every stencil point defines the output; otherwise the entry is `undefined` (e.g. no event at a stencil point, attractor switches inside the stencil, right-censored duration) or `unstable` (near a switching boundary). Outputs on a log scale give elasticities d ln(out)/d ln θ; λ is reported as dλ/d ln θ [1/month]. **Event indicators are never differentiated**: for "event by 60/120/240 months" the smallest multiplicative change of each input (either direction, within a factor 10) that flips the indicator is located by scan plus bisection. Inputs: r, q, a, δ_D, β, ρ, η, P₀ (D₀ fixed), D₀ (P₀ fixed), N_detect. Reference points: R_early (β = 0.020, first crossing 149 months), **R_baseline (β = 0.040, 219.9 months, close to the horizon)**, R_late (β = 0.0567, no event by 240 months, ~375 months). All other inputs at baseline. Stencil outcomes: 167 ok, 48 zero, 22 undefined, 3 unstable.

**Results at the baseline reference** (elasticities; `p3_local_sensitivity.csv`; Fig p3_01 tornado, p3_01b map, p3_01c flip margins):

| input | d ln N(120) | d ln N(240) | d ln T (first crossing) | dλ_TF/d ln θ [1/mo] |
|---|---|---|---|---|
| r | +10.3 | +9.9 | **undefined** (no event at r·e^{−5%}) | +0.111 |
| β | −3.46 | −4.07 | +0.99 | −0.037 |
| q | −2.72 | −2.47 | +0.66 | −0.033 |
| ρ | −0.03 | −1.07 | +0.18 | 0 |
| a | +0.99 | +0.54 | −0.15 | +0.003 |
| η | 0.00 | +0.21 | −0.03 | 0 |
| δ_D | −0.03 | −0.03 | +0.007 | 0 |
| D₀ | +0.88 | +0.31 | −0.09 | 0 |
| P₀ | +0.09 | +0.03 | −0.009 | 0 |
| N_detect | **0 (exact)** | **0 (exact)** | +0.23 | 0 |

Observations. (i) The ranking by |elasticity| of N(240) is r ≫ β > q > ρ > a > D₀ > η > P₀ ≈ δ_D; the burden is exponentially sensitive to the net growth rate because 240 months is many e-folding times of λ₊. (ii) D₀ matters ≈10× more than P₀ because the initial state is 99% dormant: N(t) is proportional to the projection on the growing mode (Phase 2, ℓ = (0.86, 0.12, 0)), and D₀ dominates it; total initial burden, not composition, is what a global scan sees (§3). (iii) N_detect cannot change the trajectory and is reported as an exact zero (value 0 in the map, zero-length bar in the tornado), not omitted; it acts only through event time and indicators. (iv) The equilibrium burden N* is independent of P₀, D₀ (monostable here) with elasticities to β of −0.39 and to η of +0.30.

**Flip margins at the baseline reference, event by 240 months** (factor by which the input must be multiplied; `inf` = no flip within ×10): β ×1.084 up, r ×0.967 down, q ×1.128 up, ρ ×1.33 up, N_detect ×1.415 up, a ×0.523 down, D₀ ×0.355 down; η, δ_D, P₀: none. **A 3.3% decrease in r removes the event by 240 months**: the binary outcome at this reference is fragile, and the event-time derivative with respect to r is correctly reported as undefined instead of extrapolated. R_early is more robust (β needs ×2.17, r ×0.78); R_late acquires an event with β ×0.77 or r ×1.13.

---

## 3. Global sensitivity (Task 2)

**Design.** Ten inputs on explicit hypothetical ranges (roughly baseline ×[1/5, 5] except initial burden; `configs/phase3.yaml`): r ∈ [0.03, 0.30], q ∈ [0.01, 0.16], a ∈ [0.001, 0.025], δ_D ∈ [2×10⁻⁴, 0.01], β ∈ [0.01, 0.20], ρ ∈ [0.01, 0.16], η ∈ [0.02, 0.32] (all log-uniform), N₀ ∈ [10⁻⁶, 5×10⁻³] (log-uniform), f_D = D₀/N₀ ∈ [0.5, 0.999] (logit-uniform), N_detect ∈ [0.02, 0.30] (log-uniform). *Check required for meaningful indices:* the event indicator has P(event by 240) = **0.304** over the design (0.300–0.306 for n = 250–4000), inside [0.2, 0.8], so its variance is not degenerate. Regime composition (n = 4000): NZ only 63.5%, TF only 32.4%, NZ+NZ 2.8%, TF+NZ 1.4%. Samples with r_eff ≤ 0 (no growth possible) are a regime (TF only), not failures; solver failures: 0.

**Outputs.** log₁₀N(240); restricted time min(T, 240); event by 240 months; λ_TF; highest stable equilibrium burden. Sobol indices are **not** computed on T conditional on an event (the sample set would change with the parameters); conditional T gets rank correlations among events only, labelled as such.

**Estimators.** Spearman and PRCC (inverse of the rank-correlation matrix) with 400-resample bootstrap intervals, n = 4000 LHS; Morris (40 trajectories, 4 levels, 440 runs); Sobol with the Saltelli (2010) first-order, Jansen total-order and Saltelli (2002) second-order estimators, N = 2048 base rows from one scrambled Sobol sequence (45,056 runs), 200-resample bootstrap. **Validation on the Ishigami function** (a unit test): S₁ = (0.314, 0.444, 0.000), S_T = (0.557, 0.442, 0.244), S₁₃ = 0.245, S₁₂ ≈ S₂₃ ≈ 0, matching the analytic values.

**Results** (`p3_global_*`; Figs p3_02, 02b, 02c, 02d, 03):

| input | log₁₀N(240): Spearman / PRCC | min(T,240): Spearman / PRCC | min(T,240): Sobol S₁ / S_T | event: Sobol S₁ / S_T |
|---|---|---|---|---|
| r | +0.67 / +0.83 | −0.66 / −0.73 | 0.53 / **0.78** | 0.42 / **0.78** |
| β | −0.42 / −0.69 | +0.33 / +0.49 | 0.12 / **0.33** | 0.13 / **0.43** |
| q | −0.18 / −0.38 | +0.17 / +0.28 | 0.01 / 0.12 | 0.02 / 0.20 |
| ρ | −0.08 / −0.19 | +0.11 / +0.18 | 0.03 / 0.08 | 0.02 / 0.16 |
| N₀ | +0.32 / +0.58 | −0.09 / −0.15 | 0.01 / 0.06 | 0.01 / 0.10 |
| N_detect | 0.00 / −0.02 | +0.09 / +0.17 | 0.02 / 0.08 | 0.01 / 0.11 |
| a | +0.05 / +0.11 | −0.05 / −0.08 | ≈0 / 0.01 | ≈0 / 0.05 |
| η, δ_D, f_D | |ρ_s| ≤ 0.06 | ≤ 0.03 | S_T ≤ 0.01 | S_T ≤ 0.03 |

(Sobol S₁ intervals are wide, e.g. S₁(r) for min(T,240) is 0.53 with 95% interval [0.39, 0.69]; S_T intervals are narrower, [0.73, 0.84].) The event indicator's AUC for r is 0.90 (β: 0.29, i.e. low β ⇒ event). Among events only, the conditional crossing time correlates with r (−0.61), N₀ (−0.35), ρ (−0.14), and not with β (+0.05, interval contains 0): conditioning on an event removes much of β's role.

**Interpretation.**
* Over these ranges, simulated relapse timing is controlled by **net proliferation (r) and immune killing (β)**, then by dormancy entry q, immune stimulation ρ, the detection threshold, and total initial burden N₀. Reactivation a, dormant loss δ_D, immune suppression η and the *composition* f_D are close to non-influential for burden and timing (S_T ≤ 0.05).
* **Interactions are large**: S_T − S₁ ≈ 0.25 for r and 0.21 for β (restricted time), 0.36 and 0.30 for the event. Second-order indices (Fig p3_03; cells whose 95% bootstrap interval contains 0 are blank, to avoid presenting estimator noise) resolve **one pair only: r×β for the event indicator (S_rβ ≈ 0.10)**. For log₁₀N(240) and min(T, 240) **no individual pair is resolved** although S_T − S₁ is large for r and β: the interaction variance is spread over many pairs and/or higher orders, so the claim "r and β interact" is supported by S_T − S₁ and by the single resolved pair for the event, not by a dominant two-factor term. (An early unmasked heatmap showed a uniform 0.04 row for q; that was estimator noise inherited from the wide S₁(q) interval and is not reported.)
* Spearman and PRCC agree in sign and in the top-3 ordering with Sobol; PRCC magnitudes exceed Spearman's (partial out the dominant r). **PRCC assumes monotone responses**, and the system has saturation and bistability. The binned scatterplots (Fig p3_02b) show what a rank statistic can and cannot summarise: for log₁₀N(240) the binned medians increase with r, decrease with β, decrease weakly with q and are nearly flat in ρ (Spearman −0.08), with very wide 10–90% bands; for min(T, 240) the median is pinned at the 240-month censoring limit except at large r, so only the lower part of the distribution carries information. The relations are noisy and thresholded (event / no event), so the rank correlations rank the inputs but do not describe a smooth response.
* **Convergence.** The ordering of the top inputs (r, β, then q/ρ/N₀) is stable from n ≈ 2000 (Kendall τ vs n = 4000: 0.82 at n = 2000, 0.64–0.82 at n = 1000 because the near-zero inputs reorder among themselves); the median Spearman bootstrap half-width falls 0.123 → 0.028 for n = 250 → 4000. Sobol S_T(r) stays within 0.76–0.81 for N = 128…2048 and S_T(β) within 0.27–0.33 from N = 256 (`p3_global_sobol_convergence.csv`: not yet converged to better than ±0.03); S₁ estimates converge slowly (S₁(β): 0.22 → 0.12 from N = 128 to 2048). Runtime was not limiting for the LHS; the 45,056-run Sobol design took 5–16 minutes on the available 8-core machine depending on concurrent load and is cached by config hash.
* **Sensitivity ≠ uncertainty.** These indices are properties of the model over the screening ranges. They say which inputs *would matter most if they were uncertain over those ranges*; they do not say how uncertain any input is (§5, §8).

---

## 4. Basins of attraction (Task 3)

Three hypothetical parameter sets (`experiments/basin_analysis.py`), each *verified by Phase 2 equilibrium enumeration and eigenvalues*:

| case | parameters (normalised, others baseline) | stable states | saddle |
|---|---|---|---|
| **A: tumour-free + nonzero** | η = 0.6, β = 0.2288 (k = 0.572, Phase 2 B2) | TF; NZ at N* = 0.737 | N = 0.1495 |
| **B1: two nonzero, baseline-adjacent** | η = 0.164, β = 0.108 (k = 0.27) | NZ-low N* = 0.0199 (focus); NZ-high 0.502 (node); TF is a *saddle* | N = 0.2085 |
| **B2: two nonzero, well separated** | a = 0.306, q = 1.58×10⁻³, ρ = 0.240, η = 1.637, β = 0.060 | NZ-low 0.0028; NZ-high 0.982; TF saddle | N = 0.111 |

**Method.** For each case two 60×60 initial-condition slices: (P₀, D₀) log-spaced in [10⁻⁶, 0.9] at e₀ = 1, and (N₀, e₀) with N₀ log-spaced, e₀ ∈ [0.05, 1.6] linear, f_D at the equilibrium ratio q/(q+a+δ_D). Each initial state is classified by the attractor reached (long-horizon LSODA, matched to the enumerated stable equilibria) **and** simulated for 240 months for the event outcome at N_detect = 0.1 plus extra thresholds (0.3, 0.5 for A, B1; 0.01, 0.3 for B2). 3,000 random 3-D initial states per case (P₀, D₀ log-uniform, e₀ uniform) give attractor frequencies. Boundaries are located by bisection along lines (bracket 2⁻⁴⁰).

**Numerical verification.** Boundary cells (120–203 per slice) and 200 interior cells per slice were reclassified with tighter numerics (DOP853 10⁻¹²/10⁻¹⁵, Radau 10⁻¹¹/10⁻¹⁴ for the long run, horizons ×5): **label agreement 100% in all six slice-by-case combinations**, 0 unresolved; a 10⁻⁶ relative jitter of the initial state flipped **no** labels, including in boundary cells. (Boundary cells are cells *adjacent* to a label change on the 60×60 grid, not points on the separatrix: the bisection brackets the separatrix itself.)

**Findings** (`p3_basin_*`; Figs p3_04, p3_05, p3_05b).
1. **A locally stable tumour-free state does not guarantee control: its basin is bounded.** In case A the nonzero state attracts 21.6% of the (P₀, D₀) slice and 13.5% of the (N₀, e₀) slice, **all with N₀ ≳ 0.15** (bisection: boundary at N₀ = 0.165 at e₀ = 1 and 0.145 at e₀ = 0.3; the saddle burden is 0.1495; the saddle's own immune level and composition differ from the slice, so the boundary need not equal the saddle burden). So the nonzero basin is reached only from initial burdens well above any plausible *residual-disease* start; every initial state below N_detect = 0.1 (n = 2,475 on the slice, 100% of those) ends at TF, but 14 of 2,071 random below-threshold starts register an event at θ = 0.1 and fall back (transient exceedance, §5). With a higher threshold (0.5), 290 of 389 below-threshold starts in the NZ basin cross it: the outcome of the relapse definition depends on where the threshold lies relative to the saddle.
2. **Two-nonzero case B1 (baseline-adjacent):** NZ-high is reached only from N₀ ≳ 0.2 (bisection 0.216 at e₀ = 1, 0.181 at e₀ = 0.3; saddle 0.2085): 21.1% / 11.5% of the two slices. At θ = 0.1, 51 of 2,071 random below-threshold starts register an event although they belong to the NZ-low basin (attractor 0.0199 < θ): 19 are followed by a downward crossing within 240 months, and the other 32 are still above θ at 240 months (overshoot that has not yet relaxed). Immune depletion lowers the boundary slightly but not enough to matter for small residual starts.
3. **Case B2 (strong immune stimulation) shows basin selection from below-threshold starts.** NZ-high attracts 35.5% / 20.6% of the slices; **8.5% of random 3-D starts with N₀ < 0.1 reach the high attractor and all of them (177/177) register an event at θ = 0.1** (at θ = 0.3, 585/585 starts below 0.3 in the high basin do). Depletion of the immune level lowers the boundary below the saddle burden: bisection gives N₀ = 0.071 at e₀ = 1 and **0.039 at e₀ = 0.3** (saddle 0.111). The NZ-low attractor (N* = 0.0028) is below θ = 0.01 too, yet 262 of 1,325 starts below 0.01 cross 0.01 and then fall back (all transient).
4. **Event outcome ≠ attractor.** The 240-month event and the attractor differ in both directions: starts in the *low* basin can cross the threshold transiently (B2: 8 of 1,895 at θ = 0.1; 19 of 2,071 in B1), and starts in the *high* basin counted as events only if they start below θ (otherwise undefined: all but one of case A's NZ-basin random starts have N₀ > 0.1). Figure p3_05b shows event fraction by attractor for starts strictly below each threshold; groups with n < 10 are omitted.
5. A transient-exceedance **fixture** (used as a unit test) found numerically: B2 parameters, N₀ = 0.0328, f_D = 0.493, e₀ = 0.606, θ = 0.05: first crossing at 8.75 months, downward crossing at 18.68 months, N(240) = 0.0022.

---

## 5. Monte Carlo relapse simulation (Tasks 4, 5, 8)

### Uncertainty distributions (hypothetical; `configs/phase3.yaml`)
**P1_base (primary).** r = 0.12, μ, h fixed; log-normal with median and 95% factor f (σ = ln f / 1.96): q (0.04, f = 3), a (0.005, 5), δ_D (0.001, 5), β (0.04, 3), ρ (0.04, 3), η (0.08, 3); N₀ log-normal (10⁻⁴, 30) truncated to [10⁻⁷, 5×10⁻³]; f_D ~ Beta(10, 1) (mean 0.91). Detection threshold 0.1, e₀ = 1. **Alternatives (model-assumption scenarios):** P2 wide (σ×1.5, f_D ~ Beta(4,1)); P3 log-uniform over the global ranges; P4/P5 Gaussian-copula rank correlation ±0.5 on (q, a) and (β, ρ) (induced Spearman ρ_s = 0.48 for +0.5, tested); P6 r also uncertain (log-normal, median 0.12, f = 2); P7 initial state fixed at baseline (parameter uncertainty only). The **regime composition under P1** (a restatement of the prior): TF unstable (β < β_c) in 95.9% of samples; NZ only 94.4%, TF only 3.8%, NZ+NZ 1.6%, TF+NZ 0.2%.

### Estimands (all simulated, N_detect = 0.1)
P(first crossing by 60/120/240 months) with 95% Wilson (and Clopper–Pearson) intervals; one-sided 95% upper bound for rare counts; empirical cumulative distribution (Kaplan–Meier equals the empirical CDF because the horizon is the only censoring); restricted mean event-free time RMST(τ) = mean(min(T, τ)) exactly; conditional median among simulations that cross. **The unconditional median event time is not reported when P(event by 240) < 0.5** (it is `NaN`; here 0.424).

### Results for P1_base (`p3_mc_P1_summary.csv`, `p3_mc_P1_pooled.csv`; Figs p3_06, 07, 08, 10, 12)

Three **independent runs of the same prior** (different seeds) are available: run A (seed `mc_base`, n = 5000; its first 1000 samples are the "initial run"), run B (seed `mc_assumptions`, n = 5000) and run C (five batches of 1000, seeds `mc_base`+1000…+1004). **The headline uses all three pooled (n = 15,000)**, with Clopper–Pearson intervals, because a single run is too small for the 5-year event: run A has 0 crossings by 60 months, runs B and C each have 4.

| quantity (N_detect = 0.1) | run A (n = 5000) | run B (n = 5000) | run C (n = 5000) | **pooled (n = 15,000)** |
|---|---|---|---|---|
| P(by **60 mo / 5 y**) | 0/5000 | 4/5000 | 4/5000 | **8/15,000 = 0.00053 [0.00023, 0.00105]** |
| P(by **120 mo / 10 y**) | 0.0854 | 0.0966 | 0.0910 | **0.0910 [0.0864, 0.0957]** |
| P(by **240 mo / 20 y**) | 0.4240 | 0.4418 | 0.4432 | **0.4363 [0.4284, 0.4443]** |
| earliest first crossing | 61.6 mo | 55.4 mo | 51.9 mo | **51.9 mo** |
| conditional median T among events | 153.2 | 154.1 | 153.0 | **153.5 mo** (IQR 124.8–186.7; 6,545 events) |
| RMST(240) | 204.2 | 203.3 | 202.9 | **203.5 mo** |
| unconditional median | not defined | not defined | not defined | **not defined (P < 0.5)** |

Homogeneity of the three runs (χ² test on counts): p = 0.15 (120 mo), p = 0.10 (240 mo); the 60-month counts (0, 4, 4) are compatible with a single rate of 0.00053 (probability of 0 events in 5000 at that rate ≈ 7%). **The one-sided upper bound 0.0006 obtainable from run A alone is therefore not a reliable bound for the 5-year probability**, and run A's "0/5000" should not be quoted by itself. The "initial run" (first 1000 of run A) gives P(60) 0/1000, P(120) 0.081 [0.066, 0.100], P(240) 0.423 [0.393, 0.454]. Run A is also the run with long-horizon and regime information (per-sample file `p3_mc_P1_samples.csv`); runs B and C are in `p3_mc_P1_base_run_B_samples.csv` and `p3_mc_P1_batches_run_C_samples.csv`, and every prior scenario has its own per-simulation file (`p3_mc_P2_wide_samples.csv` …), each with an explicit `event_indicator`, `time_or_censor` (240 where there is no event) and `censored_at_horizon` column. No solver failures and no start-above samples in any run.

The remaining statements in this subsection use run A (n = 5000). Of its 2,120 events, 29 (1.4%) are transient within 240 months (cross, then cross back) and 98.6% are still above the threshold at 240 months; 28 events (1.3%) end at an attractor below the threshold. Event probability by regime (P1): NZ only 0.446 (n = 4718), NZ+NZ 0.177 (n = 79), TF only 0.000 (n = 192), TF+NZ 0.000 (n = 11); by the attractor *reached*: NZ 0.446, NZ-low 0.145 (n = 76; 36% of these events transient), TF 0.000 (n = 203). The long-horizon event probability is 0.799 vs 0.424 by 240 months: 37 percentage points of the "eventual" relapses are censored by the 20-year horizon.

**Persistent-exceedance definition (sensitivity analysis, `p3_mc_persistent_definition.csv`).** If an event is counted only when N never crosses back below the threshold within 240 months, P(by 120 mo) and P(by 240 mo) change from 0.0966 to 0.0942 and from 0.4418 to 0.4370 for P1 (run B); the largest effects among the seven priors are −0.013 at 240 months (P3 log-uniform: 0.3168 → 0.3038; P2 wide 0.3808 → 0.3694) and −0.008 at 120 months (P2 0.1510 → 0.1436; P3 0.1198 → 0.1120). The fraction of first crossings that are transient is 0.3% (P4) to 4.1% (P3). The 5-year probability changes by at most 0.002 (P6: 0.0362 → 0.0340). The conclusions do not depend on this choice.

### Convergence and numerical verification
* **Sample size.** Running P(by 240) with 95% Wilson band (run A prefixes) stabilises by n ≈ 2000 (±0.014) and is 0.424 ± 0.014 at n = 5000; the pooled n = 15,000 estimate has half-width 0.008. The five n = 1000 batches give 0.442, 0.460, 0.423, 0.443, 0.448 (binomial SE at n = 1000 is 0.016; at n = 5000 0.007). The 120-month probability: batches 0.079–0.099. **The 60-month probability is a rare event** (pooled 8/15,000): individual 5000-runs give 0, 4 and 4 crossings, so a single run cannot give a trustworthy bound; report pooled counts with Clopper–Pearson intervals, not a ratio of tiny probabilities.
* **Seeds.** Runs A (0.424), B (0.442) and C (0.443) of the same prior differ at the second decimal (χ² homogeneity p = 0.10): seeds matter at the level of the interval width, so single-run values should be quoted with that in mind (A vs B: 1.8 standard errors of a difference).
* **Solver tolerances** (600 samples, five settings: rtol 10⁻⁶ … 10⁻¹², Radau, LSODA): **no event flips**; max |ΔT| = 0.0067 months at the loosest setting (1.2×10⁻⁵ at the tightest).
* **Reproducibility:** same seed ⇒ identical events and times (max |ΔT| = 0); parallel (8 workers) = serial (150 samples), identical.
* **Numerical failures:** 0 in 5,000 (P1) and in all seven scenario runs; no negative states (min state ≥ −10⁻¹⁰ in the tested random samples).
* **Right-censored simulations:** 2,880 of 5,000 are censored at 240 months in P1; the long-horizon run shows that 1,873 (65%) of them reach the threshold later and 1,007 never do (so event-free at 240 months is not the same as relapse-free).

### Parameter uncertainty vs event-time uncertainty (Fig p3_11)
Sobol indices **under the P1 prior** (propagated uncertainty, N = 1024 base rows, 10,240 runs, first and total order) for min(T, 240): β S_T = 0.63 (S₁ 0.38), q 0.34 (0.20), ρ 0.18 (0.03), N₀ 0.14 (0.08), a 0.05; η, δ_D, f_D ≈ 0. For the event indicator: β 0.72, q 0.39, ρ 0.28, N₀ 0.16, a 0.10. **Under the P1 prior r is fixed**, so β replaces r as the dominant source: the ranking depends on which quantities are given distributions. Within strata of β/β_c the crossing times still spread widely (0.15–0.25: IQR ~117–165 months).

---

## 6. Equilibrium-based predictions vs simulated events (Task 6)

For each P1 sample (n = 5000) the equilibria were **recomputed per sample from the exact cubic** (never from the baseline k_det formula). Predictors against simulated outcomes (`p3_mc_equilibrium_predictions_vs_events.csv`):

| predictor | vs event by 240 mo: accuracy (sens / spec) | vs long-horizon event: accuracy (sens / spec) |
|---|---|---|
| TF unstable (β < β_c; the k_c-type rule) | 0.465 (1.00 / 0.07) | **0.839 (1.00 / 0.20)** |
| highest stable N* ≥ θ | 0.609 (0.99 / 0.33) | 0.983 (0.99 / 0.94) |
| lowest stable N* ≥ θ | 0.618 (0.99 / 0.35) | 0.992 (0.99 / 0.99) |
| attractor *reached* has N_inf ≥ θ | 0.620 (0.99 / 0.35) | **0.994 (0.99 / 1.00)** |

* **k_c is a local stability boundary, not a relapse classifier.** It has perfect sensitivity but specificity 0.20 against the long-horizon event: **804 of 5000 samples (16%)** have an unstable tumour-free state (β < β_c) and yet never reach N = 0.1 because the stable nonzero attractor lies below the threshold (a low-burden "dormant" regime). Against the 240-month event its accuracy is 0.465 (1,872 samples have an attractor above θ but no event within 20 years: finite-horizon censoring).
* The attractor-based predictor agrees with the long-horizon event in 99.4% of samples; its 29 false negatives (0.6%) are trajectories that overshoot or cross transiently and end at an attractor below θ; **0 false positives**. The "event_long ⇒ N_max_long ≥ θ" internal consistency check has 0 violations; 0 attractors unresolved.
* A tumour-free-stable parameter set producing an event: 0 of 5000 in P1 (basin effects occur in the explicitly bistable cases of §4, not in this prior).
* **The apparent detection boundary is sample-specific.** For 300 P1 samples, β was scanned over [0.05, 20]×β (80 values) with a long-horizon run at each; the boundary β_det above which no relapse occurs, relative to the sample's β_c (`p3_threshold_beta_boundary.csv`): at θ = 0.1 median **0.72**, IQR 0.58–0.86, range 0.07–1.04 (simulated), with the equilibrium-based boundary (lowest stable N* = θ) equal to the simulated one in the median (ratio 1.000; IQR 1.000–1.000, range 0.43–26 in the few samples with bistability/multiple roots). No sample had a re-entrant event (a relapse reappearing at larger β). **The Phase 2 baseline value k_det/k_c = 0.7548 (θ = 0.1) is not a universal number** (it is close to the sample median 0.72, but individual samples range from 0.07 to 1.04); it is verified for the baseline at all seven thresholds (event at 0.98 k_det, none at 1.02 k_det, `p3_threshold_baseline_kdet_check.csv`). Maximum values slightly above 1 (1.035) are grid resolution (β grid step ≈ 5%).

---

## 7. Detection-threshold robustness (Task 7; Fig p3_09)

Common random numbers: the same 5000 P1 samples (run A) and dense paths at every threshold (0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5; model-defined, **not imaging limits**). Simulated probabilities (`p3_threshold_event_probabilities.csv`):

| θ | P(60 mo) | P(120 mo) | P(240 mo) | cond. median T | events transient | long-horizon P(reach θ) |
|---|---|---|---|---|---|---|
| 0.01 | 0.040 | 0.352 | 0.689 | 118.6 | 0.2% | 0.949 |
| 0.02 | 0.013 | 0.267 | 0.631 | 129.2 | 0.5% | 0.930 |
| 0.05 | 0.002 | 0.155 | 0.521 | 142.5 | 1.5% | 0.870 |
| **0.1** | **0.000** | **0.085** | **0.424** | **153.2** | 1.4% | 0.799 |
| 0.2 | 0.000 | 0.043 | 0.337 | 166.3 | 0.3% | 0.731 |
| 0.3 | 0.000 | 0.024 | 0.284 | 175.0 | 0.1% | 0.696 |
| 0.5 | 0.000 | 0.008 | 0.214 | 185.6 | 0.0% | 0.615 |

* These are run-A values (n = 5000); at θ = 0.1 the 5-year probability is 0/5000 in this run while the pooled three-run rate is 0.0005 (§5), so P(60 mo) = 0.000 entries for θ ≥ 0.1 mean "below about 1/5000", not "zero".
* The threshold assumption changes P(20 y) by a factor of ~3 (0.689 at θ = 0.01 vs 0.214 at 0.5), P(10 y) by a factor of ~45 (0.352 vs 0.008) and the conditional median crossing time by 67 months (118.6 vs 185.6). It is **as influential as the prior assumptions** (§8), and the 5-year probability is essentially zero for θ ≥ 0.1 and 4% at θ = 0.01.
* **Transient crossings do not change the conclusions**: ≤ 1.5% of events are transient at any threshold, with a maximum near θ = 0.05–0.1, where N* of low-burden samples is comparable to the threshold. Persistent crossing (above θ at 240 months) accounts for the rest.
* **Equilibrium predictions agree with long-horizon events** at every threshold: the attractor-based rule has accuracy 0.994–1.000 vs the long-horizon event, whereas the "TF unstable" rule's accuracy against the 240-month event falls from 0.73 (θ = 0.01) to 0.26 (θ = 0.5).
* **The apparent boundary shifts with θ**: simulated β_det/β_c falls from 0.95 (θ = 0.01) to 0.50 (θ = 0.5), tracking the baseline k_det/k_c formula (0.957 → 0.559) at baseline-like parameters but with a wide inter-sample spread (IQR 0.32–0.67 at θ = 0.5).

---

## 8. Separating the sources of uncertainty (Task 8)

| source | what it is | magnitude for P(240 mo) |
|---|---|---|
| **Monte Carlo sampling error** | finite n for a given model and prior | 95% Wilson half-width 0.014 at n = 5000; 0.030 at n = 1000 |
| **Prior-assumption uncertainty** | seven hypothetical distributions (P1–P7) | range **0.317–0.468** (P3 log-uniform 0.317; P2 wide 0.381; P1 0.442; P5 0.419; P7 0.427; P6 0.452; P4 0.468); max − min = 0.151, **5.5× the median MC interval width (0.027)** |
| **Threshold assumption** | N_detect = 0.01 … 0.5 | 0.689 … 0.214 |
| **Numerical** | solver tolerance, event bracketing | event flips 0/600; |ΔT| ≤ 0.007 mo |

For P(120 mo): range 0.034 (P7: no initial-condition uncertainty) to 0.232 (P6: r uncertain) vs MC width 0.016: ratio 12×; for P(60 mo): 0.000 to 0.036 (ratio 19×; the 5-year probability ranges across priors from 0/5000 (P7) to 181/5000 (P6), while run-to-run variation of one prior is 0–4/5000). **The priors, not the number of simulations, dominate the uncertainty**; this is why these probabilities cannot be read as predictions. Correlation sensitivity (P4 vs P5, ±0.5 on (q, a) and (β, ρ)): P(240) 0.468 vs 0.419, a 0.049 effect, larger than the MC interval width. Fixing the initial condition (P7) moves P(120) from 0.097 to 0.034, so uncertainty in the initial residual burden and composition is what generates early events in P1.

---

## 9. Testing (Task 11)

111 tests pass: the 71 Phase 1/2 tests unchanged plus 40 new (`tests/test_phase3_uncertainty.py`, `tests/test_phase3_sensitivity.py`): reproducible and nested sampling; LHS stratification (each marginal's CDF value falls once in each of n strata, with a power check that plain Monte Carlo fails it); range/family validation; log-normal factor and truncation; copula rank correlation; probability and interval bounds (including k = 0 and k = n); rule-of-three bound; censoring, RMST and median rules; dense-runner vs verified solver; no-event trajectories; **start-above-threshold excluded, not counted**; threshold monotonicity; persistent exceedance; **transient exceedance with a numerically found fixture**; regime classification of the Phase 2 cases; solver-tolerance invariance; nonnegativity over the wide prior; Monte Carlo 1/√n convergence and parallel = serial; Sobol estimators on Ishigami and on an additive model; log-derivative exactness and undefined/unstable flags; event-time derivative not computed across the event boundary; flip margins; spearman/PRCC (confounder removal) and bootstrap; Morris on a linear model; basin classification reproducing the Phase 2 η = 0.6 bistability (small → TF, large → NZ; tight numerics agree); two-nonzero classification; boundary mask.

---

## 10. Limitations and unresolved issues

1. **All ranges and priors are hypothetical** and only loosely anchored (baseline ×[1/5, 5]); the conclusions about *which inputs matter* inherit the ranges (r's range spans a factor of 10, β's 20, N₀'s 5,000). No claim about patient risk.
2. **First-order Sobol intervals are wide**; rankings of the near-zero inputs (a, δ_D, η, f_D) are not resolved, only that they are small. Second-order indices below ±0.03 are noise and are masked.
3. **PRCC/Spearman assume monotone responses**; checked graphically for the top inputs only.
4. **Detection by dense grid** (0.05 month) misses excursions shorter than the grid.
5. **Basins** were mapped on two-dimensional slices (60×60) of the three-dimensional initial-condition space for three parameter sets; fractions of the slice are properties of the chosen log-spaced box and not probabilities of anything; boundaries are not tracked for other parameter combinations. A time-to-event curve of *basin-crossing* initial conditions was not constructed.
6. **Bistable regimes are rare under P1 (1.8%)**, so the prior does not stress basin selection; their event probabilities are estimated from n = 79 + 11 samples.
7. **The long-run classifier clips the RHS** to the orthant (pole of P/(h+P) at P = −h) and assumes convergence to an equilibrium (no limit cycle: supported by Phase 2 scans, not proved).
8. **Model assumptions are untouched:** constant parameters, no treatment, binary dormancy, mass-action killing, deterministic dynamics, no stochastic extinction. Event times here are deterministic first-passage times for each sampled parameter vector, i.e. the Monte Carlo distribution represents *parameter and initial-condition uncertainty*, not demographic stochasticity.
9. **Equilibrium-prediction table uses N_detect = 0.1 for the headline row**; other thresholds are in §7.
10. **The k_det per-sample boundary scan** uses 300 of the 5,000 samples and an 80-point β grid (≈ 5% resolution).

## 11. Files

`configs/phase3.yaml`; `src/dormancy/uncertainty.py`, `sensitivity.py`, `basins.py`; `experiments/_p3common.py`, `local_sensitivity.py`, `global_sensitivity.py`, `basin_analysis.py`, `monte_carlo_relapse.py`, `threshold_robustness.py`; tests as above; tables `results/tables/p3_*.csv`; figures `results/figures/fig_p3_01 … 12*.png` (01 local; 02 global ranking, 02b scatter, 02c convergence, 02d Morris; 03 interactions; 04/05 basins; 05b attractor vs event; 06 time distribution; 07 cumulative probability; 08 5/10/20-year probabilities; 09 thresholds; 10 convergence; 11 parameter vs event-time uncertainty; 12 by regime). Tables begin with a `# phase3 config sha256[:12]=…` provenance comment; read with `pandas.read_csv(..., comment="#")`.
