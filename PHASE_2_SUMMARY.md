# Phase 2 Summary: Formal Bifurcation Analysis and Mathematical Verification

*2026-11-04. Hypothetical P–D–E model, illustrative parameters. Mathematical results only: no biological or clinical validation is claimed. Full derivations: `report/phase_2_mathematical_analysis.md`.*

## Tests
**71 passed, 0 failed** (`.venv/bin/python -m pytest`): the 16 Phase 1 tests unchanged plus 55 new. Phase 1 output `results/tables/baseline_scenarios.csv` regenerates byte-identical.

## Files created
- `src/dormancy/continuation.py` pseudo-arclength continuation with adaptive steps, event bisection, multi-seed branch tracing
- `src/dormancy/bifurcation.py` Sotomayor conditions, branch-graph scans (folds, Hopf function, bistability), degenerate η*/ρ*
- `src/dormancy/relapse.py` relapse-time estimators A0–A3, k_det
- `src/dormancy/normalized.py` normalised model (M1), gauge transform
- `configs/normalized.yaml`
- `experiments/bifurcation_analysis.py`, `relapse_time_approximation.py`, `normalized_model_comparison.py`, `_style.py`
- `tests/test_equilibria_phase2.py`, `test_stability_phase2.py`, `test_bifurcation_continuation.py`, `test_normalized.py`, `test_relapse.py`
- `report/phase_2_mathematical_analysis.md`, this file
- Tables (`results/tables/`): `bifurcation_branches_baseline/eta0p6.csv`, `bifurcation_events_*.csv`, `transcritical_verification.csv`, `transcritical_asymptotics.csv`, `continuation_convergence.csv`, `equilibrium_enumeration_crosscheck.csv`, `parameter_scan_summary.csv`, `parameter_scan_examples.csv`, `bistability_map_eta_rho.csv`, `bistability_simulations.csv`, `oscillation_search*.csv`, `relapse_time_approximation.csv`, `relapse_time_accuracy_summary.csv`, `relapse_time_attribution.csv`, `normalized_model_equivalence.csv`
- Figures (`results/figures/`): `fig_A` … `fig_H` (A branches, B stability diagram, C eigenvalue, D dormant fraction, E phase portraits, F additional bifurcations, G relapse time, H normalised vs original)

## Parameter documentation
The normalised parameter table (β = ks/μ, e = μE/s, units, β_c = r_eff) is in `report/phase_2_mathematical_analysis.md` §5 and `configs/normalized.yaml`; `parameter_evidence.md` carries a dated status note on the gauge proposal. `configs/baseline.yaml` is unchanged.

## Files modified (additive only)
`src/dormancy/model.py` (+`jacobian_k`, `drhs_dk`, `hessian_form`), `equilibria.py` (+exact branch curve, cubic, `enumerate_equilibria`, `multistart_equilibria`, `fold_P`, `residual`), `stability.py` (+Routh–Hurwitz, Hopf function, labels), `parameters.py` (+read-only properties `beta`, `r_eff`, `phi`), `__init__.py` (exports), `README.md` (Phase 2 usage). No existing function, parameter, config or Phase 1 document was changed.

## Mathematical findings
1. **Exact structure.** All nonzero equilibria lie on the explicit graph k = A(P)B(P)/s, D = qP/(a+δ_D), E = s/B(P). They satisfy a **cubic** in P, so there are **at most 3** nonzero equilibria for any parameters, and P < r_eff K/(rφ). The dormant fraction at every nonzero equilibrium equals q/(q+a+δ_D) = 0.8696 identically (not an outcome).
2. **Transcritical bifurcation at k_c = 0.28333 (baseline) is verified** via the Sotomayor conditions: simple zero eigenvalue (others −0.05, −0.0393), F_k = 0, transversality α = −0.0610, nondegeneracy β = −0.620, centre-manifold branch slope −0.0984467294 equal to the exact branch slope (16 digits), stability exchange at all 12 sampled k ≠ k_c. Qualification: the exchanging branch has P < 0 for k > k_c, so physically it is a **forward (supercritical) branching of a stable low-burden state out of the tumour-free state**; the nonzero state is admissible only for k < k_c.
3. **The direction is parameter dependent.** k′(0) > 0 (backward bifurcation) for η > η* = 0.314 or ρ < ρ* = 0.0041 (others baseline). At η* the transcritical point is degenerate (β = 0).
4. **Saddle-node bifurcations and two kinds of bistability exist** (continuation, exact cubic roots, eigenvalue labels, simulation). (a) *Tumour-free + nonzero*: e.g. η = 0.6 (backward transcritical): fold at k_SN = 0.8608, both stable for k ∈ (0.2833, 0.8608). (b) *Two stable nonzero states* (tumour-free unstable) when a fold pair lies below k_c: e.g. η = 0.12 (folds 0.2293, 0.2302) and η = 0.164 (window (0.2561, 0.2833); only (0.2833, 0.2840) also has stable tumour-free). Fig. F2 distinguishes them.
5. **Two stable nonzero equilibria at one k** (kind (b)) occur in ≈ 18% of the 14,055 sampled valid parameter sets, always in a window below k_c (large ρ/μ): low-burden focus + high-burden node + saddle; confirmed by simulation from 30 random initial conditions.
6. **No Hopf bifurcation, no oscillation found** (0 Hopf candidates in 14,055 parameter sets; 0 persistent oscillations in 263 long simulations; 15 initially flagged cases were slow relaxations). One Hopf candidate occurs only on the nonphysical P < 0 continuation (η = 0.6) and is not reported as a model bifurcation. This is "none in the scanned ranges", not a proof of absence.
7. **Baseline is structurally simple:** one transcritical point, no fold, no Hopf, no bistability; about a factor 2 (in η) away from the onset of bistability.
8. **M1 normalisation:** e = μE/s, β = ks/μ gives ė = μ(1−e) + ρP/(h+P)e − ηNe with k, s eliminated; β_c = r_eff = 0.11333. Equivalence to the original model verified across a 1000-fold range of the E gauge (relapse-time difference ≤ 8.4×10⁻⁷ months, state differences ≤ 10⁻¹¹), and with a negative control.
9. **Relapse-time discrepancy explained** (k = 0.5 k_c, N_det = 0.1; 375.5 vs 237.5 months): projection of the 99%-dormant initial state onto the growing mode +44.0 months; logistic and other no-feedback nonlinearity +11.5; immune feedback during growth (E reaches 1.40 E_f) +82.5 (order-dependent split). Estimator performance: A1 (eigen-projected, closed form) is the best *closed-form* estimate at N_det = 0.1 (10% median, 42% max error) and is accurate to ≤ 3% at N_det = 10⁻³; A2 (second-order manifold, closed form) is good only for N_det ≲ 10⁻² (≤ 3.4% at 10⁻²) and returns "no relapse" at 0.1 (P_det/h = 0.65). **A3 is semi-analytic**: closed form to N = 10⁻³ plus a numerical tail that is ≈ 67–81% of T at N_det = 0.1 (τ = 252.046 months, identical for all initial compositions), so its 0.007% error is not independent evidence of a closed form; its value is the structural result T = T_lin(x₀) + τ(parameters).
10. **Validity boundary is k_det, not k_c** (k_det/k_c = 0.7548 at N_det = 0.1): at baseline, 0 of 54 fine-scan runs with k > 1.001 k_det produced an event and T grows as k → k_det⁻ (≈ 16,300 months at k_det). Verified only at baseline parameters and initial conditions; overshoot of a focus-type equilibrium could shift the boundary elsewhere. A0 and A1 still return finite times for k > k_det.

## Critical parameter values (baseline, illustrative)
| quantity | value |
|---|---|
| k_c (transcritical) | 0.28333 (continuation: 0.283333333333) |
| β_c = r_eff | 0.11333 |
| η* (k′(0) = 0) / ρ* | 0.3138 / 0.00415 |
| first fold appears (η) / bistability onset (η) | ≈ 0.116 / ≈ 0.164 |
| η = 0.6: fold k_SN | 0.86083 |
| k_det/k_c at N_det = 10⁻³, 10⁻², 0.1 | 0.9954, 0.9568, 0.7548 |
| k′(0), dP*/dk at k_c | −10.158, −0.098447 |

## Remaining uncertainties
- All parameters are hypothetical; only the structure transfers to the real problem.
- No global proof (boundedness, absence of limit cycles); scans are bounded and sampled.
- Transcritical label is for the extended (P ∈ (−h, ∞)) system; only quadratic normal-form coefficients were computed.
- Event-location accuracy at steep branch points is 10⁻⁹–10⁻⁶ in k (analytic values used for conclusions).
- A2/A3 validity thresholds measured at baseline; they scale with λ₊/|γ| and φh. A3 is semi-analytic. N₀ was not varied. The no-event-beyond-k_det statement was checked only at baseline.
- M1 fixes the unit of E but does not make ρ, h, η identifiable.
- Phase 1 documents were left untouched; where their wording is baseline-specific (e.g. no unstable interior state) §8 of the report records the exception.

## Recommendations for Phase 3
1. Local/global sensitivity (LHS + rank correlation) of **relapse occurrence and time**, using β/β_c, k_det, φ, λ₊ and N_det/φh as natural, scale-free coordinates; the structure shows relapse hinges on a few combinations.
2. Define the Monte Carlo output as P(event before horizon) *and* T; account for the censoring induced by k_det (events never occur for β beyond it).
3. Place priors in the normalised model (β, not k, s) to avoid the gauge.
4. Decide explicitly whether the multistable regimes (large η, small ρ, or large ρ/μ) are in the prior support; in them the outcome depends on the initial burden, not only on parameters.
5. Use the **numerical first-crossing time** as the production event time. Keep A1 as an interpretable closed-form check (accurate for N_det ≲ 10⁻², within ~10% at 0.1 only in the median) and A3 as a diagnostic of the structure T = T_lin(x₀) + τ.
6. Defer clinical comparison, Streamlit and parameter fitting as planned.
