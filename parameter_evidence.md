# Parameter Evidence Table

*Phase 1 deliverable, 2026-11-04. The equations and baseline values in `src/dormancy/` and `configs/baseline.yaml` are unchanged; this document only records how well each parameter is supported.*

**Model (unchanged):**
dP/dt = rP(1 − (P+D)/K) − qP + aD − kEP; dD/dt = qP − (a+δ_D)D; dE/dt = s + ρ·P/(h+P)·E − μE − η(P+D)E.

## 1. Status vocabulary

| Status | Meaning |
|---|---|
| **Measured** | A value for the same quantity in the same system appears in a source we read. |
| **Estimated (other system)** | Estimated by fitting in a different system (e.g. mouse tumour). Cited for order of magnitude only, not transferable. |
| **Hypothetical** | A mathematical choice with no sourced value. All baseline values currently fall here. |

**Result of this audit: all 11 parameters plus the model-definition constants are Hypothetical for human post-treatment residual disease.** The only numerical parameter values we actually retrieved from a source are from a murine fibrosarcoma model (Wilkie & Hahnfeldt 2013, Tables 1–2), in different units and structure. Values from other papers (Kuznetsov 1994, de Pillis 2005, Page & Uhr 2005) were not retrieved because only their abstracts were read.

## 2. Parameter table

Units: time in months; P, D, E normalised (K = 1; E in arbitrary effector units). "Proposed range" is **our exploratory design choice for sensitivity/Monte Carlo work, not a literature-derived interval**, and must be replaced by sourced priors where possible.

| Symbol | Meaning | Units | Baseline | Evidence actually found | How it could be measured | Proposed exploratory range (design choice) | Status |
|---|---|---|---|---|---|---|---|
| **r** | Proliferation rate of P | month⁻¹ | 0.12 | Wilkie & Hahnfeldt 2013 (CR, FT): murine 3-MCA fibrosarcoma, μ = 0.16 d⁻¹, α = 0.72 (Table 1), i.e. intrinsic rate μα ≈ 0.115 d⁻¹ ≈ 3.5 month⁻¹, about **30× faster** than our r and a different growth law. Not transferable. Our r gives a doubling time ln2/r ≈ 5.8 months. | Tumour volume doubling time from serial imaging; Ki-67/mitotic index; ctDNA growth kinetics | 0.06–0.25 (doubling 2.8–11.6 months) | Hypothetical |
| **K** | Carrying capacity | normalised burden | 1.0 | Wilkie & Hahnfeldt treat K_C as constant or dynamic (dK_C/dt = pC − qK_C C^{2/3}, via Hahnfeldt 1999). Ours is constant and normalised, so absolute scale is unconstrained. | Not directly measurable; maximal tumour burden | Fixed at 1 (units choice) | Hypothetical (units convention) |
| **q** | P → D switching rate | month⁻¹ | 0.04 | None quantitative. Blath et al. 2025 (ABS) model spontaneous and treatment-induced switching into dormancy but numerical values were not read. | Cell-cycle reporter lineage tracing (e.g. label retention) in models; not available in patients | 0.005–0.1 | Hypothetical |
| **a** | D → P reactivation rate | month⁻¹ | 0.005 | Aguirre-Ghiso 2007 (ABS): latency "from years even to decades". With a = 0.005 and δ_D = 0.001, the mean dormant dwell time 1/(a+δ_D) ≈ 167 months ≈ 13.9 years, consistent in order of magnitude with that qualitative statement. This is a plausibility check, not a measurement. | Fraction of label-retaining cells re-entering cycle; indirectly via late-recurrence hazard | 0.001–0.02 | Hypothetical (qualitatively constrained) |
| **δ_D** | Dormant-cell loss rate | month⁻¹ | 0.001 | None found. Note it enters the stability boundary via qδ_D/(a+δ_D); for δ_D = 0 the dormancy term vanishes (§4). With baseline values, 1/6 of dormant-cell exits are losses. | Apoptosis rates of quiescent cells in vitro/in vivo models | 0–0.01 | Hypothetical |
| **k** | Immune killing effectiveness | (E·month)⁻¹ | 0.1 (config); critical k_c = 0.2833 | Wilkie & Hahnfeldt use a saturating ratio-dependent predation Ψ with θ = 2.5, β = 0.5, φ = 50 (Table 1), a different functional form (ours is mass-action kEP). de Pillis et al. 2005 (ABS) fit lysis and tumour-growth data, but values were not read. | Chromium-release/lysis assays, in vivo killing assays | 0.02–0.6 (straddles k_c; see §5) | Hypothetical |
| **s** | Baseline immune replenishment | E·month⁻¹ | 0.02 | None. Wilkie & Hahnfeldt note immune size/growth parameters were "chosen to approximately match those of the tumor" because they are hard to measure. | Circulating antigen-specific T-cell frequency and turnover | Not independent of k (see §3) | Hypothetical |
| **μ** | Immune turnover | month⁻¹ | 0.05 | None retrieved. Mean lifetime 1/μ = 20 months in the baseline; no source checked. | T-cell lifespan/turnover studies | 0.02–0.2 | Hypothetical |
| **ρ** | Tumour-stimulated immune expansion | month⁻¹ | 0.04 | None retrieved; Wilkie & Hahnfeldt use logistic recruitment with λ = 0.22 d⁻¹ and r = 10⁻³ (different form). | T-cell proliferation assays | 0.01–0.2 | Hypothetical |
| **h** | Half-saturation of immune stimulation | burden | 0.02 | None. | Dose–response of T-cell activation vs antigen load | 0.005–0.1 | Hypothetical |
| **η** | Burden-associated immune suppression | (burden·month)⁻¹ | 0.08 | None. | Immunosuppression assays (e.g. Treg/MDSC effects) | 0.01–0.3 | Hypothetical |

### Model-definition constants (not biological parameters)

| Constant | Value | Basis |
|---|---|---|
| N_detect | 0.1 (fraction of K) | Design assumption set during implementation. Wilkie & Hahnfeldt's operational dormancy criteria (tumour diameter 1–30 mm; (1/C)dC/dt < 0.04) show that such thresholds are always operational choices. |
| Horizon H | 240 months | From the project plan; overlaps the 20-year follow-up scale of the late-recurrence literature (Pan 2017; Rueda 2019, ABS). |
| P(0), D(0) | 10⁻⁶, 10⁻⁴ | Hypothetical residual burden; E(0) = s/μ = 0.4. Residual-burden sizes are not measured in the sources read. |

## 3. Structural identifiability observations (derived from our equations, not from literature)

These follow algebraically from the implemented model and are independent of data quality.

1. **Scale redundancy between k and s (and E).** Substituting E′ = cE, s′ = cs, k′ = k/c leaves all three equations unchanged in form. Effector units are therefore arbitrary: k and s cannot be identified separately even with perfect data on P and D. Any reported "immune killing strength" is meaningful only relative to a chosen normalisation of E.
2. **At the tumour-free state only the product ks/μ appears** (via the tumour-subsystem matrix with entry r − q − k s/μ). Early post-treatment dynamics, when P and D are small and E ≈ s/μ, depend on (r, q, a, δ_D, ks/μ). The parameters ρ, h and η affect behaviour only once the tumour is large enough to change E. A relapse-timing study based on the low-burden regime therefore cannot constrain them.
3. **If only total burden N = P + D is observed**, the linearised tumour subsystem has two eigenvalues and two amplitudes, which cannot determine its four matrix entries uniquely. Plan: confirm by a profile-likelihood experiment in Phase 9 rather than assert it here.
4. **Group-level clinical event times** (Kaplan–Meier curves) add only a distribution of threshold-crossing times. They cannot measure P(t), D(t), E(t) (see `data/dataset_feasibility.md`).

## 4. Dormancy-related quantities implied by the baseline (computed with the existing code)

| Quantity | Value | Formula |
|---|---|---|
| Tumour-free immune level | 0.4 | s/μ |
| Critical killing k_c | 0.2833 | (r − qδ_D/(a+δ_D)) · μ/s |
| Mean dormant dwell time | ≈ 167 months | 1/(a+δ_D) |
| Equilibrium D/P ratio (nonzero equilibria) | 6.67 (dormant fraction 0.87) | q/(a+δ_D) |
| Dormant fraction at every nonzero equilibrium found | 0.87 | D*/(P*+D*) |

## 5. Observation from the existing code that affects parameter priors

Using the existing `nonzero_equilibria` routine (no code changed) at baseline parameters with k varied:

| k / k_c | Stable nonzero equilibrium (P*, D*, E*) | N* = P*+D* | Max real eigenvalue at tumour-free state |
|---|---|---|---|
| 0.5 | (0.0756, 0.5042, 0.309) | 0.580 | +0.0290 |
| 0.9 | (0.0034, 0.0226, 0.432) | 0.026 | +0.00225 |
| 0.99 | (2.8×10⁻⁴, 1.9×10⁻³, 0.403) | 0.0022 | +1.8×10⁻⁴ |
| 0.999 | (2.8×10⁻⁵, 1.9×10⁻⁴, 0.400) | 2.1×10⁻⁴ | +2×10⁻⁵ |
| 1.001, 1.05 | none found | n/a | negative |

The nonzero equilibrium is stable, appears only for k < k_c, and shrinks continuously to the tumour-free state as k → k_c⁻. This is the numerical signature of a transcritical bifurcation. It is an **observation to be verified in Phase 2 (bifurcation continuation), not a proven result.** Consequences for priors:

- For k in (k_c·f*, k_c) where N*(k) < N_detect = 0.1, the model predicts persistent low-burden disease that **never crosses the detection threshold**, with 87% of cells dormant. This is the model's tumour-mass dormancy. Relapse in the model therefore requires N* > N_detect or slow transients.
- Wilkie & Hahnfeldt's transcritical bifurcations (in θ) occur at the boundaries of an immune-efficacy window. Ours occurs at one boundary only (k_c). We have not yet looked for a second boundary or for a saddle-type interior equilibrium; their analytic result (an unstable interior dormant state) has no counterpart in our numerical scan so far.

## 6. Priority for obtaining real values

1. **Fix the gauge** (e.g. set s/μ = 1) before any fitting, so k has a defined meaning. *Proposed change; not implemented.* *(Status 2026-11-04: implemented in Phase 2 as a separate normalised model, `src/dormancy/normalized.py` + `configs/normalized.yaml`, with beta = k s/mu; `configs/baseline.yaml` is unchanged.)*
2. Source sourced priors for r (tumour doubling times / Ki-67) and for the aggregate late-recurrence hazard (Pan 2017 full text; METABRIC KM curves in `results/tables/metabric_km_overall.csv`).
3. Treat q, a, δ_D as population-level nuisance parameters: constrain a·(dormant dwell) by the observed *late hazard*, but expect weak identifiability.
4. Read full texts of Kuznetsov 1994 and de Pillis 2005 to extract effector-cell turnover and killing ranges, with explicit species/units caveats.
