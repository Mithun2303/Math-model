# Research Gap, Contribution, Objectives and Scope

*Phase 1 deliverable, 2026-11-04. Read together with `literature_review.md`, `parameter_evidence.md` and `data/dataset_feasibility.md`.*

## 1. Model under evaluation (frozen in Phase 1)

P (proliferating), D (dormant), E (effective immune activity):

- dP/dt = rP(1 − (P+D)/K) − qP + aD − kEP
- dD/dt = qP − (a+δ_D)D
- dE/dt = s + ρ·P/(h+P)·E − μE − η(P+D)E

Implemented and tested: analytic Jacobian, tumour-free equilibrium (0, 0, s/μ), nonzero equilibria by reduction to one scalar equation, tumour-free stability boundary r < ks/μ + qδ_D/(a+δ_D), first-upward-crossing relapse event with right-censoring at the horizon.

## 2. Evaluation of the model against the literature

| Question | Assessment | Basis |
|---|---|---|
| Is the P/D separation biologically defensible? | **Defensible as a simplification, not as a measured partition.** Cellular dormancy (G0–G1 arrest) is an established mechanism (`Aguirre2007`, ABS) and state switching is modelled elsewhere (`Blath2025`, `Page2005`, quiescence augmentation in `Wilkie2013Review`). But the binary form is contested: `Baffoe2026` (ABS) argues dormancy is "a continuous and heterogeneous phenotype". | Abstract-level evidence. |
| Is the immune interaction reasonable? | **Reasonable as an aggregate, with two caveats.** (i) Mass-action killing kEP has no saturation; Wilkie & Hahnfeldt use a saturating, ratio-dependent predation term with an extra innate-immunity term (FT). (ii) Our constant source `s` gives an immune-surveilled tumour-free state; their I(0)=0 treatment simulations do not. | Equations compared with `Wilkie2013CancerRes` (FT). |
| Should immune killing act only on P? | **Contested.** `Mehdizadeh2021` (ABS) describes a dormant-like, *more immunogenic* population; `Saudemont2005` (ABS, mouse AML) finds persistent dormant cells resist CTL lysis but are killed by NK cells. Our assumption 6 is one of several literature-consistent options. | Abstract-level. |
| What do the transition parameters represent? | q and a are first-order state-switching rates (mean dwell times 1/q and 1/(a+δ_D)); δ_D is dormant-cell loss. They are **population-averaged mathematical rates**, not measured quantities. | `parameter_evidence.md` |
| Which mechanisms are omitted? | Environmental feedback on carrying capacity or awakening (`Wilkie2013CancerRes` FT; `YanezFeliu2026` ABS); immunoediting/resistance (`Wilkie2013Interface` ABS); angiogenic dormancy (`Aguirre2007`); spatial structure, cell-death-driven repopulation (`Enderling2009`); treatment; stochastic extinction of small populations; time-varying parameters; delays. | Compared to listed sources. |
| Which parameters have support? | **None quantitatively for human post-treatment disease.** Only mouse-model values were retrieved (different units, functional forms). One qualitative constraint: a·dwell-time consistent with "years to decades" latency. | `parameter_evidence.md` |
| Identifiability? | Structural redundancy k↔s (gauge in E); at low burden only ks/μ matters; ρ, h, η unobservable from low-burden data; total-burden-only data cannot identify the P/D split (to be tested). Wilkie & Hahnfeldt (FT) independently show that four parameter sets fit growth data equally well yet predict different treatment outcomes. | Algebra on our equations; `Wilkie2013CancerRes`. |

**An observation from the existing code, flagged for Phase 2 verification.** Below k_c a stable nonzero equilibrium exists with dormant fraction 0.87, and its burden N* → 0 as k → k_c⁻ (`parameter_evidence.md` §5). Consequently "relapse" in the current model needs either N* > N_detect or a long transient escape from the near-unstable tumour-free state. The latter is what produces delayed relapse in the baseline runs (relapse at ≈171 months for k = 0.25 k_c; none within 240 months for k = 0.5 k_c, where the equilibrium burden 0.58 exceeds N_detect so a crossing occurs later, ≈375 months in a longer run). Whether dormancy-then-relapse is *transient dormancy in the sense of Wilkie & Hahnfeldt* (a saddle's slow manifold) or *slow linear growth from a small residual* is an open analytic question.

## 3. Research gap

Based on the sources read (a bounded, non-systematic search):

1. Dormancy–immune ODE models are analysed mainly for **qualitative regimes** (elimination, dormancy, escape) and **parameter-window bifurcations**. We found no abstract-level evidence that they are turned into **model-based time-to-event distributions** and compared with **right-censored recurrence data**.
2. Clinical recurrence work (`Rueda2019`, `Pan2017`) characterises late hazards statistically; it does not propose a **mechanistic generative model** of the dormant compartment.
3. **Identifiability and endpoint-definition mismatch** (competing risks, time origin, composite events) is seldom confronted explicitly when linking mechanistic models to survival data.

## 4. Scientific contribution (what this project can legitimately claim)

1. An **interpretable P–D–E framework** with a closed-form tumour-free stability boundary that separates the roles of immune control and the dormancy-loss mechanism, plus numerically verified equilibria and a documented structural-identifiability analysis.
2. A **protocol for converting a deterministic dormancy ODE into a relapse-time distribution** (first-crossing events, explicit parameter uncertainty, censoring at the horizon), reported as *simulated event probabilities*.
3. A **transparent group-level comparison** with a public breast-cancer cohort under a verified endpoint definition, with the mismatches (time origin, composite event, competing risks) quantified rather than hidden.
4. A set of **negative or limiting results** as first-class outputs, for example if the data cannot constrain dormant-cell parameters.

*Not claimed:* clinical validity, patient-level prediction, identification of dormant-cell counts, or novelty of the quiescent compartment itself (Wilkie 2013 review, Page & Uhr 2005, Blath 2025, Baffoe 2026 all include such states).

## 5. Objectives (refined from the project plan)

| # | Objective | Success criterion | Phase |
|---|---|---|---|
| O1 | Verify and extend the equilibrium/stability analysis | Analytic boundary matches eigenvalues (done); locate and classify any second boundary or interior saddle; document each bifurcation with continuation | 2 |
| O2 | Define operational dormancy and relapse | Thresholds documented as illustrative (N_low, D/N, |dlogN/dt|, N_detect) with sensitivity of conclusions to them | 2–3 |
| O3 | Quantify parameter influence | Local and global (LHS + rank correlation) sensitivity of relapse occurrence and time | 3 |
| O4 | Relapse-time distributions | Cumulative simulated event probability at 5/10/20 years with Monte Carlo intervals under at least two prior scenarios | 3–4 |
| O5 | Structural identifiability | Synthetic recovery with only N observed; document which combinations are recoverable (e.g. ks/μ) | 4 |
| O6 | Group-level clinical comparison | Compare model-generated event-time curves with METABRIC RFS Kaplan–Meier curves at 5, 10 (and cautiously 15) years, with endpoint caveats stated | 5 |
| O7 | Ablation | Logistic vs P/D vs P/D/E (and optionally +M) compared on the same targets | 4–5 |

## 6. Research questions and falsifiable hypotheses

| RQ | Hypothesis (can be refuted) | Test |
|---|---|---|
| RQ1 What bounds immune control of residual disease in the model? | H1: the tumour-free state is locally stable iff r < ks/μ + qδ_D/(a+δ_D) | Done at the eigenvalue level; extend to nonlinear basins |
| RQ2 Does dormancy alone prevent escape? | H2: if δ_D = 0 the boundary reduces to r < ks/μ, so dormancy switching alone does not change tumour-free stability | Analytic (done in algebra); verify numerically across (q, a) |
| RQ3 How is relapse time controlled? | H3: for k < k_c the relapse time is dominated by the leading unstable eigenvalue of the tumour-free Jacobian and by the residual burden, approximately ln(N_detect/N(0))/λ₊. | Quantify the discrepancy. A first check at k = 0.5 k_c gives λ₊ ≈ 0.029 per month and a crude estimate of ≈ 240 months vs ≈ 375 months simulated; H3 in this crude form is **not supported** and needs the eigenvector-projection and nonlinear corrections. |
| RQ4 Can relapse likelihood be estimated? | H4: Monte Carlo over parameter priors yields a cumulative event curve F̂(t) whose shape depends strongly on the prior (e.g. 5-year probability varies by more than a factor of 2 across two reasonable priors) | Run two prior scenarios; report spread |
| RQ5 Can the model reproduce group-level recurrence patterns? | H5: a prior can be found such that F̂(t) lies within the METABRIC RFS Kaplan–Meier 95% band at 5 and 10 years; this is not evidence of mechanism | Feasible only with the conditional GO conditions in `data/dataset_feasibility.md` |
| RQ6 Is the P/D split identifiable from N alone? | H6: no; at least two distinct (q, a, δ_D) triples fit the same N(t) within noise | Synthetic recovery |

(The H3 check above uses the existing eigenvalue and simulation results and shows the naive hypothesis fails. That is a useful finding, recorded rather than concealed.)

## 7. Scope

**In scope:** three-state ODE with constant parameters; post-treatment residual-disease initial conditions; deterministic dynamics plus Monte Carlo parameter uncertainty; group-level comparison with the METABRIC RFS endpoint.

**Out of scope for the baseline:** patient-level prediction or any clinical use; treatment scheduling; spatial structure; stochastic extinction (planned extension); parameter estimation of individual patients; combining datasets with different endpoints.

## 8. Proposed model improvements (recorded only; none implemented)

| # | Proposal | Motivation (source) | Priority |
|---|---|---|---|
| M1 | Fix the E-gauge (e.g. s/μ = 1) | Structural redundancy between k and s | High, cheap |
| M2 | Saturating killing: replace kEP by kEP/(1 + P/P_sat) or a ratio-dependent form | `Wilkie2013CancerRes` (FT) | Medium |
| M3 | Immune action on dormant cells with separate coefficient k_D | `Mehdizadeh2021`, `Saudemont2005` (ABS) | Medium |
| M4 | Environmental signal M(t) driving a(M) | `YanezFeliu2026` (ABS); plan §19 | Medium (planned) |
| M5 | Time-varying immune efficacy or supply | `Wilkie2013Interface` (ABS) | Low |
| M6 | Stochastic birth–death–switching version for small populations | `Blath2025` (ABS) | Medium |
| M7 | Continuous dormancy depth | `Baffoe2026` (ABS) | Low |
| M8 | Treatment-start offset τ as an explicit nuisance parameter when comparing with RFS from diagnosis | Time-origin mismatch (`data/dataset_feasibility.md`) | High for the clinical track |
