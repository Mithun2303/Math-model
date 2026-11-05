# Phase 1 Summary: Literature Review and Dataset Feasibility

*2026-11-04. Existing model, solver, stability code and the 16 tests were not modified.*

## What was investigated
- The five required papers plus 14 further literature sources (all bibliographic records verified against Europe PMC/publisher metadata; see `report/references.bib`).
- Fit of the implemented P–D–E equations to that literature, with a per-parameter evidence audit.
- The METABRIC clinical export (2,509 patients) and, as an alternative, the TCGA-BRCA PanCancer Atlas clinical export (1,084 patients), both downloaded and analysed.

## Verified findings
1. **Access limits.** Only Wilkie & Hahnfeldt 2013 (*Cancer Res*) was read in full. The other four required papers were verified at abstract/metadata level (PMC/PubMed block automated access; Nature paywalls Rueda 2019; the npj 2026 article body was not served). Equation-level claims are made only for the full-text paper.
2. **Literature position.** Quiescent states in immune-tumour models are not novel (Wilkie 2013 review, Page & Uhr 2005, Blath 2025, Baffoe 2026). A link from a dormancy ODE to censored time-to-event data was not found in the abstracts read (bounded, non-systematic search).
3. **Parameters.** All 11 model parameters are hypothetical for human post-treatment disease. The only numbers retrieved (Wilkie & Hahnfeldt, mouse fibrosarcoma) are ~30× faster in growth rate and structurally different.
4. **Structural identifiability (derived from our equations).** (k, s, E) share a scale redundancy; at low burden only ks/μ matters; ρ, h, η are invisible at first order in the tumour-free regime.
5. **Observation from existing code.** A stable nonzero equilibrium (dormant fraction 0.87) exists for k < k_c = 0.2833 and shrinks to the tumour-free state as k → k_c⁻ (transcritical signature; not yet formally verified). Hypothesis H3 in its crude form (relapse time ≈ ln(N_det/N0)/λ₊) failed a first check (≈240 vs ≈375 months at 0.5 k_c).
6. **METABRIC export.** Has `RFS_STATUS`/`RFS_MONTHS` (not DFS/PFS). RFS is a composite (loco-regional, distant, disease-specific death) with non-disease deaths censored; origin presumed diagnosis (RFS ≤ OS in all 1,981 cases) but not stated; no treatment dates; no loco-regional/distant split.
7. **Counts.** 2,380 usable records, 964 events, median potential follow-up 157.7 months. KM RFS: 73.7% (5 y), 62.4% (10 y), 54.4% (15 y), 49.5% (20 y; 162 at risk). Landmark after 5 y, ER+: 1,225 patients, 312 later events.

## Dataset decision
**CONDITIONAL GO**: group-level, exploratory comparison at 5 and 10 years (15 and 20 years cautiously). **NO-GO** for patient-level prediction and for any claim of recovering dormant-cell parameters. TCGA-BRCA is unsuitable for long-term relapse (10 DFS events after 60 months). Datasets are not combined. Conditions are listed in `data/dataset_feasibility.md` §8.

## Biological assumptions and support
| Assumption | Support | Contest |
|---|---|---|
| Dormancy exists as cellular state and delays relapse | Aguirre-Ghiso 2007 (abs.) | |
| P↔D switching is a reasonable model | Blath 2025, Page & Uhr 2005 (abs.) | Continuous phenotype (Baffoe 2026, abs.) |
| Immune control matters for dormancy | Wilkie & Hahnfeldt 2013 (FT) | |
| Immune acts only on P | none | Mehdizadeh 2021, Saudemont 2005 (abs.) |
| Constant carrying capacity | none | dynamic K_C (Wilkie & Hahnfeldt, FT) |
| Late relapse hazard is non-zero | Pan 2017, Rueda 2019 (abs.); METABRIC data | |

## Model limitations (recorded, not fixed)
Binary dormancy; mass-action killing without saturation; no environmental feedback or immunoediting; constant parameters; no treatment; deterministic; E-gauge redundancy; hypothetical parameters; N_detect = 0.1 is a design choice. Proposals M1–M8 are in `research_gap_and_objectives.md` §8.

## Recommended next phase
**Phase 2 (formal analysis):** (a) continuation-based verification of the transcritical branch and search for a second boundary or interior saddle; (b) fix the E-gauge in a *new* config rather than altering the baseline; (c) re-examine H3 (eigenvector projection and nonlinear corrections); (d) then local/global sensitivity and Monte Carlo (Phases 3–4). Clinical comparison (Phase 5) should include the offset τ.

## Remaining uncertainties and blockers
- Full texts of Rueda 2019, Wilkie & Hahnfeldt 2013 (*Interface Focus*), Yáñez Feliú 2026, Blath 2025, Kuznetsov 1994, de Pillis 2005 are unread; institutional access or manual download is needed before the final report quotes beyond abstracts.
- Time origin of METABRIC RFS is inferred, not documented; the original METABRIC/EGA documentation was not consulted.
- `survival` tooling (lifelines or scikit-survival) is not installed; Phase 1 used a hand-written Kaplan–Meier, checked only for internal consistency, not against a library.
- TCGA-CDR and EGA METABRIC data not examined.

## Files created in Phase 1
`literature_review.md`, `research_gap_and_objectives.md`, `parameter_evidence.md`, `PHASE_1_SUMMARY.md`, `report/references.bib`, `data/dataset_feasibility.md`, `data/clinical_data_dictionary.csv`, `data/README.md`, `experiments/metabric_audit.py`, `experiments/tcga_audit.py`, `results/tables/metabric_*.csv|txt`, `results/tables/tcga_audit_log.txt`, raw downloads in `data/raw/`.
