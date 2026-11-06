# Tumour dormancy and relapse: P-D-E ODE model

Hypothetical research simulations only; parameters are illustrative, not fitted, and nothing here is a clinical prediction.

    python3.12 -m venv .venv && .venv/bin/pip install -r requirements.txt && .venv/bin/pip install -e .
    .venv/bin/python -m pytest
    .venv/bin/python experiments/baseline.py

- `src/dormancy/`: model, solver (first-upward-crossing relapse event, censored at horizon), equilibria, stability
- `configs/baseline.yaml`: parameters with units
- Tumour-free stability: `r < k s/mu + q delta_D/(a+delta_D)`; `critical_k(p)` gives the boundary in k.

Phase 2 (bifurcation analysis, M1 normalisation, relapse-time approximations; see `PHASE_2_SUMMARY.md`, `report/phase_2_mathematical_analysis.md`):

    .venv/bin/python experiments/bifurcation_analysis.py
    .venv/bin/python experiments/relapse_time_approximation.py
    .venv/bin/python experiments/normalized_model_comparison.py

- `src/dormancy/continuation.py` (pseudo-arclength continuation), `bifurcation.py` (Sotomayor conditions, folds, branch scans), `relapse.py` (relapse-time estimators), `normalized.py` (e = mu E/s, beta = k s/mu; `configs/normalized.yaml`)


Phase 3 (global sensitivity, basins, Monte Carlo relapse simulation; see `PHASE_3_SUMMARY.md`, `report/phase_3_sensitivity_and_uncertainty.md`; design in `configs/phase3.yaml`):

    .venv/bin/python experiments/local_sensitivity.py
    .venv/bin/python experiments/global_sensitivity.py      # ~20 min the first time (Sobol design is cached)
    .venv/bin/python experiments/basin_analysis.py          # --replot regenerates figures from saved tables
    .venv/bin/python experiments/monte_carlo_relapse.py
    .venv/bin/python experiments/threshold_robustness.py    # needs the Monte Carlo run first

- `src/dormancy/uncertainty.py` (priors, sampling, dense trajectory runner, long-run classifier, aggregation), `sensitivity.py`, `basins.py`. All simulated probabilities are model-defined and hypothetical.


Phase 4 (identifiability, synthetic validation, endpoint mapping; see `PHASE_4_SUMMARY.md`, `report/phase_4_identifiability_and_validation.md`; design and pre-registered decision criteria in `configs/phase4.yaml`):

    .venv/bin/python experiments/identifiability_analysis.py      # structural (exact) + Fisher + profile likelihoods (~10-30 min)
    .venv/bin/python experiments/synthetic_recovery.py            # constrained multistart recovery (--figures-only to replot)
    .venv/bin/python experiments/model_reduction.py && .venv/bin/python experiments/model_reduction_combinations.py
    .venv/bin/python experiments/time_origin_sensitivity.py
    .venv/bin/python experiments/survival_validation.py
    .venv/bin/python experiments/population_identifiability_nuisance.py && .venv/bin/python experiments/population_identifiability_nuisance2.py
    .venv/bin/python experiments/metabric_endpoint_investigation.py   # descriptive audit only; no model is fitted

- `src/dormancy/identifiability.py` (exact series/AD identifiability), `observation.py`, `synthetic_validation.py`, `endpoint_mapping.py`. All Phase 4 results are SYNTHETIC; no clinical calibration was performed.

Not yet built: clinical (METABRIC) comparison, parameter fitting, stochastic (birth-death-switching) model, Streamlit app.


Phase 5 (descriptive METABRIC survival analysis and pre-registered model-vs-METABRIC consistency; see `PHASE_5_SUMMARY.md`, `report/phase_5_clinical_comparison.md`; design and hash in `configs/phase5.yaml`):

    python experiments/phase5_metabric_km.py      # METABRIC KM, subgroups, log-rank, landmark, late hazard (Figs 1-3)
    python experiments/model_consistency.py       # needs Phase 3 Monte Carlo samples; refuses to run if configs/phase5.yaml hash changed (Figs 4-9)
    python experiments/phase5_figures.py          # Fig 10 summary

- `src/dormancy/km_descriptive.py`, `clinical_comparison.py`. No parameter is fitted to METABRIC; all results are descriptive.
