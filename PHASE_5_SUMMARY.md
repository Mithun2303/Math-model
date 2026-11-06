# Phase 5 Summary: Clinical Survival Analysis and Descriptive Model Consistency (METABRIC)

*Hypothetical model, illustrative parameters. The clinical part is a descriptive analysis of the exported METABRIC RFS variable (composite event, undocumented time origin). The model-versus-clinical part is a pre-registered consistency check. **No biological parameter was fitted to METABRIC and no simulation was selected or tuned to resemble it.** Nothing here is validation or a patient risk. Full text: `report/phase_5_clinical_comparison.md`. Design and hash: `configs/phase5.yaml` (sha256 `9081048b88634da27432cb42035d87042f68b674ed1a570ae62e6199583ab39f`, recorded 2026-11-05 18:55:30 UTC before METABRIC was loaded for any comparison).*

## Tests
**164 passed, 0 failed** (`python -m pytest`): the 138 Phase 1–4 tests unchanged + 16 new clinical-survival tests (`tests/test_phase5_clinical.py`) + 10 new model-consistency tests (`tests/test_phase5_model_consistency.py`). One earlier run, made while other processes were using the CPU, stopped at a single failure in `test_phase3_uncertainty.py::test_monte_carlo_estimate_converges_and_parallel_equals_serial`; the clean full run passed. The cause was not investigated, so treat that test as possibly timing/CPU-sensitive on this Windows machine.

## Files
- **Created:** `configs/phase5.yaml`; `src/dormancy/clinical_comparison.py`; `experiments/model_consistency.py`, `experiments/phase5_figures.py`; `tests/test_phase5_clinical.py`, `tests/test_phase5_model_consistency.py`; `results/p5_rfs_time_origin_note.md`; 18 tables `results/tables/p5_*.csv` plus `p5_prereg_hash.txt`; 10 figures `results/figures/fig_p5_01 … 10_*.png`; `report/phase_5_clinical_comparison.md` (Sections 1–16); this file.
- **Modified:** `experiments/phase5_metabric_km.py` (existing script, now run); the Phase 4 `src/dormancy/km_descriptive.py` is reused.
- Phase 1–4 sources, tables, figures and configs were not changed.

## Clinical results (METABRIC, `RFS_STATUS`/`RFS_MONTHS` as exported)
- **Records:** 2,509 patients in the export; **2,380 usable, 964 recorded RFS events**, 1,416 censored (15 of 15 coding checks pass; 0 patients with RFS > OS; 4 with time 0 kept).
- **Kaplan–Meier cumulative incidence (95% log–log CI; at risk):** 5 y **0.263** (0.245–0.281; 1,581); 10 y **0.376** (0.355–0.397; 984); 15 y **0.456** (0.433–0.480; 521; no reliability criterion); 20 y **0.505** (0.478–0.532; 162).
- **ER status:** ER-positive 0.225 / 0.358 / 0.515 and ER-negative 0.390 / 0.442 / 0.497 at 5 / 10 / 20 y; log-rank $p<0.001$ (5, 10 y) and $p=0.0002$ (20 y). ER-negative events are early (10.3 per 100 person-years in years 0–5, 1.8 in years 5–10); ER-positive rates are higher after year 5 and still 2.6 per 100 person-years in years 15–20.
- **PAM50:** only strata meeting the criterion (≥ 50 events, ≥ 100 at risk at 5 and 10 y, ≥ 30 at 20 y) are interpreted. Luminal A is reliable to 20 y (0.156 / 0.267 / 0.461), Luminal B and claudin-low to 10 y, Her2, Basal and Normal only to 5 y (Basal 0.409, Her2 0.391 at 5 y). Holm-adjusted pairwise differences at 5 y: Luminal A differs from Luminal B, Her2, Basal and Normal.
- **Late recurrence (descriptive):** interval rate 6.14, 3.35, 2.82, 2.09 per 100 person-years in years 0–5, 5–10, 10–15, 15–20. Landmark at 5 y (1,581 event-free): conditional incidence 0.153, 0.262, 0.329 over the next 5, 10, 15 years.
- **Death versus recurrence:** death can be the first recorded event in at most 10.5% (known vital status) to 27.2% (if all missing-status events) of events; 72.8% have documented relapse before death or no death. No competing-risk analysis (event-cause times unavailable).
- **Time origin:** undocumented; only a shared clock with OS is supported (745 of 745 and 432 of 432 exact equalities). Every time is "since the unknown RFS origin".

## Model-versus-clinical consistency (pre-registered grid)
- **Grid:** 7 priors × 7 thresholds × 4 offsets $\tau$ (0, 6, 12, 24 mo) × 3 delays $d$ (0, 6, 12 mo) × 3 death-first shares (0, 0.105, 0.272) = **1,764 scenarios**; reference scenario fixed in advance (P1_base, threshold 0.1, $\tau=d=0$). Simulated curves reuse the Phase 3 first-crossing samples (no new sweep). The `exploratory` flag is false for all grid rows.
- **Reference scenario:** $F_{\mathrm{sim}}$ = 0.0005 / 0.0910 / 0.3083 / 0.4363 at 5 / 10 / 15 / 20 y against 0.263 / 0.376 / 0.456 / 0.505; outside the METABRIC CI at all four horizons; supremum distance over 0–240 months 0.310.
- **Share of scenarios inside the METABRIC CI:** 0.0% (5 y), 1.1% (10 y), 7.3% (15 y), 11.1% (20 y).
- **Robust across all 1,764 scenarios:** S1 (simulated incidence below the CI at 5 y) and S5 (a larger share of simulated events after 10 y than observed: 0.791 against 0.256 in the reference scenario).
- **Assumption-dependent:** S2 (below CI at 10 y) 98.8%; S3 (inside CI at 20 y) 11.1%; S4 (later late increment) 99.6%; S6 (ratio $F(240)/F(120)>2$) 94.3%.
- **Detection threshold:** for P1 the 20-year probability goes from 0.699 (θ = 0.01) to 0.218 (θ = 0.5), the 10-year from 0.360 to 0.008, the median crossing time from 118 to 183 months. Share of scenarios inside the CI at 20 y: 11.1%, 11.5%, 55.2% for θ = 0.01, 0.02, 0.05 and 0% for θ ≥ 0.1.
- **Time origin and delay:** $\tau+d$ translates the curve exactly; 36 months lowers the P1 reference $F(10\,\mathrm y)$ from 0.091 to 0.010 and $F(20\,\mathrm y)$ from 0.436 to 0.370, i.e. it widens the early gap. The death-first mapping changes $F(20\,\mathrm y)$ by at most 0.012.
- **Alignment diagnostic (exploratory, descriptive only):** the minimum is at the boundary (−60 months) for all seven priors, so there is no interior minimum and no offset estimate.

## Counts
10 figures (`fig_p5_01 … 10`), 18 CSV tables + 1 hash text file, 8 numbered tables in the report. Report length about 10,000 words (Sections 1–16).

## Limitations
Undocumented RFS origin; composite event with censored other-cause deaths; no treatment dates; deterministic model with no treatment or early-recurrence mechanism; hypothetical priors and an arbitrary detection threshold; limited identifiability and exact confounding of offset, threshold and initial burden; subgroup strata overlap and are unadjusted across families; few patients at risk beyond 15 years; single cohort. The comparison supports only descriptive statements about timing and shape.

## Remaining research questions
1. Does a stochastic birth–death–switching extension (M6) change the early-event deficit, or is a treatment-era or early-recurrence mechanism needed?
2. Which data (documented time origin, event-cause times, treatment dates, paired proliferating/dormant or immune biomarkers) would make mechanistic parameters identifiable?
3. Can the composite endpoint be split into loco-regional, distant and death events in another cohort?
4. Would a restricted phenomenological scope-B analysis (location, scale, cure fraction) be informative, and does it need separate approval?
5. How should the threshold be linked to a clinical detection process rather than a burden fraction?
