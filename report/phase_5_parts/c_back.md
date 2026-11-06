## 13. Discussion

### 13.1 What the clinical analysis shows

The clinical analysis characterises the *exported RFS variable*, an event that combines loco-regional relapse, distant relapse and disease-specific death, in 2,380 usable patients with 964 events. The headline estimates, the ER contrast and the late-event pattern are given in Section 7; the key point here is their meaning. The cumulative incidence of the RFS event is 0.263 at 5 years, 0.376 at 10 years, 0.456 at 15 years and 0.505 at 20 years (Table 7.1), with 162 patients still at risk at 20 years. ER-negative patients accumulate events much faster early (0.390 against 0.225 at 5 years), but the curves converge by about 15 years and the 20-year values overlap (0.497 against 0.515, Table 7.2); after year 5 the interval rate is higher in ER-positive patients in every window. The interval event rate falls from 6.14 per 100 person-years in years 0–5 to 3.35, 2.82 and 2.09 in the following 5-year windows, and is still 2.61 (95% CI 1.74–3.77) per 100 person-years for ER-positive patients in years 15–20 (Table 7.4). A non-zero event rate beyond ten years in the ER-positive group is descriptively compatible with the late-recurrence literature (Pan et al. 2017; Rueda et al. 2019), but the composite event and undocumented origin prevent reading it as reactivation of dormant cells.

### 13.2 Biological meaning of the compartments

The model's three state variables are mathematical aggregates, not measured quantities.

- **$P$ (proliferating burden).** Cells that divide, can be killed by immune effectors, and can enter dormancy at rate $q$. It is the population whose growth rate $r$ and immune killing $\beta=ks/\mu$ control the stability of the tumour-free state.
- **$D$ (dormant burden).** A non-dividing reservoir that returns to $P$ at rate $a$ and is lost at rate $\delta_D$. Cellular dormancy is an established concept (Aguirre-Ghiso 2007), and state-switching has been modelled elsewhere (Page and Uhr 2005; Blath et al. 2025). Continuous dormancy depth is an alternative view (Baffoe and Bouchnita 2026), so binary dormancy is a simplification.
- **$E$ (immune activity).** An aggregate effector level that is replenished at rate $s$, stimulated by $P$ and suppressed by total burden. Whether the immune system acts on dormant cells is contested (Saudemont et al. 2005; Mehdizadeh et al. 2021). The model assumes that it acts only on $P$.
- A structural property of the model is that the dormant fraction at every non-zero equilibrium equals $q/(q+a+\delta_D)=0.8696$ at baseline. It is an identity of the switching rates and not an outcome of immune control.

No parameter has quantitative support for human post-treatment disease. Only mouse values with different functional forms were found (Wilkie and Hahnfeldt 2013).

### 13.3 How simulated relapse differs from clinical recurrence

Simulated relapse is a deterministic first crossing of $N_{\mathrm{detect}}$ for one sampled parameter vector, whereas clinical RFS is a composite event recorded under treatment, imaging and follow-up schedules. A simulated curve describes parameter and initial-condition uncertainty under a prior, not individual patient risk. Section 4.4 lists the five main non-equivalences. The curves can therefore only be compared in shape (timing, early versus late events, plateau), and then only conditionally on $\tau$, $d$ and the threshold.

### 13.4 Interpreting the model-versus-clinical comparison

In the reference scenario the simulated cumulative incidence lies outside the METABRIC 95% CI at all four horizons. Across the 1,764 prespecified scenarios the share inside the CI is 0.0% at 5 years, 1.1% at 10 years, 7.3% at 15 years and 11.1% at 20 years (Table 10.1). The detection threshold is the strongest single factor: at 20 years the share inside the CI is 0% for thresholds of 0.1 and above, but 11.1%, 11.5% and 55.2% for thresholds 0.01, 0.02 and 0.05. The priors give reference values of 0.317 to 0.468 at 20 years, and the offset $\tau$ changes the shares by a few percentage points (Sections 10.3, 11, 12). The simulation has a robust deficit of early events (statement S1 holds in all 1,764 scenarios) and a robustly later, more back-loaded pattern (S5: 0.791 of the 240-month events occur after 120 months in the reference scenario, against 0.256 in METABRIC).

The following interpretation framework is fixed in advance and does not depend on the numbers.

- If the observed curve lies *inside* the envelope at a horizon, the model family is not contradicted there. It does not follow that the mechanism is correct, because Phase 4 showed that very different parameter vectors give nearly the same event-time curve and that an offset of up to about 12 months cannot be separated from a factor-2 change in initial burden at METABRIC sample size.
- If the observed curve lies *outside* the envelope at a horizon, that statement applies only to the declared hypothetical prior, threshold and nuisance grid, and to the structure of the model (no treatment, no early-recurrence mechanism). It cannot be attributed to any single biological cause.
- A finding counts as robust only if it holds across the whole pre-specified grid. Any change made after viewing the comparison is labelled exploratory.

### 13.5 Detection threshold and time origin

The threshold alone moves the simulated 20-year event probability by a factor of about three (Phase 3), and the nuisance offset $\tau+d$ shifts the curve exactly (Phase 4). The METABRIC RFS origin is undocumented, so Phase 5 treats $\tau$ as a range, not a value. An offset $\tau+d$ of 36 months lowers the reference $F_{\mathrm{sim}}$ at 10 years from 0.091 to 0.010 and at 20 years from 0.436 to 0.370 (Table 12.1), and it moves the simulation further from METABRIC at the early horizons. Under prior P1_base the simulated 20-year probability ranges from 0.699 at threshold 0.01 to 0.218 at threshold 0.5, against 0.505 observed (Table 11.1).

### 13.6 Identifiability and why clinical calibration was excluded

Phase 4 pre-registered three decisions in `configs/phase4.yaml`, before the results were known.

- **A. Descriptive population-level comparison: yes.** The endpoint is documented, the cohort is large (2,380 usable, 964 events) and subgroup reliability can be stated. This is the scope of Phase 5.
- **B. Restricted exploratory population calibration: not supported.** A finite profile interval for a population combination, with $\tau$ and the threshold or $N_0$ scale free, was obtained only in an idealised case with one known heterogeneity source. With two varying parameters the location is a ridge of about one log unit, and $\tau=12$ months is indistinguishable from a rescaling of initial burden by 0.5 at METABRIC size.
- **C. Mechanistic calibration: unsuitable.** An event time is one scalar, so every parameter lies on a nine-dimensional constant-time manifold. From total burden alone only the initial burden $N_0$ and the dominant eigenvalue $\lambda_+$ of the low-burden block are practically identifiable. The parameters $\delta_D$ and $\eta$ are not identifiable even from $(P,D,e)$. METABRIC supplies event times only.

Phase 5 respects these decisions: no biological parameter is fitted, and no simulation is chosen or tuned to resemble the clinical curve. Good agreement of event-time curves would not distinguish between mechanisms, and Phase 4 showed that good fits can coexist with parameters wrong by factors of 5–10.

### 13.7 Unknown clinical time origin

The RFS origin is the largest single barrier to aligning model and data. The documentation supports only a shared clock for RFS and OS. Candidate origins (diagnosis, surgery, end of therapy) differ by an unrecorded interval that may be months. Because the effect is an exact shift in the simple case, it is explored by sensitivity analysis and not estimated.

### 13.8 Future work

1. **Stochastic birth–death–switching model (M6)** for small residual populations (Blath et al. 2025), because all present variability comes from parameters and initial conditions.
2. **Data that could support mechanistic work.** Paired proliferating-like and dormant-like (or immune) biomarker series on an absolute scale, visits through the nonlinear regime, and an anchor for the time origin. Phase 4 showed these would make $r$, $q$, $a$, $\beta$, $\rho$ and the dormant fraction identifiable.
3. **Model-structure alternatives** recorded in Phase 1: saturating killing (M2), immune action on dormant cells (M3), an environmental signal driving reactivation (M4; Yáñez Feliú et al. 2026), continuous dormancy depth (M7) and an explicit treatment-offset parameter (M8).
4. **Competing-risk and endpoint-resolved analyses** if a dataset with event causes, dates and treatment timing becomes available.
5. **A restricted scope-B analysis** only on separate approval, using a phenomenological event-time family (location, scale, cure fraction) and not mechanistic parameters.

## 14. Limitations

1. **Undocumented RFS time origin.** Only a shared clock with OS is supported. The offset to the model clock is unknown and is treated as a nuisance.
2. **Composite endpoint.** RFS mixes loco-regional relapse, distant relapse and disease-specific death. These cannot be separated, and death may be the first recorded event in up to 10.5–27.2% of events, so the observed curves are not a measure of tumour reactivation.
3. **Other-cause deaths are censored.** One minus the Kaplan–Meier estimate therefore overstates recurrence incidence relative to a competing-risk estimate, and no competing-risk analysis is possible without event-cause times.
4. **No treatment dates or durations.** Treatment is available only as flags, so the diagnosis-to-treatment interval and treatment effects on early events are not modelled or adjusted.
5. **Deterministic model without treatment or an early-recurrence mechanism.** The model starts from post-treatment residual disease. It has no mechanism for early recurrences or for treatment effects, and its variability is parameter and initial-condition uncertainty only.
6. **Hypothetical parameters and priors.** No value has an empirical basis for human post-treatment disease. The simulated probabilities depend on the prior (a spread of 0.15 in the 20-year probability across seven priors) and are model-defined quantities, not patient risks.
7. **Arbitrary detection threshold.** The threshold is a design choice and not an imaging limit. Its effect on simulated probabilities is as large as that of the prior.
8. **Limited identifiability.** Mechanistic parameters are not practically identifiable from burden or event data, and the time offset, threshold and initial burden are confounded in event times. Agreement or disagreement cannot be attributed to a specific mechanism.
9. **Scope B not pursued; no calibration.** By design, no parameter was fitted, optimised or selected to match the clinical curve.
10. **Descriptive consistency only.** A comparison of curves is not validation. Any change to the scenario grid or criteria after the comparison is seen is exploratory.
11. **Overlapping subgroups and multiplicity.** ER and PAM50 strata overlap, tests are not adjusted across families or horizons, and subgroup differences are associations.
12. **Small late at-risk numbers.** At 20 years only 162 patients remain at risk overall, and only the whole cohort, ER-positive, ER-negative and Luminal A pass the reliability criterion. Estimates beyond 15 years are therefore wide and are interpreted cautiously.
13. **Missing fields.** About 21% of usable patients lack the extended subgroup fields, and 161 events lack vital status.
14. **Single cohort.** METABRIC is one dataset with its own selection and treatment era. Nothing here generalises to other populations or settings.

## 15. Conclusion

Phase 5 connects the nonlinear dormancy model to observed clinical recurrence patterns only at a descriptive level. The clinical analysis reports Kaplan–Meier estimates of the exported RFS event for 2,380 patients, with the origin, composite structure and censoring convention stated as limitations (cumulative incidence 0.263, 0.376, 0.456 and 0.505 at 5, 10, 15 and 20 years, from 964 events). The model-based analysis converts the existing Phase 3 and 4 simulations to event-free survival curves and compares them with the observed curves over a pre-specified grid of priors, thresholds, time offsets and delays (the model robustly produces too few early events and a more back-loaded late pattern than METABRIC (statements S1 and S5 hold in all 1,764 scenarios), whereas agreement at a single 20-year horizon occurs in only 11.1% of scenarios and depends mainly on the assumed threshold).

Three conclusions are independent of the numbers.

1. The comparison cannot validate the mechanism. Event times carry one scalar per patient, the nuisances are confounded, and Phase 4 showed that the mechanistic parameters are not identifiable.
2. The unknown RFS time origin and the arbitrary detection threshold are first-order sources of uncertainty, comparable to or larger than the prior assumptions, and they must accompany every statement about agreement.
3. The defensible use of the model at present is as a hypothesis-generating description of how dormancy, immune control and a detection threshold can produce delayed events. It is not a predictor of any patient's outcome.

Making progress beyond this point needs different data (biomarker time series, a documented time origin, event-cause resolution) or a stochastic extension of the model, as listed in Section 13.8.

## 16. References

References are limited to entries in `report/references.bib`. Entries marked (abstract) were verified at abstract or metadata level only (Phase 1).

- Aguirre-Ghiso, J. A. (2007). Models, mechanisms and clinical evidence for cancer dormancy. *Nature Reviews Cancer*. (abstract)
- Baffoe, E. N. and Bouchnita, A. (2026). A multiscale model reveals how ERK/p38-regulated dormancy shapes tumor-immune dynamics and immunoediting outcomes. *Frontiers in Immunology*. (abstract)
- Blath, J., Kraut, A., Paul, T. and Töbiás, A. (2025). A stochastic population model for the impact of cancer cell dormancy on therapy success. *Journal of Theoretical Biology*. (abstract)
- cBioPortal for Cancer Genomics datahub (2026). Breast Cancer (METABRIC, Nature 2012 & Nat Commun 2016), study `brca_metabric`.
- cBioPortal for Cancer Genomics datahub (2026). Breast Invasive Carcinoma (TCGA, PanCancer Atlas), study `brca_tcga_pan_can_atlas_2018`.
- Curtis, C., Shah, S. P., Chin, S.-F., Turashvili, G., Rueda, O. M., Dunning, M. J., et al. (2012). The genomic and transcriptomic architecture of 2,000 breast tumours reveals novel subgroups. *Nature*.
- Mehdizadeh, R., Shariatpanahi, S. P., Goliaei, B., Peyvandi, S. and Rüegg, C. (2021). Dormant tumor cell vaccination: a mathematical model of immunological dormancy in triple-negative breast cancer. *Cancers*. (abstract)
- Page, K. and Uhr, J. W. (2005). Mathematical models of cancer dormancy. *Leukemia & Lymphoma*. (abstract)
- Pan, H., Gray, R., Braybrooke, J., Davies, C., Taylor, C., McGale, P., et al., EBCTCG (2017). 20-year risks of breast-cancer recurrence after stopping endocrine therapy at 5 years. *New England Journal of Medicine*. (abstract)
- Pereira, B., Chin, S.-F., Rueda, O. M., Vollan, H.-K. M., Provenzano, E., Bardwell, H. A., et al. (2016). The somatic mutation profiles of 2,433 breast cancers refines their genomic and transcriptomic landscapes. *Nature Communications*.
- Rueda, O. M., Sammut, S.-J., Seoane, J. A., Chin, S.-F., Caswell-Jin, J. L., Callari, M., et al. (2019). Dynamics of breast-cancer relapse reveal late-recurring ER-positive genomic subgroups. *Nature*.
- Saudemont, A., Jouy, N., Hetuin, D. and Quesnel, B. (2005). NK cells that are activated by CXCL10 can kill dormant tumor cells that resist CTL-mediated lysis and can express B7-H1 that stimulates T cells. *Blood*. (abstract)
- Tolaney, S. M., Garrett-Mayer, E., White, J., Blinder, V. S., Foster, J. C., Amiri-Kordestani, L., et al. (2021). Updated standardized definitions for efficacy end points (STEEP) in adjuvant breast cancer clinical trials: STEEP version 2.0. *Journal of Clinical Oncology*. (abstract)
- Wilkie, K. P. and Hahnfeldt, P. (2013). Tumor-immune dynamics regulated in the microenvironment inform the transient nature of immune-induced tumor dormancy. *Cancer Research*. (full text)
- Yáñez Feliú, G., Rossato, G., Valleriani, A. and Cipitria, A. (2026). Dynamic cancer dormancy and awakening emerge from tumor microenvironment feedback in a minimal theoretical model. *npj Systems Biology and Applications*. (abstract)
