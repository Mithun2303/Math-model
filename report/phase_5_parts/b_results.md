## 7. Clinical Survival Results

All results in this section describe the exported METABRIC RFS variable: a composite of loco-regional relapse, distant relapse and disease-specific death, with an undocumented time origin. No model is involved in this section. "Cumulative incidence" means $1-\hat S(t)$ from the Kaplan–Meier estimator, with other-cause deaths censored (Section 6).

### 7.1 Cohort and event coding checks

The coding was re-verified in `p5_clinical_coding_checks.csv`. All 15 checks pass.

- The export has 2,509 patients. `RFS_STATUS` has 1,486 "Not Recurred", 1,002 "Recurred" and 21 missing. `RFS_MONTHS` is missing for 121 patients.
- The usable cohort, with both time and status, is **2,380 patients with 964 events**. Four usable patients have a time of exactly zero and are kept, as in Phase 1.
- No usable patient has a negative RFS time, and no patient has RFS longer than OS (0 violations among patients with both times).
- The maximum RFS time among usable patients is 389.33 months, and 9 events occur after 240 months.
- The cumulative incidence at 5, 10, 15 and 20 years reproduces the Phase 1 table (Table 7.1).

### 7.2 Overall recurrence-free survival

**Table 7.1.** Kaplan–Meier cumulative incidence of the RFS event, all usable patients (source: `results/tables/p5_km_overall.csv`). CI is the 95% log–log Greenwood interval. "At risk" counts patients with observed time at or beyond the horizon. Events are counted up to and including the horizon.

| Horizon | At risk | Events by horizon | Cumulative incidence | 95% CI | Event-free $\hat S$ | 95% CI of $\hat S$ |
|---|---|---|---|---|---|---|
| 5 y | 1,581 | 604 | 0.2628 | 0.2453–0.2814 | 0.7372 | 0.7186–0.7547 |
| 10 y | 984 | 819 | 0.3757 | 0.3554–0.3968 | 0.6243 | 0.6032–0.6446 |
| 15 y | 521 | 923 | 0.4559 | 0.4330–0.4795 | 0.5441 | 0.5205–0.5670 |
| 20 y | 162 | 955 | 0.5050 | 0.4783–0.5324 | 0.4950 | 0.4676–0.5217 |

The 15-year estimate has no Phase 4 reliability criterion and is shown for the whole cohort only. The 20-year estimate rests on 162 patients still at risk.

![Figure 7.1. METABRIC recurrence-free survival event, all usable patients](../results/figures/fig_p5_01_metabric_km_overall.png)

**Figure 7.1.** METABRIC recurrence-free survival event, all usable patients. *Setting:* 2,380 usable patients, 964 events, Kaplan–Meier estimate of the cumulative incidence of the exported RFS event with a 95% log–log Greenwood band. *Axes/panels:* x, years since the (undocumented) RFS time origin; y, cumulative incidence $1-\hat S$; markers at 5, 10, 15 and 20 years carry the values 0.263, 0.376, 0.456 and 0.505; the strip below lists numbers at risk (2,380, 1,581, 984, 521, 162, 18 at 0, 5, 10, 15, 20, 25 years). *What it shows:* a steep rise in the first 5 years that slows steadily; the curve keeps rising after 10 years. *Mathematical interpretation:* the interval hazard falls from 6.14 to 2.09 events per 100 person-years between the first and the last 5-year window (Table 7.4), so the curve is concave but not flat. The band widens after 20 years as the risk set shrinks. *Biological/clinical interpretation:* a non-zero event rate after 10 years is descriptively compatible with the late-recurrence literature, but the endpoint includes disease-specific death and the origin is unknown, so no process is inferred. *Limitations:* composite event, other-cause deaths censored, undocumented origin, wide band beyond 20 years.

### 7.3 ER-positive versus ER-negative

**Table 7.2.** Cumulative incidence of the RFS event by ER (IHC) status at the horizons where each group meets the reliability criterion (source: `results/tables/p5_km_subgroups.csv`, `p5_logrank_subgroups.csv`).

| Group | n (events) | Horizon | At risk | Cumulative incidence | 95% CI |
|---|---|---|---|---|---|
| ER-positive | 1,735 (689) | 5 y | 1,225 | 0.2251 | 0.2058–0.2459 |
| | | 10 y | 738 | 0.3584 | 0.3347–0.3831 |
| | | 20 y | 110 | 0.5154 | 0.4824–0.5493 |
| ER-negative | 573 (258) | 5 y | 314 | 0.3897 | 0.3505–0.4317 |
| | | 10 y | 222 | 0.4423 | 0.4014–0.4854 |
| | | 20 y | 49 | 0.4973 | 0.4519–0.5447 |

Log-rank tests between the two groups, administratively censored at each horizon, give $\chi^2=74.5$ at 5 years ($p<0.001$), $\chi^2=32.0$ at 10 years ($p<0.001$) and $\chi^2=14.0$ at 20 years ($p=0.0002$). These are unadjusted associations and the strata overlap with the PAM50 groups.

![Figure 7.2. RFS event by ER status](../results/figures/fig_p5_02_er_positive_vs_negative.png)

**Figure 7.2.** METABRIC RFS event by ER (IHC) status. *Setting:* ER-positive (n = 1,735; 689 events) against ER-negative (n = 573; 258 events), with 95% log–log bands; no adjustment for treatment, stage or grade. *Axes/panels:* x, years since the undocumented origin; y, cumulative incidence; the table below gives numbers at risk (ER-positive 1,735, 1,225, 738, 369, 110; ER-negative 573, 314, 222, 143, 49 at 0, 5, 10, 15, 20 years). *What it shows:* the ER-negative curve rises much faster in the first 3 to 4 years and then flattens; the ER-positive curve rises more slowly but keeps rising, and the two curves converge at about 15 years with overlapping bands. *Mathematical interpretation:* interval hazards (Table 7.4) cross over: ER-negative 10.35 per 100 person-years in years 0–5 against 1.84 in years 5–10, while ER-positive falls from 5.07 to 3.79; after year 5 the ER-positive rate is higher in every interval. *Biological/clinical interpretation:* the pattern is descriptively the same kind of contrast reported in the late-recurrence literature (Rueda et al. 2019), but this analysis cannot say why. *Limitations:* composite event, no covariate adjustment, the ER-negative curve stops at its last event (about 17 years) and only 49 patients remain at risk at 20 years, and the log-rank p-values are descriptive.

### 7.4 PAM50 and claudin-low subgroups

**Table 7.3.** Reliability of the PAM50 and claudin-low strata under the Phase 4 criterion, and the estimates where it holds (source: `results/tables/p5_km_subgroups.csv`). Cumulative incidence at the horizon, 95% CI in brackets. "–" means the criterion is not met and the group is not interpreted.

| Group | n (events) | 5 y | 10 y | 20 y |
|---|---|---|---|---|
| Luminal A | 699 (231) | 0.156 [0.131–0.186] (535 at risk) | 0.267 [0.233–0.304] (339) | 0.461 [0.408–0.517] (61) |
| Luminal B | 475 (227) | 0.303 [0.263–0.348] (299) | 0.461 [0.414–0.511] (171) | – |
| Her2 | 224 (116) | 0.391 [0.330–0.460] (126) | – | – |
| Basal | 209 (93) | 0.409 [0.345–0.479] (116) | – | – |
| claudin-low | 218 (68) | 0.205 [0.156–0.266] (158) | 0.282 [0.224–0.351] (102) | – |
| Normal | 148 (65) | 0.259 [0.195–0.340] (102) | – | – |

Global log-rank tests among the reliable PAM50 strata are $\chi^2=101.2$ at 5 years (six strata) and $\chi^2=48.6$ at 10 years (Luminal A, Luminal B, claudin-low), both $p<0.001$. At 5 years, with Holm adjustment, the Luminal A stratum differs from Luminal B, Her2, Basal and Normal (adjusted $p\le0.021$), and claudin-low differs from Her2 and Basal ($p\le0.0004$) but not from Luminal A (adjusted $p=0.25$). Luminal A against claudin-low at 10 years gives adjusted $p=0.42$, while Luminal A against Luminal B gives $p<0.001$. At 20 years only Luminal A is reliable, so no between-group test is possible.

![Figure 7.3. RFS event by PAM50 and claudin-low subtype](../results/figures/fig_p5_03_molecular_subgroups.png)

**Figure 7.3.** METABRIC RFS event by PAM50 and claudin-low subtype, reliable strata and horizons only. *Setting:* six strata, each drawn only up to the longest horizon (5, 10 or 20 years) at which it meets the criterion; 6 patients with NC subtype and patients with missing subtype are omitted. *Axes/panels:* x, years since the undocumented origin; y, cumulative incidence; the table below lists numbers at risk at 0, 5, 10, 15, 20 years, with "–" where the group is not shown. *What it shows:* Basal and Her2 rise fastest in the first 5 years (0.41 and 0.39); Luminal B is intermediate and reaches 0.46 by 10 years; Luminal A is lowest early (0.16 at 5 years) but is the only stratum followed to 20 years, where it reaches 0.46. *Mathematical interpretation:* the ordering at 5 years is supported by the Holm-adjusted pairwise tests; later ordering cannot be assessed for most strata because the criterion fails. *Biological/clinical interpretation:* the figure only describes how the exported event accumulates in each stratum; it is not a statement about dormancy in any subtype. *Limitations:* unadjusted, overlapping with the ER family, truncated curves, small strata at long horizons, no multiplicity correction across families.

### 7.5 Late-recurrence pattern

**Table 7.4.** Interval-specific event rates of the RFS event, per 100 person-years, with exact Poisson 95% intervals (source: `results/tables/p5_late_recurrence_hazard.csv`). The event counts in the four windows (601, 215, 104, 32) agree with the Phase 1 window table (`p5_late_vs_phase1_windows.csv`).

| Group | 0–5 y | 5–10 y | 10–15 y | 15–20 y |
|---|---|---|---|---|
| All usable (events) | 6.14 [5.66–6.65] (601) | 3.35 [2.92–3.83] (215) | 2.82 [2.30–3.41] (104) | 2.09 [1.43–2.94] (32) |
| ER-positive (events) | 5.07 [4.57–5.61] (376) | 3.79 [3.26–4.37] (186) | 3.30 [2.65–4.07] (89) | 2.61 [1.74–3.77] (28) |
| ER-negative (events) | 10.35 [9.01–11.82] (215) | 1.84 [1.19–2.72] (25) | 1.31 [0.68–2.28] (12) | 0.91 [0.25–2.34] (4) |

The rate declines with time in every group, but it does not fall to zero within 20 years in the whole cohort or in the ER-positive group. The ER-negative rate is high early and low later, and its late intervals rest on 25, 12 and 4 events.

**Landmark analysis.** Among patients event-free and observed at 5 years (1,581 patients), the cumulative incidence of a later event is 0.153 (95% CI 0.135–0.173) over the next 5 years, 0.262 (0.237–0.289) over 10 years and 0.329 (0.297–0.363) over 15 years. For ER-positive patients (1,225 at the landmark) the 5-year and 15-year conditional values are 0.172 (0.151–0.196) and 0.375 (0.336–0.416); for ER-negative patients (314 at the landmark) they are 0.086 (0.059–0.125) and 0.176 (0.130–0.237). For Luminal A (535 at the landmark) they are 0.131 (0.104–0.165) and 0.361 (0.304–0.425). These are conditional on being event-free at 5 years and are not comparable to the curves in Figures 7.1 to 7.3.

### 7.6 Death versus recurrence in the composite event

The composite endpoint cannot be split, only bounded (source: `results/tables/p5_death_vs_recurrence_bounds.csv`).

- At most **10.5%** of the 964 events can be death-first among events with known vital status. If all 161 events with missing vital status were also death-first, the bound is **27.2%**; the lower bound is 0.
- **72.8%** of events have a documented relapse before death or no death.
- Among 545 patients who died of disease with RFS earlier than OS, the difference OS − RFS has a median of 15.5 months (IQR 6.0–34.1). This is conditional on death and is not a delay distribution.

No competing-risk analysis is performed, because the times of distinct causes are not available.

---

## 8. Model-Generated Event-Time Distributions

### 8.1 Source of the simulated curves

No new simulation was run for Phase 5. The simulated event-free curves are built from the per-sample first-crossing times already produced in Phase 3 (`p3_mc_*_samples.csv`): the validated dense DOP853 runs with relative tolerance $10^{-9}$ and horizon 240 months, with the same samples used at all seven thresholds. The seven hypothetical prior scenarios are those of Phase 3: P1_base (three pooled runs, n = 15,000) and the six alternatives P2 wide, P3 log-uniform, P4 and P5 (rank-correlated), P6 (also uncertain $r$) and P7 (fixed initial state), each with n = 5,000. The model time unit is months, as in the clinical data. The simulated event-free fraction at time $t$ is the fraction of samples whose first crossing time exceeds $t$; because the horizon is the only censoring, this is the empirical survival function. Every simulated probability is a model-defined quantity under a hypothetical prior, not a patient risk.

![Figure 8.1. Model-generated event-free survival](../results/figures/fig_p5_04_model_event_free_survival.png)

**Figure 8.1.** Model-generated event-free survival under seven hypothetical priors. *Setting:* detection threshold 0.1 of $K$, no time offset (reference mapping), sample sizes as listed in the legend. *Axes/panels:* (a) x, months since the residual state (with a years scale on top); y, simulated event-free fraction $1-F(t)$; (b) histogram of first-crossing times among the 6,545 of 15,000 P1 samples that cross within 240 months, with a dashed line at the median, 153 months; the remaining 8,455 are censored at the horizon. *What it shows:* all curves stay at 1.0 for the first 50 to 60 months and then decline; at 240 months they range from about 0.68 (P3) to about 0.53 (P4); P6 declines earliest and P7 latest. In panel (b) crossing times rise from about 50 months, peak between about 110 and 160 months and then decrease towards the horizon, because the horizon truncates the right tail. *Mathematical interpretation:* the curves are convex-then-concave, which reflects a time to cross a threshold from a small initial burden by exponential growth; the spread among priors at 240 months is about 0.15 (Phase 3). The earliest crossing in P1 is 51.9 months. *Biological/clinical interpretation:* the model produces no event before about 4 years under any prior, because it starts from microscopic, mostly dormant residual disease with no treatment-era mechanism. *Limitations:* hypothetical priors, deterministic first-passage times, a threshold that is not an imaging limit, and event-free fraction truncated at 240 months.

---

## 9. Comparison Methodology

### 9.1 Pre-registration

The comparison design is fixed in `configs/phase5.yaml`. Its SHA-256 is `9081048b88634da27432cb42035d87042f68b674ed1a570ae62e6199583ab39f`; it was recorded at 2026-11-05 18:55:30 UTC, **before the METABRIC curve was loaded for any comparison** (`p5_prereg_hash.txt`), and the comparison scripts refuse to run if the file hash differs. The hash prefix `9081048b8863` appears in the footer of every comparison figure.

### 9.2 Scenario grid

The grid is a full factorial of 1,764 scenarios: 7 priors × 7 thresholds (0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5) × 4 offsets $\tau\in\{0,6,12,24\}$ months × 3 delays $d\in\{0,6,12\}$ months × 3 death-first fractions $\{0,0.105,0.272\}$. The death-first fractions are the Phase 4 bounds, with an exponential lag of mean 24 months. The mapping is

$$
t_{\mathrm{clin}}=\tau+t_{\mathrm{model}}+d \quad(\text{recurrence}),\qquad
t_{\mathrm{clin}}=\tau+t_{\mathrm{model}}+\mathrm{Lag},\ \mathrm{Lag}\sim\mathrm{Exp}(24\ \mathrm{mo})\quad(\text{death-first}).
\tag{9.1}
$$

The reference scenario was fixed in advance as the Phase 3 headline definition: prior P1_base, threshold 0.1, $\tau=d=0$, death-first 0. The initial-burden scale is not in the grid, because in the linear regime it is exactly confounded with a time shift (Phase 4) and is therefore represented by the $\tau$ range. Horizons are 60 and 120 months (primary), and 180 and 240 months (secondary; 15 years has no reliability criterion and 20 years has 162 at risk). No comparison horizon exceeds 240 months, where model events are exact for any $\tau,d\ge0$.

### 9.3 Metrics and criterion

For each scenario and horizon $h$ the signed difference $D(h)=F_{\mathrm{sim}}(h)-F_{\mathrm{KM}}(h)$ is computed, together with whether $F_{\mathrm{sim}}(h)$ lies inside the METABRIC 95% CI, the supremum-norm distance over a one-month grid, the late increment $F(240)-F(120)$, the conditional late risk $[F(240)-F(120)]/[1-F(120)]$, the ratio $F(240)/F(120)$ and the fraction of 240-month events after 120 months.

**Criterion.** A scenario is *consistent at horizon $h$* if $F_{\mathrm{sim}}(h)$ lies inside the METABRIC 95% CI of $F_{\mathrm{KM}}(h)$. This means only that the two numbers are not distinguishable at that horizon. A statement is **robust** only if it holds in all 1,764 scenarios (and for every prior); otherwise it is **assumption-dependent** and the fraction of scenarios is reported.

### 9.4 Pre-specified statements

| Statement | Content |
|---|---|
| S1 | $F_{\mathrm{sim}}(60)$ is below the METABRIC 95% CI at 60 months |
| S2 | $F_{\mathrm{sim}}(120)$ is below the METABRIC 95% CI at 120 months |
| S3 | $F_{\mathrm{sim}}(240)$ is inside the METABRIC 95% CI at 240 months |
| S4 | the simulated late increment $F(240)-F(120)$ is larger than the METABRIC late increment |
| S5 | the simulated fraction of 240-month events after 120 months exceeds the METABRIC fraction |
| S6 | $F_{\mathrm{sim}}(240)/F_{\mathrm{sim}}(120)>2$, while METABRIC's ratio is below 2 |

The configuration states the expectation before the comparison: the model has no treatment phase and no early-recurrence mechanism, starts from a mostly dormant microscopic residual state, and its events are threshold crossings, whereas RFS is a composite with an undocumented origin. A deficit of early simulated events is therefore expected, and is neither evidence for nor against the dormancy mechanism.

### 9.5 Exploratory elements and no-tuning rules

No biological parameter is fitted, no prior is adjusted, and no "best" simulation is selected; all 1,764 scenarios are reported (the `exploratory` flag is false for every row of `p5_consistency_grid.csv`). The **alignment diagnostic** (the supremum distance between the METABRIC curve and the reference simulated curve shifted horizontally by $s\in[-60,+60]$ months) was declared a descriptive diagnostic only, excluded from the criterion, and any value taken from it is labelled exploratory. Binomial Monte Carlo (Wilson) intervals for $F_{\mathrm{sim}}$ are shown for transparency and are not used in the criterion.

---

## 10. Model-versus-Clinical Results

### 10.1 Observed against simulated curves

![Figure 10.1. Model event-free curves against METABRIC recurrence-free survival](../results/figures/fig_p5_05_model_vs_metabric.png)

**Figure 10.1.** Model event-free curves against METABRIC RFS (descriptive overlay, no fitting). *Setting:* seven priors at threshold 0.1 with $\tau=d=0$, plus P1_base with $\tau=12$ and 24 months; METABRIC Kaplan–Meier with its 95% band; the grey region from 180 months marks the part with no reliability criterion and few patients at risk. The METABRIC curve ends at about 230 months. *Axes/panels:* upper panel, event-free fraction against months (years scale on top); lower panel, signed difference $F_{\mathrm{sim}}-F_{\mathrm{KM}}$ with the METABRIC 95% band around zero. *What it shows:* METABRIC falls immediately while all simulated curves stay near 1.0 for the first 60 months; the P1 difference reaches about −0.31 between roughly 90 and 110 months and then shrinks to −0.069 at 240 months. P6 (uncertain $r$) is closest to METABRIC and P3 (log-uniform) is farthest at late times. *Mathematical interpretation:* the simulated curves are concave-up at early times while the clinical curve is concave-down, so the discrepancy changes sign of slope, which is a statement about shape; a time offset moves the simulated curves further to the right (the $\tau=12$ and 24 months curves) and does not reduce the early gap. *Biological/clinical interpretation:* the early clinical events (years 0–5 of RFS) have no counterpart in this model, consistent with the stated expectation. *Limitations:* thresholds, priors and offsets are assumptions; an agreement or disagreement in this plot says nothing about the mechanism.

### 10.2 Consistency envelope

**Table 10.1.** Consistency of the 1,764 scenarios with the METABRIC 95% CI (source: `results/tables/p5_consistency_envelope.csv`, group "all").

| Horizon | Scenarios inside CI | Below CI | Above CI | Min / median / max of $D$ | Range of $F_{\mathrm{sim}}$ | $F_{\mathrm{KM}}$ [95% CI] |
|---|---|---|---|---|---|---|
| 5 y (60 mo) | 0 (0.0%) | 100.0% | 0.0% | −0.263 / −0.263 / −0.124 | 0–0.139 | 0.263 [0.245–0.281] |
| 10 y (120 mo) | 20 (1.1%) | 98.8% | 0.1% | −0.376 / −0.310 / +0.041 | 0–0.416 | 0.376 [0.355–0.397] |
| 15 y (180 mo)* | 129 (7.3%) | 79.1% | 13.6% | −0.447 / −0.198 / +0.157 | 0.009–0.613 | 0.456 [0.433–0.480] |
| 20 y (240 mo) | 196 (11.1%) | 65.3% | 23.6% | −0.394 / −0.128 / +0.223 | 0.111–0.728 | 0.505 [0.478–0.532] |

*15 years has no Phase 4 reliability criterion.

The reference scenario (P1_base, threshold 0.1, $\tau=d=0$) is outside the CI at every horizon: $F_{\mathrm{sim}}$ is 0.0005, 0.0910, 0.3083 and 0.4363 against 0.2628, 0.3757, 0.4559 and 0.5050 (differences −0.262, −0.285, −0.148, −0.069; supremum distance over 0–240 months 0.310; `p5_horizon_comparison.csv`).

![Figure 10.2. Consistency envelope](../results/figures/fig_p5_06_consistency_envelope.png)

**Figure 10.2.** Consistency envelope of the model-generated event-free curves against METABRIC. *Setting:* all 1,764 prespecified scenarios (7 priors, 7 thresholds, $\tau$ 0–24, $d$ 0–12, death-first 0–27%). *Axes/panels:* (a) event-free fraction against months with the minimum–maximum band (light blue), the 10th–90th percentile band (dark blue), the METABRIC curve and its CI; (b) for 5, 10, 15 and 20 years, the range (bar) and median (dot) of $F_{\mathrm{sim}}-F_{\mathrm{KM}}$, the METABRIC CI as a grey band around zero and the percentage of scenarios inside the CI (0%, 1%, 7%, 11%). *What it shows:* the METABRIC curve lies below the whole envelope for roughly the first 8 years, inside the min–max band from roughly 8 to 9 years on, and inside the 10th–90th percentile band only from roughly 13 to 15 years on (read from the figure). *Mathematical interpretation:* the width of the envelope (a spread of roughly 0.4 to 0.6 in event-free fraction at late times) is much larger than the METABRIC CI (about ±0.03), so late-horizon agreement is easy to obtain and uninformative about shape. *Biological/clinical interpretation:* only a small fraction of assumption sets give late-horizon agreement, and none gives early agreement. *Limitations:* the grid is hypothetical, scenarios are not weighted by plausibility, and the 15-year and 20-year cells are weakly supported clinically.

### 10.3 Pre-specified statements: robust and assumption-dependent findings

**Table 10.2.** Verdicts on S1–S6 (source: `results/tables/p5_consistency_robust_findings.csv`).

| Statement | Scenarios true | Verdict | Range across the seven priors |
|---|---|---|---|
| S1 early deficit at 5 y | 1,764 / 1,764 (100%) | **robust** | 100% in all priors |
| S2 deficit at 10 y | 1,742 / 1,764 (98.8%) | assumption-dependent | 93.3% (P6) to 100% (P3, P7) |
| S3 agreement at 20 y | 196 / 1,764 (11.1%) | assumption-dependent | 9.5% (P4) to 12.3% (P6) |
| S4 later late increment | 1,757 / 1,764 (99.6%) | assumption-dependent | 97.2% (P5) to 100% |
| S5 back-loaded events | 1,764 / 1,764 (100%) | **robust** | 100% in all priors |
| S6 acceleration ratio above 2 | 1,664 / 1,764 (94.3%) | assumption-dependent | 78.9% (P6) to 100% |

The two robust findings are shape statements. (S1) At 5 years the model is below the METABRIC CI in every scenario. (S5) A larger share of the simulated 240-month events occurs after 120 months than in METABRIC: 0.791 in the reference scenario against 0.256 observed, and a grid median of 0.829. The late increment is 0.345 simulated (reference) against 0.129 observed, the conditional late risk 0.380 against 0.207, and the ratio $F(240)/F(120)$ 4.79 against 1.34. S2, S4 and S6 hold in 94% to 99.6% of scenarios but not in all, so they are assumption-dependent by the pre-specified rule; the per-prior fractions are in Table 10.2. S3 shows that agreement of a single late number is not robust.

**What changes the fraction inside the CI.** Averaging over the grid, the fraction inside the CI at 240 months is 11.1%, 11.5% and 55.2% for thresholds 0.01, 0.02 and 0.05, and 0% for thresholds of 0.1 and above. At 120 months the fraction is 7.5% for threshold 0.01 and 0.4% for 0.02, and 0% above. By offset $\tau$ the fraction inside at 120 months is 2.9%, 1.1%, 0.5% and 0% for $\tau=0,6,12,24$ months; at 240 months it is 11.1%, 13.2%, 13.4% and 6.8%.

![Figure 10.3. Observed and simulated cumulative incidence at supported horizons](../results/figures/fig_p5_09_horizon_estimates.png)

**Figure 10.3.** Observed and simulated cumulative incidence at the horizons supported by the data. *Setting:* METABRIC Kaplan–Meier estimates at 5, 10, 15 and 20 years with 95% CIs (the 15-year point is drawn as an open square because it is flagged); coloured dots, the reference mapping (threshold 0.1, $\tau=d=0$) for each of the seven priors; the grey bar is the range over all 1,764 scenarios. *Axes/panels:* x, horizon with numbers at risk and events (1,581 and 604; 984 and 819; 521 and 923; 162 and 955); y, cumulative incidence. *What it shows:* at 5 years all dots are at or near zero (largest, P6, 0.036) against 0.263; at 10 years they lie between 0.034 and 0.232 against 0.376; at 15 years between 0.249 and 0.369 against 0.456; at 20 years between 0.317 and 0.468 against 0.505, and only the P4 dot is close to the CI. *Mathematical interpretation:* the gap closes as the horizon grows (the reference P1 difference is −0.262, −0.285, −0.148, −0.069), and the grey range bars include the observed value only at 10 years and later. *Biological/clinical interpretation:* the model's timing is later than the observed composite endpoint. *Limitations:* dots are single reference scenarios, not fitted values; the 15-year and 20-year points are weak clinically.

### 10.4 Alignment diagnostic (exploratory)

The alignment diagnostic shifts the reference simulated curve horizontally by $s\in[-60,+60]$ months and reports the supremum distance to the METABRIC curve over the window where both exist (`p5_timeorigin_alignment_curve.csv`, `p5_timeorigin_alignment_minima.csv`). For each of the seven priors the minimum is at the **edge of the search range, $s=-60$ months**: 0.181 (P1), 0.117 (P2), 0.148 (P3), 0.181 (P4), 0.186 (P5), 0.039 (P6) and 0.231 (P7), against 0.310, 0.275, 0.281, 0.312, 0.311, 0.228 and 0.350 at no shift. The fixed 0–180 month window, a column added to avoid the window shrinking at negative shifts, gives identical values. Because the minimum is at the boundary there is **no interior minimum**, and the diagnostic does not provide an alignment offset. A negative shift moves the simulated events earlier; its only message is that the simulated events are late relative to the observed ones, as found in S1 and S5. Since the offset $\tau$ is non-negative by definition, this diagnostic is not an estimate of the time origin, and it is excluded from the consistency criterion.

### 10.5 Summary of the comparison

![Figure 10.4. Summary of the clinical comparison](../results/figures/fig_p5_10_clinical_comparison_summary.png)

**Figure 10.4.** Summary of the clinical comparison. *Setting:* METABRIC estimates and the P1 reference scenario at 5, 10, 15 and 20 years, and the verdict on the six pre-specified statements. *Axes/panels:* (a) bars for METABRIC (black, with CI caps) and the reference model (blue) with the orange line showing the model range over the P1 threshold, $\tau$, $d$ and death-first grid; (b) text of the statements with the fraction of scenarios in which each holds (S1 and S5 robust; S2 98.8%; S3 11.1%; S4 99.6%; S6 94.3%) and a reading note. *What it shows:* the reference model is far below METABRIC at 5 and 10 years and closer at 15 and 20 years; the P1 range covers the observed value only at 15 and 20 years. *Mathematical interpretation:* the findings are about timing and shape (a deficit of early events and a back-loaded late increment), and they hold across assumptions; level agreement at 20 years does not. *Biological/clinical interpretation:* the figure supports the descriptive statement that this model family, with no treatment or early-recurrence mechanism, does not generate the early recurrences seen in the composite RFS endpoint. It does not support a statement about dormancy as a mechanism. *Limitations:* descriptive consistency only; no calibration; every number depends on the declared grid.

---

## 11. Detection-Threshold Uncertainty

**Table 11.1.** Simulated cumulative incidence under prior P1_base (n = 15,000) for each threshold (source: `results/tables/p5_threshold_sensitivity.csv`). Inside = within the METABRIC 95% CI.

| Threshold | $F_{\mathrm{sim}}(5\,\mathrm y)$ | $F_{\mathrm{sim}}(10\,\mathrm y)$ | $F_{\mathrm{sim}}(15\,\mathrm y)$ | $F_{\mathrm{sim}}(20\,\mathrm y)$ | Median crossing time (events) | Inside CI at 10 y / 20 y |
|---|---|---|---|---|---|---|
| 0.01 | 0.0421 | 0.3595 | 0.5915 | 0.6993 | 118.4 | yes / no |
| 0.02 | 0.0149 | 0.2697 | 0.5211 | 0.6431 | 129.8 | no / no |
| 0.05 | 0.0027 | 0.1615 | 0.4025 | 0.5329 | 143.0 | no / no |
| 0.1 | 0.0005 | 0.0910 | 0.3083 | 0.4363 | 153.5 | no / no |
| 0.2 | 0.0000 | 0.0429 | 0.2173 | 0.3421 | 165.2 | no / no |
| 0.3 | 0.0000 | 0.0247 | 0.1627 | 0.2923 | 174.0 | no / no |
| 0.5 | 0.0000 | 0.0083 | 0.1027 | 0.2181 | 183.3 | no / no |
| METABRIC | 0.2628 | 0.3757 | 0.4559 | 0.5050 | – | – |

For P1 the threshold alone moves $F_{\mathrm{sim}}(20\,\mathrm y)$ from 0.699 to 0.218 and $F_{\mathrm{sim}}(10\,\mathrm y)$ from 0.360 to 0.008, and the conditional median crossing time from 118 to 183 months. At no threshold does P1 get close to the 5-year METABRIC value (the largest is 0.042). At 20 years P1 at threshold 0.05 gives 0.533, just above the CI upper limit 0.532. Two other priors are inside the 20-year CI at one threshold each: P3 (log-uniform) at 0.01 (0.507) and P6 at 0.05 (0.522).

![Figure 11.1. Detection-threshold sensitivity](../results/figures/fig_p5_07_detection_threshold_sensitivity.png)

**Figure 11.1.** Detection-threshold sensitivity of the simulated cumulative incidence. *Setting:* seven priors, seven thresholds (0.01 to 0.5 of $K$), $\tau=d=0$; METABRIC estimate with 95% CI as a horizontal line and grey band. *Axes/panels:* (a), (b), (c) simulated cumulative incidence at 5, 10 and 20 years against the assumed threshold (log scale); (d) P1 event-free curves for each threshold against the METABRIC curve. *What it shows:* every curve decreases monotonically with the threshold; at 5 years all values are far below 0.26 (P6 at 0.01 reaches 0.139); at 10 years only the lowest thresholds of P1, P2, P4, P5 and P6 approach 0.376; at 20 years the curves cross the METABRIC line near thresholds 0.05 to 0.1 (P3 near 0.01). In (d) a lower threshold moves the whole curve left, so the early part of METABRIC is still not reached. *Mathematical interpretation:* the threshold shifts the first-crossing time by about $\ln(\theta_2/\theta_1)/\lambda_+$ in the exponential regime, so it acts like a time shift that depends on the individual growth rate. *Biological/clinical interpretation:* agreement at 20 years can be produced or removed by choosing the threshold, so such agreement cannot be used as evidence. *Limitations:* thresholds are model-defined burden fractions, not imaging limits; priors are hypothetical.

---

## 12. Time-Origin and Observation-Delay Uncertainty

**Exact translation.** Because $\tau$ and $d$ enter only through $\tau+d$, a nuisance offset translates the whole simulated curve to the right by $\tau+d$ months and leaves its shape unchanged (Equation 9.1). For P1 with threshold 0.1 and no death-first share (`p5_timeorigin_sensitivity.csv`):

**Table 12.1.** Effect of the offset $\tau+d$ on the simulated cumulative incidence (P1_base, threshold 0.1, death-first 0).

| $\tau+d$ (months) | $F_{\mathrm{sim}}(5\,\mathrm y)$ | $F_{\mathrm{sim}}(10\,\mathrm y)$ | $F_{\mathrm{sim}}(15\,\mathrm y)$ | $F_{\mathrm{sim}}(20\,\mathrm y)$ |
|---|---|---|---|---|
| 0 | 0.0005 | 0.0910 | 0.3083 | 0.4363 |
| 6 | 0.0001 | 0.0706 | 0.2908 | 0.4259 |
| 12 | 0.0000 | 0.0539 | 0.2707 | 0.4157 |
| 24 | 0.0000 | 0.0266 | 0.2273 | 0.3933 |
| 36 | 0.0000 | 0.0103 | 0.1805 | 0.3697 |

An offset up to 36 months lowers $F_{\mathrm{sim}}(10\,\mathrm y)$ from 0.091 to 0.010 and $F_{\mathrm{sim}}(20\,\mathrm y)$ from 0.436 to 0.370. Because simulated incidence is already below METABRIC at 5 and 10 years, **any non-negative offset moves the model further from the observed curve at these horizons**. At 20 years an offset widens the gap as well (the METABRIC estimate is 0.505), and for P1 at threshold 0.1 no offset in the grid brings the 20-year value inside the CI. The share of scenarios inside the CI at 10 years falls from 2.9% at $\tau=0$ to 0% at $\tau=24$ months, and at 20 years it changes from 11.1% to 13.4% (at $\tau=12$) and 6.8% (at $\tau=24$): the 20-year result depends non-monotonically on the offset because it interacts with threshold and prior.

**Death-first mapping.** Assigning a death-first share of 10.5% or 27.2% with a 24-month mean lag changes P1 (reference) $F_{\mathrm{sim}}(10\,\mathrm y)$ from 0.0910 to 0.0854 and 0.0781, and $F_{\mathrm{sim}}(20\,\mathrm y)$ from 0.4363 to 0.4307 and 0.4243. This is a small effect compared with the threshold and offset.

![Figure 12.1. Time-origin and observation-delay uncertainty](../results/figures/fig_p5_08_time_origin_delay_uncertainty.png)

**Figure 12.1.** Time-origin and observation-delay uncertainty. *Setting:* P1_base at threshold 0.1 with $\tau+d$ from 0 to 36 months; death-first shares 0, 10.5% and 27.2%; alignment diagnostic over seven priors. *Axes/panels:* (a) P1 event-free curves shifted by $\tau+d=0,6,12,24,36$ months against METABRIC; (b) simulated cumulative incidence at 5, 10 and 20 years against $\tau+d$ (0 to 48 months), with METABRIC CI bands; (c) simulated $F(10\,\mathrm y)$ and $F(20\,\mathrm y)$ against the death-first share, with METABRIC values as black lines; (d) supremum distance to METABRIC against the horizontal shift applied to the simulated curve (−60 to +60 months), with minima marked by triangles. *What it shows:* (a) a parallel shift to the right; (b) the 10-year curve falls from 0.09 towards 0, the 20-year curve from 0.44 to about 0.34, and the 5-year curve stays at 0; (c) bars barely change with the death-first share; (d) all seven curves increase with the shift and have their minima at −60 months, the edge of the range. *Mathematical interpretation:* (a)–(c) confirm exact translation and the small effect of the death-first mapping; (d) is a descriptive diagnostic with no interior minimum, so it contains no estimate of $\tau$ (Section 10.4). *Biological/clinical interpretation:* since the RFS origin and any detection delay are undocumented, the data cannot tell whether the simulated clock is early or late; they only show that no non-negative offset removes the early deficit. *Limitations:* the offset grid ends at $\tau+d=36$ months, the lag distribution is a single assumption, and panel (d) is exploratory.
