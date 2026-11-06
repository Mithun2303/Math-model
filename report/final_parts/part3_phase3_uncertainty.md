## 5. Phase 3: Sensitivity, Basins of Attraction and Monte Carlo Relapse Simulation

*Everything in this section describes a **hypothetical** mathematical model, **illustrative** parameter ranges and **assumed** uncertainty distributions. "Relapse" means the first upward crossing of a model-defined burden threshold. No simulated probability is a patient risk, a clinical prediction or a calibrated estimate, and no threshold is a clinical imaging limit.*

**Why this phase exists.** Phase 2 established the exact structure of the P–D–E model at one baseline parameter set. That baseline is illustrative, so the natural next questions are: (i) which parameters actually control relapse occurrence and timing when all of them are uncertain; (ii) can a locally stable tumour-free state still be left by a large enough disturbance (basins of attraction); (iii) what cumulative relapse probabilities does the model generate when parameters and the initial residual burden are drawn from stated distributions; and (iv) how much of that output is driven by numerics, by finite sampling, by the choice of prior, or by the choice of detection threshold. These map to objectives O3 and O4 of the research plan (sensitivity analysis and relapse-time distributions). Design file: `configs/phase3.yaml` (sha256[:12] `c7e70720ebfc`, printed in every Phase 3 figure footer).

---

### 5.1 Design principles and definitions

The following choices were fixed before any result was inspected.

1. **Gauge.** All runs use the normalised model of Phase 2 (Section 4): tumour-free immune level $e=1$ and gauge-invariant killing $\beta = ks/\mu$. Fixed: $K=1$, $\mu=0.05$, $h=0.02$. Initial immune level $e_0=1$ unless a basin experiment varies it.
2. **One trajectory runner.** `uncertainty.run_trajectory` integrates once, densely and without terminal events on $[0,240]$ months (DOP853, rtol $10^{-9}$, atol $10^{-12}$), and derives from that path the first upward crossing for every threshold, $N(60/120/240)$, downward crossings, and the dormancy-duration metric. It was cross-checked against the previously verified terminal-event solver on every Monte Carlo sample: 5,000/5,000 agree, max $|\Delta T| = 8\times10^{-11}$ months.
3. **Long-horizon outcome kept separate.** A second classifier integrates to $T_\text{long}=\max(10^5, 50/\text{slowest stable rate})$ months and records which attractor is reached. A label is accepted only if the end state is within $10^{-3}$ (relative) of an enumerated stable equilibrium; otherwise it is "unresolved" (count: 0 in every experiment). The right-hand side is clipped to the non-negative orthant *only inside this classifier*, to stay away from the pole of $P/(h+P)$ at $P=-h$.
4. **Samples drawn up front** from seeded generators (Latin hypercube, scrambled Sobol, Gaussian copula) and mapped over a process pool, so results do not depend on scheduling (tested: parallel = serial exactly).
5. **Start-above-threshold trap.** The legacy solver returns "event at $t=0$" when $N_0\ge N_\text{detect}$. In Phase 3 such samples have *no defined event* and are excluded and counted. Supports are chosen so $N_0\le5\times10^{-3}$ is below every threshold $\ge0.01$; the excluded count was 0 in all runs.

**Event definition.** For threshold $\theta$, relapse is the first time $N=P+D$ crosses $\theta$ upward from below, right-censored at the horizon (240 months by default; 60/120/240 months are reported as 5/10/20 years). After the first crossing an exceedance is **transient** if $N$ crosses back below $\theta$ within the horizon and **persistent** otherwise; a trajectory that later settles at an attractor below $\theta$ is reported through the long-run attractor. Crossings are bracketed on a 0.05-month grid and refined by Brent's method on the dense interpolant (an excursion shorter than 0.05 month would be missed).

**Operational dormancy duration (hypothetical descriptor).** Total time before the first crossing (or horizon) during which $N<\theta$, $|d\ln N/dt|<0.01\ \text{month}^{-1}$ and $D/N\ge0.5$. At the baseline reference it is 18.1 months, and it depends on the definition's own cut-offs (10.2 months for $|d\ln N/dt|<0.005$, 32.9 for 0.02; the dormant-fraction cut-off hardly matters; `p3_dormancy_definition_sensitivity.csv`). It is a descriptor of the transient, not a measurement.

---

### 5.2 Local sensitivity

**Method.** For each of ten inputs ($r,q,a,\delta_D,\beta,\rho,\eta,P_0,D_0,N_\text{detect}$) central differences in log-parameter at steps of 5% and 2.5% give elasticities $d\ln(\text{output})/d\ln\theta$. The finer estimate is reported **only if** the two agree within 20% and every stencil point defines the output; otherwise the entry is *undefined* (e.g. the event is absent at one stencil point, or an attractor switches) or *unstable*. Event *indicators* are never differentiated: for "event by 240 months" the smallest multiplicative change (within a factor 10, either direction) that flips the indicator is found by scan plus bisection. Three reference points were used: $R_\text{early}$ ($\beta=0.020$, first crossing at 149 months), $R_\text{baseline}$ ($\beta=0.040$, 219.9 months, close to the horizon) and $R_\text{late}$ ($\beta=0.0567$, no event by 240 months, about 375 months in a longer run). Stencil outcomes: 167 ok, 48 zero, 22 undefined, 3 unstable.

**Table 5.1.** Elasticities at the baseline reference $R_\text{baseline}$ (source: `results/tables/p3_local_sensitivity.csv`).

| input | $d\ln N(120)$ | $d\ln N(240)$ | $d\ln T_\text{first crossing}$ | $d\lambda_{TF}/d\ln\theta$ [1/month] |
|---|---|---|---|---|
| $r$ | +10.3 | +9.9 | undefined | +0.111 |
| $\beta$ | −3.46 | −4.07 | +0.99 | −0.037 |
| $q$ | −2.72 | −2.47 | +0.66 | −0.033 |
| $\rho$ | −0.03 | −1.07 | +0.18 | 0 |
| $a$ | +0.99 | +0.54 | −0.15 | +0.003 |
| $\eta$ | 0.00 | +0.21 | −0.03 | 0 |
| $\delta_D$ | −0.03 | −0.03 | +0.007 | 0 |
| $D_0$ | +0.88 | +0.31 | −0.09 | 0 |
| $P_0$ | +0.09 | +0.03 | −0.009 | 0 |
| $N_\text{detect}$ | 0 (exact) | 0 (exact) | +0.23 | 0 |

![Figure 5.1. Local sensitivity tornado plots](../results/figures/fig_p3_01_local_sensitivity_tornado.png)

**Figure 5.1.** Local sensitivity tornado plots. *Setting:* central log-derivatives (steps 5% and 2.5% must agree within 20%) at three reference points, all other inputs at baseline. *Axes/panels:* top row, elasticity of $N(240)$ for $R_\text{early}$, $R_\text{baseline}$, $R_\text{late}$; bottom row, elasticity of the first-crossing time $T$ for the same three points. Bars (orange positive, blue negative) are dimensionless elasticities; italic labels mark undefined entries. *What it shows:* $r$ dominates $N(240)$ at every reference (about +2.8, +9.9, +15.6 moving from early to late), followed by $\beta$ and $q$ with negative signs. For $T$ at $R_\text{baseline}$ the bar for $r$ is "undefined"; at $R_\text{late}$ the whole panel is undefined because no event exists. *Mathematical interpretation:* $N(t)$ is proportional to $e^{\lambda_+ t}$ with $\lambda_+$ depending on $r-\beta-q$, so 240 months spans many e-folds and elasticities grow as the reference approaches the stability boundary. *Biological interpretation:* in the model, net proliferation outweighs every other process. *Limitations:* local, three reference points, and meaningless where the event disappears inside the stencil.

![Figure 5.2. Local sensitivity map](../results/figures/fig_p3_01b_local_sensitivity_map.png)

**Figure 5.2.** Local sensitivity map. *Setting:* the same three references; eight output columns. *Axes/panels:* one heat map per reference; rows are the ten inputs; columns are $N(60)$, $N(120)$, $N(240)$, first-crossing time $T$, equilibrium burden $N^*$, dormancy duration, $\lambda_{TF}$ and $\lambda_\text{att}$ (the last two as $d\lambda/d\ln\theta$ in 1/month); colour is scaled within each column, "n/d" = undefined and "~" = unstable. *What it shows:* the $r$ row is dark red for burdens (e.g. 10 for $N(120)$ and 9.9 for $N(240)$ at baseline, 16 at $R_\text{late}$); $\beta$ gives −3.5 and −4.1; $D_0$ elasticity (0.92–0.95 for $N(60)$) is much larger than $P_0$ (0.05–0.11). The $N_\text{detect}$ row is exactly 0 for all burdens and non-zero only for $T$ (0.15, 0.23). $N^*$ depends on $\beta$ (−0.39 at baseline) and $\eta$ (+0.3). *Mathematical interpretation:* a start that is 99% dormant puts most of the initial burden on $D_0$, which then dominates the projection on the growing eigenmode (Phase 2). The detection threshold cannot change the trajectory, only the time at which it is intercepted. *Biological interpretation:* the composition of residual disease matters locally, the detection rule does not alter biology. *Limitations:* a local statement; the global analysis shows composition becomes unimportant once total burden is varied.

![Figure 5.3. Flip margins of the 240-month event](../results/figures/fig_p3_01c_local_flip_margins.png)

**Figure 5.3.** Distance of the baseline reference to the 240-month event boundary. *Setting:* $R_\text{baseline}$; for each input the smallest relative change that flips "event by 240 months". *Axes:* x = input; y = size of the smallest relative change [%] (orange: increase, blue: decrease); "no flip (×10)" marks inputs for which no change within a factor 10 flips the outcome. *What it shows:* $r$ −3%, $\beta$ +8%, $q$ +13%, $\rho$ +33%, $N_\text{detect}$ +41%, $a$ −48%, $D_0$ −65%; $\delta_D$, $\eta$, $P_0$ do not flip. *Mathematical interpretation:* the event indicator is a discontinuous function of the inputs, so a derivative is not defined across the boundary; the margin is the honest local measure. *Biological interpretation:* a reference whose first crossing lies near the 20-year horizon is fragile: a 3.3% reduction in net proliferation removes the simulated event within 240 months. *Limitations:* margins are one-at-a-time, and the baseline reference was deliberately chosen close to the horizon so it is a pessimistic case for robustness; $R_\text{early}$ is more robust ($\beta$ needs ×2.17, $r$ ×0.78).

---

### 5.3 Global sensitivity analysis

**Design.** Ten inputs on hypothetical ranges (roughly baseline ×[1/5, 5] except initial burden): $r\in[0.03,0.30]$, $q\in[0.01,0.16]$, $a\in[0.001,0.025]$, $\delta_D\in[2\times10^{-4},0.01]$, $\beta\in[0.01,0.20]$, $\rho\in[0.01,0.16]$, $\eta\in[0.02,0.32]$ (log-uniform), $N_0\in[10^{-6},5\times10^{-3}]$ (log-uniform), dormant fraction $f_D=D_0/N_0\in[0.5,0.999]$ (logit-uniform), $N_\text{detect}\in[0.02,0.30]$ (log-uniform). Sobol indices of an indicator are only meaningful if its variance is not degenerate; the event-by-240 probability over the design is 0.304 (0.300–0.306 for $n=250$–4000), inside the required [0.2, 0.8]. Regime composition ($n=4000$): nonzero only 63.5%, tumour-free only 32.4%, two nonzero 2.8%, tumour-free + nonzero 1.4%; samples with $r_\text{eff}\le0$ are a regime, not failures; solver failures 0.

**Outputs analysed:** $\log_{10}N(240)$; restricted time $\min(T,240)$; the event indicator; $\lambda_{TF}$; highest stable equilibrium burden. Sobol indices are *not* computed on $T$ conditional on an event (the conditioning set would change with the parameters).

**Estimators.** Spearman and PRCC with 400-resample bootstrap intervals ($n=4000$ LHS); Morris screening (40 trajectories, 4 levels, 440 runs); Sobol with Saltelli (2010) first-order, Jansen total-order and Saltelli (2002) second-order estimators, $N=2048$ base rows (45,056 runs), 200 bootstrap resamples. The estimators were validated on the Ishigami function (analytic $S_1=(0.314,0.444,0)$, $S_T=(0.557,0.442,0.244)$, $S_{13}=0.245$), reproduced in a unit test.

**Table 5.2.** Global sensitivity of the main outputs (source: `results/tables/p3_global_rank_correlations.csv`, `p3_global_sobol.csv`).

| input | $\log_{10}N(240)$ Spearman / PRCC | $\min(T,240)$ Spearman / PRCC | $\min(T,240)$ Sobol $S_1$ / $S_T$ | event Sobol $S_1$ / $S_T$ |
|---|---|---|---|---|
| $r$ | +0.67 / +0.83 | −0.66 / −0.73 | 0.53 / 0.78 | 0.42 / 0.78 |
| $\beta$ | −0.42 / −0.69 | +0.33 / +0.49 | 0.12 / 0.33 | 0.13 / 0.43 |
| $q$ | −0.18 / −0.38 | +0.17 / +0.28 | 0.01 / 0.12 | 0.02 / 0.20 |
| $\rho$ | −0.08 / −0.19 | +0.11 / +0.18 | 0.03 / 0.08 | 0.02 / 0.16 |
| $N_0$ | +0.32 / +0.58 | −0.09 / −0.15 | 0.01 / 0.06 | 0.01 / 0.10 |
| $N_\text{detect}$ | 0.00 / −0.02 | +0.09 / +0.17 | 0.02 / 0.08 | 0.01 / 0.11 |
| $a$ | +0.05 / +0.11 | −0.05 / −0.08 | ≈0 / 0.01 | ≈0 / 0.05 |
| $\eta,\delta_D,f_D$ | $\lvert\rho_s\rvert\le0.06$ | ≤0.03 | $S_T\le0.01$ | $S_T\le0.03$ |

![Figure 5.4. Global sensitivity ranking](../results/figures/fig_p3_02_global_sensitivity_ranking.png)

**Figure 5.4.** Global sensitivity ranking. *Setting:* ten hypothetical inputs, $n=4000$ LHS for rank statistics and Sobol $N=2048$ (45,056 runs). *Axes/panels:* three columns (outputs $\log_{10}N(240)$, restricted time $\min(T,240)$, event by 240 months). Top row: Spearman (blue) and PRCC (orange) rank correlations with 95% bootstrap bars. Bottom row: Sobol first-order $S_1$ (green) and total $S_T$ (purple) with bootstrap intervals; inputs sorted by importance. *What it shows:* $r$ is first and $\beta$ second in every panel (restricted time $S_T$ 0.78 and 0.33; event 0.78 and 0.43); then $q$, $\rho$, $N_0$, $N_\text{detect}$; $a$, $\eta$, $\delta_D$, $f_D$ are close to zero. *Mathematical interpretation:* PRCC magnitudes exceed Spearman because they partial out the dominant $r$; large gaps $S_T-S_1$ (about 0.25 for $r$ and 0.21 for $\beta$ in restricted time, 0.36 and 0.30 for the event) indicate interactions; wide $S_1$ intervals (e.g. $S_1(r)=0.53$, 95% interval [0.39, 0.69]) are visible. *Biological interpretation:* within these ranges, timing is controlled by net proliferation and immune killing; reactivation, dormant loss, immune suppression and the *composition* of the initial burden barely matter. *Limitations:* indices describe the chosen ranges (those for $r$ span ×10, $\beta$ ×20, $N_0$ ×5000), not real uncertainty; PRCC assumes monotone responses.

![Figure 5.5. Binned scatter against the four most influential inputs](../results/figures/fig_p3_02b_global_scatter_top_inputs.png)

**Figure 5.5.** Output against the four most influential inputs. *Setting:* $n=4000$ LHS samples. *Axes/panels:* columns $r$, $\beta$, $q$, $\rho$ (log x-axes); top row $\log_{10}N(240)$, bottom row restricted time $\min(T,240)$; grey dots are samples, blue line the binned median, shaded band the 10–90% range. *What it shows:* the median of $\log_{10}N(240)$ rises steeply with $r$ (from about −4.4 to near 0), falls with $\beta$ (about −1.3 to −4.8), falls weakly with $q$ and is nearly flat in $\rho$; every band is very wide (several decades). In the restricted-time row the median sits on the 240-month ceiling except at large $r$ (it reaches about 55 months at the largest $r$), and for $\beta$, $q$, $\rho$ it is flat at 240. *Mathematical interpretation:* the outcome is thresholded and right-censored, so rank correlations rank the inputs but do not describe a smooth response; the ceiling at 240 is the censoring limit. *Biological interpretation:* faster proliferation reliably produces early relapse, but any one parameter alone leaves the outcome highly uncertain because the others vary simultaneously. *Limitations:* marginal view only; interactions are hidden in the width of the bands.

![Figure 5.6. Convergence of the sensitivity indices](../results/figures/fig_p3_02c_global_sensitivity_convergence.png)

**Figure 5.6.** Convergence of the sensitivity analysis. *Setting:* nested scrambled-Sobol subsets $N=128\ldots2048$; independent LHS at each $n$. *Axes/panels:* left and centre, total-order index $S_T$ of each input versus Sobol base size $N$ (restricted time; event indicator); right, Kendall $\tau$ between the $|$Spearman$|$ ranking at sample size $n$ and the ranking at $n=4000$, for the three outputs. *What it shows:* $S_T(r)$ stays between about 0.76 and 0.81 for restricted time (0.78–0.86 for the event), $S_T(\beta)$ within 0.27–0.33 from $N=256$; lower curves ($q,\rho,N_0,N_\text{detect}$) are flat and well below 0.2. Kendall $\tau$ is 0.82 at $n=2000$ and dips to about 0.64–0.82 at $n=1000$. *Mathematical interpretation:* the top inputs are stable; the near-zero inputs reorder among themselves, which lowers $\tau$ without changing the conclusion. $S_1$ converges more slowly than $S_T$ (e.g. $S_1(\beta)$ 0.22→0.12). *Biological interpretation:* none beyond Figure 5.4. *Limitations:* $S_T$ is stable only to about ±0.03 and the experiment is not repeated with independent designs.

![Figure 5.7. Morris screening](../results/figures/fig_p3_02d_morris_screening.png)

**Figure 5.7.** Morris elementary-effect screening. *Setting:* 40 trajectories, 4 levels, 440 runs. *Axes/panels:* for each output ($\log_{10}N(240)$, $\min(T,240)$, event) $\mu^*$ (mean absolute elementary effect) on x and $\sigma$ (spread, an indicator of interaction or non-linearity) on y; dotted diagonal $\sigma=\mu^*$. *What it shows:* $r$ has the largest $\mu^*$ (about 3.7 for $\log_{10}N(240)$, 156 months for restricted time, 0.87 for the event), then $\beta$, then $q$, $N_0$ and $N_\text{detect}$; $\eta$, $\delta_D$ lie at the origin. Most points are close to or above the diagonal (strongly so for the event indicator), meaning effects are non-linear or interaction-dominated. *Mathematical interpretation:* an independent screening method reproduces the Sobol top-three ordering, which supports it against estimator artefacts. *Biological interpretation:* the same as Figure 5.4. *Limitations:* Morris is a screening tool (few runs), not a variance decomposition.

![Figure 5.8. Second-order interaction heat map](../results/figures/fig_p3_03_parameter_interaction_heatmap.png)

**Figure 5.8.** Parameter interactions. *Setting:* Saltelli (2002) second-order estimator, $N=2048$; cells whose 95% bootstrap interval contains 0 are blank (estimator noise). *Axes/panels:* three heat maps (outputs as before); off-diagonal cells show $S_{ij}$; the diagonal shows $S_T-S_1$ (interaction share of each input). *What it shows:* the diagonals are non-zero for $r$ (0.15, 0.25, 0.36) and $\beta$ (0.12, 0.21, 0.30) and for $q$ (0.06, 0.10, 0.18); **only one off-diagonal pair is resolved: $r\times\beta$ for the event indicator, $S_{r\beta}\approx0.10$**. *Mathematical interpretation:* the interaction variance is large but spread over many pairs and/or higher orders; the claim "$r$ and $\beta$ interact" rests on $S_T-S_1$ and the single resolved pair, not on a dominant two-factor term. An earlier unmasked version showed a uniform 0.04 row for $q$; that was estimator noise and is not reported. *Biological interpretation:* the effect of immune killing depends on how fast the tumour grows (killing matters when $\beta$ is comparable with $r_\text{eff}$). *Limitations:* second-order estimators have wide intervals; third-order effects are unmeasured.

**Interpretation of §5.3.** Over these ranges, simulated relapse is controlled by net proliferation $r$ and immune killing $\beta$ first, then dormancy entry $q$, immune stimulation $\rho$, detection threshold and total initial burden $N_0$. The AUC of the event indicator for $r$ is 0.90 ($\beta$: 0.29, i.e. low $\beta$ ⇒ event). Among events only, crossing time correlates with $r$ (−0.61) and $N_0$ (−0.35) but not with $\beta$ (+0.05, interval contains 0): conditioning on an event removes much of $\beta$'s role. **Sensitivity is not uncertainty**: these indices say which inputs *would* matter if they were uncertain over these ranges; they say nothing about how uncertain any input actually is.

---

### 5.4 Basins of attraction

Phase 2 showed that two kinds of bistability exist. Phase 3 asks how large the basins are and whether "below detection threshold at the start" implies "tumour-free at the end". Three hypothetical parameter sets were used, each verified by Phase 2 equilibrium enumeration and eigenvalues.

**Table 5.3.** Basin-analysis cases (source: `report/phase_3_sensitivity_and_uncertainty.md` §4; `results/tables/p3_basin_*.csv`).

| case | parameters (others baseline) | stable states | saddle burden |
|---|---|---|---|
| A: tumour-free + nonzero | $\eta=0.6$, $\beta=0.2288$ ($k=0.572$) | TF; NZ at $N^*=0.737$ | 0.1495 |
| B1: two nonzero, baseline-adjacent | $\eta=0.164$, $\beta=0.108$ ($k=0.27$) | NZ-low 0.0199 (focus); NZ-high 0.502 (node); TF is a saddle | 0.2085 |
| B2: two nonzero, well separated | $a=0.306$, $q=1.58\times10^{-3}$, $\rho=0.240$, $\eta=1.637$, $\beta=0.060$ | NZ-low 0.0028; NZ-high 0.982; TF saddle | 0.111 |

**Method.** For each case two 60×60 slices of the 3-D initial-condition space: $(P_0,D_0)$ log-spaced in $[10^{-6},0.9]$ at $e_0=1$, and $(N_0,e_0)$ with $N_0$ log-spaced, $e_0\in[0.05,1.6]$ and the dormant fraction at its equilibrium value $q/(q+a+\delta_D)$. Each initial state is classified by the long-horizon attractor *and* simulated for 240 months for the event outcome at $N_\text{detect}=0.1$ plus extra thresholds. 3,000 random 3-D starts per case give attractor frequencies. Boundaries are located by bisection along lines (bracket $2^{-40}$).

**Numerical verification.** Boundary cells (120–203 per slice) and 200 interior cells per slice were reclassified with tighter numerics (DOP853 $10^{-12}/10^{-15}$, Radau, horizons ×5): label agreement 100% in all six slice-by-case combinations, 0 unresolved; a $10^{-6}$ relative jitter of the initial state flipped no labels. (Boundary cells are cells adjacent to a label change on the grid, not points on the separatrix.)

![Figure 5.9. Basins: tumour-free + nonzero bistability (case A)](../results/figures/fig_p3_04_basins_tumour_free_nonzero.png)

**Figure 5.9.** Basins of attraction for tumour-free + nonzero bistability, $\eta=0.6$ (case A). *Setting:* 60×60 grids; attractor from a long-horizon run. *Axes/panels:* left, slice $(P_0,D_0)$ at $e_0=1$ (log–log); right, slice $(N_0,e_0)$ with equilibrium dormant fraction. Blue = basin of tumour-free (TF), orange = basin of nonzero (NZ); stars are stable equilibria, the crossed marker is the saddle, dots mark boundary cells, the dashed curve/line is $N_0=0.1$ (the detection threshold), circles/triangles mark starts with a 240-month event at $N_\text{detect}=0.3$/0.5. *What it shows:* the NZ basin (21.6% of the left slice, 13.5% of the right) starts only at $N_0\gtrsim0.15$; the region below the dashed line is entirely TF. Bisection puts the boundary at $N_0=0.165$ at $e_0=1$ and 0.145 at $e_0=0.3$ (saddle burden 0.1495). *Mathematical interpretation:* a locally stable tumour-free state has a *bounded* basin; the saddle's own immune level and composition differ from the slice, so the boundary need not equal the saddle burden. *Biological interpretation:* in this regime control holds from every below-threshold residual state, but a large re-seeded burden escapes. With a higher threshold ($\theta=0.5$), 290 of 389 below-threshold starts in the NZ basin cross it, so the event outcome depends on where the threshold lies relative to the saddle. *Limitations:* 2-D slices of a 3-D space; fractions are properties of the chosen box, not probabilities.

![Figure 5.10. Basins: two nonzero attractors (cases B1, B2)](../results/figures/fig_p3_05_basins_two_nonzero.png)

**Figure 5.10.** Basins for two stable nonzero states. *Setting:* top row case B1 ($\eta=0.164$), bottom row case B2 (strong immune stimulation); same grid and method as Figure 5.9. *Axes/panels:* left column $(P_0,D_0)$ at $e_0=1$; right column $(N_0,e_0)$. Green = basin of the low-burden state (NZ-low), purple = basin of the high-burden state (NZ-high); circles/triangles = 240-month events at the extra thresholds (0.3, 0.5 for B1; 0.01, 0.3 for B2). *What it shows:* in B1 the NZ-high basin needs $N_0\gtrsim0.2$ (bisection 0.216 at $e_0=1$, 0.181 at $e_0=0.3$; saddle 0.2085; 21.1% / 11.5% of the two slices). In B2 the boundary is lower and slanted: 0.071 at $e_0=1$ and **0.039 at $e_0=0.3$** (saddle 0.111), so immune depletion pulls small residual states into the high basin; 35.5% / 20.6% of the slices are NZ-high, and 8.5% of random starts with $N_0<0.1$ reach the high attractor. In the lower-right panel a wedge of event markers (threshold 0.01) lies *inside* the NZ-low basin. *Mathematical interpretation:* basin boundaries move with $e_0$ and need not coincide with equilibrium burdens. *Biological interpretation:* when immune stimulation is strong, a depleted immune system can convert a sub-detection residual burden into a high-burden relapse. *Limitations:* three hypothetical parameter sets; B2 is an extreme corner.

![Figure 5.11. Event outcome versus attractor reached](../results/figures/fig_p3_05b_basin_attractor_vs_event.png)

**Figure 5.11.** Event outcome versus attractor reached. *Setting:* 3,000 random 3-D initial states per case, restricted to starts strictly below each threshold; groups with $n<10$ omitted. *Axes/panels:* three panels (A, B1, B2); x = detection threshold, y = fraction with an event by 240 months, bars coloured by the attractor reached; $n$ printed on each bar. *What it shows:* case A: starts that end at TF essentially never register an event (0.1: $n=2071$; the NZ group at 0.1 has $n=1$, omitted), whereas for threshold 0.3 and 0.5, 0.8 and 0.75 of the NZ-basin starts cross ($n=192$ and 389). B1: starts ending in NZ-low register events at 0.1 (a few percent, $n=2071$; 19 downward crossings and 32 still above $\theta$ at 240 months), 0.38 of the NZ-high group crosses at 0.3 ($n=149$). B2: **all** NZ-high starts below 0.1 (177) and 0.3 (585) register an event (fraction 1.0), and about 0.2 of the NZ-low starts cross 0.01 ($n=1325$) and fall back. *Mathematical interpretation:* event and attractor differ in both directions: a transient excursion can cross a threshold without a change of basin, and a start in the high basin only counts if it starts below $\theta$. *Biological interpretation:* a "relapse" defined by threshold crossing may therefore be a transient overshoot rather than a state change. *Limitations:* the fraction depends on the sampling box.

**Findings summarised.** (1) A stable tumour-free state does not guarantee control from every initial condition, but in the baseline-adjacent case the nonzero basin lies above any plausible residual-disease start. (2) Strong immune stimulation (B2) allows basin selection from below-threshold starts, and the basin boundary falls with immune depletion. (3) Event outcome ≠ attractor. (4) Bistable regimes are rare under the main prior P1 (1.8% of samples), so the prior hardly stresses basin effects. A transient-exceedance fixture found numerically (B2, $N_0=0.0328$, $f_D=0.493$, $e_0=0.606$, $\theta=0.05$: first crossing at 8.75 months, downward crossing at 18.68 months, $N(240)=0.0022$) is used as a unit test.

---

### 5.5 Monte Carlo relapse simulation

#### 5.5.1 Uncertainty distributions (all hypothetical)

**P1_base (primary):** $r=0.12$, $\mu$, $h$ fixed; log-normal with stated median and 95% factor $f$ ($\sigma=\ln f/1.96$): $q$ (0.04, $f=3$), $a$ (0.005, 5), $\delta_D$ (0.001, 5), $\beta$ (0.04, 3), $\rho$ (0.04, 3), $\eta$ (0.08, 3); $N_0$ log-normal (median $10^{-4}$, $f=30$) truncated to $[10^{-7},5\times10^{-3}]$; $f_D\sim\text{Beta}(10,1)$. $N_\text{detect}=0.1$, $e_0=1$. **Alternatives** (model-assumption scenarios): P2 wide ($\sigma\times1.5$, $f_D\sim\text{Beta}(4,1)$); P3 log-uniform over the global ranges; P4/P5 Gaussian-copula rank correlation ±0.5 on $(q,a)$ and $(\beta,\rho)$; P6 $r$ also uncertain (log-normal, median 0.12, $f=2$); P7 initial state fixed at baseline (parameter uncertainty only). Regime composition under P1: tumour-free unstable in 95.9% of draws; nonzero only 94.4%, tumour-free only 3.8%, two nonzero 1.6%, tumour-free + nonzero 0.2%.

**Estimands (all simulated, $N_\text{detect}=0.1$):** $P(\text{first crossing by }60/120/240\text{ months})$ with 95% Wilson and Clopper–Pearson intervals; one-sided upper bounds for rare counts; the empirical cumulative curve (equal to Kaplan–Meier since the horizon is the only censoring); restricted mean event-free time $\text{RMST}(\tau)=E[\min(T,\tau)]$; conditional median among events. The unconditional median is **not reported when $P(\text{event by }240)<0.5$** (here 0.424).

#### 5.5.2 Results for P1_base

Three independent runs of the same prior (different seeds) exist: run A (seed `mc_base`, $n=5000$), run B (seed `mc_assumptions`, $n=5000$), run C (five batches of 1000). The headline pools all three ($n=15{,}000$) with Clopper–Pearson intervals, because one run is too small for the 5-year event (0, 4, 4 crossings).

**Table 5.4.** Simulated event probabilities under prior P1 (hypothetical; source: `results/tables/p3_mc_P1_pooled.csv`).

| quantity ($N_\text{detect}=0.1$) | run A | run B | run C | **pooled (15,000)** |
|---|---|---|---|---|
| $P$(by 60 mo) | 0/5000 | 4/5000 | 4/5000 | **8/15,000 = 0.00053 [0.00023, 0.00105]** |
| $P$(by 120 mo) | 0.0854 | 0.0966 | 0.0910 | **0.0910 [0.0864, 0.0957]** |
| $P$(by 240 mo) | 0.4240 | 0.4418 | 0.4432 | **0.4363 [0.4284, 0.4443]** |
| earliest crossing | 61.6 mo | 55.4 mo | 51.9 mo | 51.9 mo |
| conditional median $T$ | 153.2 | 154.1 | 153.0 | **153.5 mo** (IQR 124.8–186.7; 6,545 events) |
| RMST(240) | 204.2 | 203.3 | 202.9 | 203.5 mo |

Homogeneity of the three runs: $\chi^2$ $p=0.15$ (120 mo), $p=0.10$ (240 mo). The 60-month counts (0, 4, 4) are compatible with a single rate 0.00053 (probability of 0 events in 5000 at that rate ≈ 7%), so the one-sided bound 0.0006 from run A alone is **not** a reliable bound for the 5-year probability. Of run A's 2,120 events, 29 (1.4%) are transient within 240 months and 28 (1.3%) end at an attractor below the threshold. 2,880 of 5,000 samples are censored at 240 months; a long-horizon run shows 1,873 of them reach the threshold later and 1,007 never do, so the long-horizon event probability is 0.799 versus 0.424 by 240 months: event-free at 240 months is not relapse-free. Counting only *persistent* exceedance changes $P(240)$ from 0.4418 to 0.4370 (P1) and by at most −0.013 across the seven priors, so conclusions do not depend on this choice.

![Figure 5.12. Distribution of simulated relapse times](../results/figures/fig_p3_06_relapse_time_distribution.png)

**Figure 5.12.** Distribution of simulated first-crossing times. *Setting:* prior P1_base, run A, $n=5000$, $N_\text{detect}=0.1$. *Axes/panels:* left, histogram of crossing times among the 2,120 events (conditional median 153 months marked); centre, probability mass per 8-month bin for all 5,000 simulations with a grey bar for "no event (censored at 240)"; right, density of crossing times for the two regimes with events (NZ+NZ, $n=79$, 14 events; NZ only, $n=4718$, 2,106 events). *What it shows:* events are unimodal and right-skewed, beginning at about 60 months and peaking at 120–160 months; the censored bar holds 0.576 of the mass (2,880/5,000); the two-nonzero regime has few events, which cluster earlier (about 65–175 months) than those of the nonzero-only regime (up to 240). *Mathematical interpretation:* the wide spread comes from log-normal parameter uncertainty mapped through $T\propto1/\lambda_+$, and the censored mass is the part of the distribution beyond 240 months. *Biological interpretation:* within the model, a late-relapse tail arises with no mechanism other than slow growth of a small residual burden. *Limitations:* the histogram of events alone is conditional and must not be read as the distribution of time to relapse for the population.

![Figure 5.13. Simulated cumulative relapse probability under seven priors](../results/figures/fig_p3_07_cumulative_relapse_probability.png)

**Figure 5.13.** Cumulative first-crossing probability $\hat F(t)$. *Setting:* seven priors, 5000 samples each, $N_\text{detect}=0.1$; the thick black line is P1_base (run A) with a 95% Wilson band (Monte Carlo error only). *Axes:* x = time [months] (dotted verticals at 60, 120, 240); y = simulated cumulative probability. *What it shows:* all curves are ≈0 before about 60 months and rise to 0.32–0.47 at 240 months; P6 (r uncertain) rises earliest, P7 (fixed initial burden) latest, P2 (wide) flattens after 150 months, P3 (log-uniform) ends lowest (≈0.32), P4 (positive correlation) ends highest (≈0.47). *Mathematical interpretation:* the shaded band is narrow relative to the separation between scenarios: assumptions dominate sampling error. *Biological interpretation:* none beyond "the model produces late relapses at 5–20 years under plausible-looking priors"; this is not a prediction. *Limitations:* each curve is a deterministic-model mixture over a hypothetical prior; the band says nothing about assumption uncertainty.

![Figure 5.14. Event probabilities at 5, 10 and 20 years](../results/figures/fig_p3_08_event_probabilities_5_10_20y.png)

**Figure 5.14.** Simulated crossing probabilities at 5, 10, 20 years for the seven priors. *Setting:* $n=5000$ each, bars with 95% Wilson intervals. *Axes/panels:* three panels (60, 120, 240 months); x = prior P1–P7; y = simulated probability. *What it shows:* at 60 months the values range from 0.000 (P4, P7) to 0.036 (P6; P1 0.001, P2 0.006, P3 0.005, P5 0.001); at 120 months 0.034 (P7) to 0.232 (P6) (P1 0.097, P2 0.151, P3 0.120, P4 0.096, P5 0.090); at 240 months 0.317 (P3) to 0.468 (P4) (P1 0.442, P2 0.381, P5 0.419, P6 0.452, P7 0.427). *Mathematical interpretation:* the scenario spread (0.15 at 240 months) is 5.5× the typical Monte Carlo interval width (0.027); at 120 and 60 months the ratios are 12× and 19×. Fixing the initial condition (P7) cuts the 10-year value from 0.097 to 0.034, so uncertainty in the initial residual burden generates the early events. *Biological interpretation:* not a risk estimate. *Limitations:* the P1 bar here is run B (0.442), not the pooled 0.436.

![Figure 5.15. Monte Carlo convergence](../results/figures/fig_p3_10_monte_carlo_convergence.png)

**Figure 5.15.** Monte Carlo convergence. *Setting:* prefix sub-samples of one seeded P1 run (run A). *Axes/panels:* running estimate of $P(\text{crossing})$ with 95% Wilson band versus sample size $n$ for 60, 120 and 240 months; orange dots are the independent batches of $n=1000$. *What it shows:* the 240-month estimate stabilises by $n\approx2000$ (0.424 ± 0.014 at $n=5000$); the 120-month band narrows from about ±0.05 to ±0.008; at 60 months there are 0 events at $n=5000$ and the plotted band is only an upper limit (one-sided bound 0.0006). The five batch estimates scatter 0.423–0.460 at 240 months and 0.079–0.099 at 120. *Mathematical interpretation:* binomial $1/\sqrt n$ convergence; seeds matter at the second decimal (runs A, B, C: 0.424, 0.442, 0.443). *Biological interpretation:* none. *Limitations:* a rare event (5-year) cannot be quantified by a single run: pooled counts with Clopper–Pearson intervals are required.

![Figure 5.16. Which parameter uncertainty drives event-time uncertainty](../results/figures/fig_p3_11_parameter_vs_event_time_uncertainty.png)

**Figure 5.16.** Parameter uncertainty versus event-time uncertainty. *Setting:* prior P1; Sobol indices with $N=1024$ (10,240 runs). *Axes/panels:* left, box plots of event times within strata of $\beta/\beta_c$ (0–0.15, $n=275$; 0.15–0.25, 970; 0.25–0.35, 1044; 0.35–0.5, 1189; 0.5–0.75, 989); centre and right, $S_1$ (green) and $S_T$ (purple) for $\min(T,240)$ and for the event indicator. *What it shows:* the median event time rises with $\beta/\beta_c$ (about 124 to 210 months) but each stratum still spans about 100 months; the number of events falls from 252/275 in the first stratum to 39/989 in the last. $S_T$ for restricted time: $\beta$ 0.63, $q$ 0.34, $\rho$ 0.18, $N_0$ 0.14, $a$ 0.05; for the event: $\beta$ 0.72, $q$ 0.39, $\rho$ 0.28, $N_0$ 0.16, $a$ 0.10; $\eta$, $\delta_D$, $f_D\approx0$. *Mathematical interpretation:* under P1 $r$ is fixed, so $\beta$ replaces $r$ as the dominant source: **the ranking depends on which inputs are given distributions**. *Biological interpretation:* even at fixed immune strength, event timing is broadly spread. *Limitations:* depends on the P1 spreads.

---

### 5.6 Equilibrium-based predictions versus simulated events

For each P1 sample ($n=5000$) equilibria were recomputed from the exact cubic (never from the baseline $k_\text{det}$ formula) and used as predictors.

**Table 5.5.** Predictors versus outcomes (source: `results/tables/p3_mc_equilibrium_predictions_vs_events.csv`).

| predictor | vs event by 240 mo: accuracy (sens/spec) | vs long-horizon event: accuracy (sens/spec) |
|---|---|---|
| TF unstable ($\beta<\beta_c$; the $k_c$-type rule) | 0.465 (1.00/0.07) | 0.839 (1.00/0.20) |
| highest stable $N^*\ge\theta$ | 0.609 (0.99/0.33) | 0.983 (0.99/0.94) |
| lowest stable $N^*\ge\theta$ | 0.618 (0.99/0.35) | 0.992 (0.99/0.99) |
| attractor reached has $N_\infty\ge\theta$ | 0.620 (0.99/0.35) | 0.994 (0.99/1.00) |

Reasoning: **$k_c$ is a local stability boundary, not a relapse classifier.** It has perfect sensitivity but specificity 0.20 against the long-horizon event: 804 of 5000 samples (16%) have an unstable tumour-free state and still never reach $N=0.1$, because the stable nonzero attractor lies below the threshold (a low-burden "dormant" regime). Against the 240-month event its accuracy is 0.465 (1,872 samples have an attractor above $\theta$ but no event within 20 years: finite-horizon censoring). The attractor-based predictor agrees with the long-horizon event in 99.4% of samples; its 29 false negatives are overshooting or transient trajectories that end below $\theta$; there are 0 false positives. No tumour-free-stable sample produced an event under P1.

**The apparent detection boundary is sample specific.** For 300 P1 samples, $\beta$ was scanned over $[0.05,20]\times\beta$ (80 values) with a long-horizon run at each. At $\theta=0.1$, $\beta_\text{det}/\beta_c$ has median 0.72, IQR 0.58–0.86, range 0.07–1.04; the equilibrium-based boundary equals the simulated one in the median (ratio 1.000). The Phase 2 baseline value $k_\text{det}/k_c=0.7548$ is therefore not a universal number (verified at baseline for all seven thresholds; event at $0.98k_\text{det}$, none at $1.02k_\text{det}$). Values slightly above 1 (1.035) are grid resolution.

![Figure 5.17. Event probability by equilibrium regime and by attractor reached](../results/figures/fig_p3_12_event_probability_by_regime.png)

**Figure 5.17.** Event probability by regime and by attractor. *Setting:* prior P1_base, $n=5000$, $N_\text{detect}=0.1$. *Axes/panels:* left, by equilibrium regime of the sampled parameter vector (TF only $n=192$, 3.8%; TF+NZ $n=11$, 0.2%; NZ only $n=4718$, 94.4%; NZ+NZ $n=79$, 1.6%); right, by the attractor actually reached (NZ $n=4718$, NZ-high $n=3$, NZ-low $n=76$, TF $n=203$). Blue = event by 240 months, orange = event by the long horizon ($\ge10^5$ months). *What it shows:* TF-only and TF+NZ have 0.000 for both; NZ only has 0.446 (240 mo) and 0.840 (long); NZ+NZ 0.177 and 0.367; by attractor, NZ-high 1.000 (n = 3, tiny), NZ-low 0.145 and 0.342, TF 0.000. *Mathematical interpretation:* stability regime alone does not give the event probability; the gap between blue and orange bars is finite-horizon censoring (for NZ only, 0.446 → 0.840). *Biological interpretation:* even where a persistent tumour state exists, most simulated relapses are beyond 20 years. *Limitations:* the TF+NZ and NZ-high groups are very small, so their bars are not estimates.

---

### 5.7 Detection-threshold robustness

Common random numbers: the same 5000 P1 samples (run A) and dense paths for thresholds 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5 (model-defined burdens, **not imaging limits**).

**Table 5.6.** Simulated probabilities versus threshold (source: `results/tables/p3_threshold_event_probabilities.csv`).

| $\theta$ | P(60 mo) | P(120 mo) | P(240 mo) | cond. median $T$ | transient events | long-horizon P |
|---|---|---|---|---|---|---|
| 0.01 | 0.040 | 0.352 | 0.689 | 118.6 | 0.2% | 0.949 |
| 0.02 | 0.013 | 0.267 | 0.631 | 129.2 | 0.5% | 0.930 |
| 0.05 | 0.002 | 0.155 | 0.521 | 142.5 | 1.5% | 0.870 |
| **0.1** | 0.000 | 0.085 | 0.424 | 153.2 | 1.4% | 0.799 |
| 0.2 | 0.000 | 0.043 | 0.337 | 166.3 | 0.3% | 0.731 |
| 0.3 | 0.000 | 0.024 | 0.284 | 175.0 | 0.1% | 0.696 |
| 0.5 | 0.000 | 0.008 | 0.214 | 185.6 | 0.0% | 0.615 |

(Zeros at $\theta\ge0.1$ mean "below about 1/5000", not exactly zero; pooled P1 gives 0.0005 at $\theta=0.1$.)

![Figure 5.18. Detection-threshold sensitivity](../results/figures/fig_p3_09_detection_threshold_sensitivity.png)

**Figure 5.18.** Detection-threshold sensitivity. *Setting:* prior P1, $n=5000$, seven thresholds with common random numbers. *Axes/panels:* top left, event probability by 60/120/240 months and long-horizon probability (dashed black) versus $N_\text{detect}$ (log axis; bands 95% Wilson); top right, restricted mean event-free time (240 mo) and conditional median crossing time with IQR; bottom left, fraction of events that are transient/persistent/end below $\theta$, and accuracy of two equilibrium predictors; bottom right, per-sample boundary $\beta_\text{det}/\beta_c$ (median and IQR) for the simulated long run, the two equilibrium rules and the Phase 2 baseline formula $k_\text{det}/k_c$. *What it shows:* $P(240)$ falls from 0.689 to 0.214 across the thresholds, $P(120)$ from 0.352 to 0.008; the conditional median rises from 118.6 to 185.6 months; transient crossings stay below about 1.5% of events; the attractor-based rule has accuracy 0.994–1.000 against the long-horizon event whereas the "TF unstable" accuracy against the 240-month event falls from 0.73 to 0.26; simulated $\beta_\text{det}/\beta_c$ falls from about 0.95 to 0.50, close to the baseline formula (0.957 → 0.559) but with a wide IQR (0.32–0.67 at $\theta=0.5$). *Mathematical interpretation:* raising $\theta$ pushes the crossing later along an exponential growth, and moves the no-relapse boundary down toward the stability boundary's complement. *Biological interpretation:* what is called "relapse" depends on the assumed detectability of disease. *Limitations:* thresholds are arbitrary model values, not imaging limits; the formula curve is valid only at baseline-like parameters.

The threshold assumption changes $P(20\text{ y})$ by a factor of about 3, $P(10\text{ y})$ by a factor of about 45, and the conditional median by 67 months: it is as influential as the prior assumptions.

---

### 5.8 Separating the sources of uncertainty

**Table 5.7.** Sources of uncertainty for the 20-year probability (source: `results/tables/p3_mc_assumption_vs_sampling_uncertainty.csv`, `p3_threshold_event_probabilities.csv`).

| source | what it is | magnitude for $P(240)$ |
|---|---|---|
| Monte Carlo sampling | finite $n$ for a given model and prior | 95% half-width 0.014 at $n=5000$; 0.030 at $n=1000$ |
| Prior assumption | seven hypothetical distributions P1–P7 | range 0.317–0.468; spread 0.151 = 5.5× median MC interval width (0.027) |
| Threshold assumption | $N_\text{detect}=0.01\ldots0.5$ | 0.689 … 0.214 |
| Numerical | solver tolerances, event bracketing | 0/600 event flips; $\lvert\Delta T\rvert\le0.007$ months |

At 120 months the prior spread (0.198) is 12× the MC width; at 60 months the spread (0.036) is 19×. Correlation sensitivity (P4 vs P5, ±0.5 on $(q,a)$ and $(\beta,\rho)$) gives 0.468 vs 0.419, a 0.049 effect larger than the MC width. Solver tolerances (rtol $10^{-6}\ldots10^{-12}$, Radau, LSODA; 600 samples) produced no event flips; same seed gives identical results; parallel equals serial. **Priors, not the number of simulations, dominate the uncertainty**, which is why these probabilities are not predictions.

---

### 5.9 Conclusions of Phase 3 and what they implied for Phase 4

1. **Which parameters matter.** In the hypothetical ranges, net proliferation $r$ and immune killing $\beta$ dominate relapse occurrence and timing (restricted time $S_T$ 0.78 and 0.33; event 0.78 and 0.43), followed by $q$, $\rho$, $N_0$ and the threshold. $a$, $\delta_D$, $\eta$ and $f_D$ are near non-influential. Interactions are strong, though only $r\times\beta$ is resolved at second order. Four independent methods (Spearman, PRCC, Morris, Sobol) agree on the top three.
2. **Fragility near the horizon.** At the baseline reference a 3.3% decrease in $r$ removes the 240-month event; the binary outcome is not differentiable there.
3. **Basins.** A locally stable tumour-free state does not protect against all disturbances, but in the baseline-adjacent case the escape basin lies above plausible residual burdens; with strong immune stimulation and immune depletion, sub-threshold starts can reach a high-burden attractor. Event outcome and attractor are different things.
4. **Simulated probabilities** (P1, pooled): 0.00053 (5 y), 0.0910 (10 y), 0.4363 (20 y), conditional median 153.5 months; long-horizon probability 0.799. These are model-defined and hypothetical.
5. **$k_c$ is not a relapse classifier.** 16% of P1 samples have an unstable tumour-free state and never reach the threshold; the attractor-based rule reaches 99.4% agreement.
6. **Assumptions dominate sampling error** (priors 5.5×, threshold factor 3 for 20-year probability).

**Verdict used to decide Phase 4.** Phase 4 could begin as identifiability and model-structure work (synthetic recovery, endpoint/time-origin mapping, nuisance priors on the threshold and initial burden) but **should not begin as calibration of these simulated probabilities to patient data**, because they are prior-dependent and threshold-dependent. The concrete recommendations carried forward were: decide the clinical event definition and time origin first; test identifiability of $(r,\beta,q,\rho,N_0)$ from time-to-event curves with the threshold unknown; treat threshold, initial-burden distribution and fixed parameters as explicit nuisance assumptions; consider a stochastic birth–death–switching extension; and deliberately stress bistability since P1 barely samples it.

**Limitations of Phase 3.** All ranges and priors are hypothetical and loosely anchored; first-order Sobol intervals are wide and near-zero inputs are not rank-resolved; PRCC/Spearman assume monotone responses (checked graphically for the top inputs only); detection uses a 0.05-month grid; basins are 2-D slices for three parameter sets; bistable regimes rest on $n=79+11$ samples under P1; the long-run classifier clips the right-hand side and assumes convergence to an equilibrium (no limit cycle: supported by Phase 2 scans, not proved); the Monte Carlo spread represents parameter and initial-condition uncertainty only, not demographic stochasticity; no treatment, constant parameters, binary dormancy, mass-action killing. Tests: 111 pass (71 from Phases 1–2 plus 40 new for Phase 3).
