## 6. Phase 4: Identifiability, Synthetic Validation and Model Reduction

*Everything in this section is SYNTHETIC: the data are generated from a hypothetical model with known ground truth and an assumed observation model. Nothing here is clinical evidence, a validated measurement model or a patient risk. The METABRIC descriptive audit and the final feasibility decision are in Section 7.*

**Question addressed.** Phases 2–3 characterised the P–D–e model and its simulated relapse-time distribution. Phase 4 asks the prerequisite question for any link to data: *if we could observe this system the way a study would, which parameters could actually be learned?* The route is (i) fix the decision rules before seeing results, (ii) define explicitly what is observed, (iii) test identifiability exactly (structural), (iv) test it under noise and truncation (practical), (v) test recovery on synthetic data with known truth, (vi) reduce the model where the data cannot support its parameters, (vii) quantify how nuisance quantities (time origin, detection threshold, initial burden) confound event-time curves, and (viii) check survival estimators on synthetic cohorts.

### 6.1 Guardrails and pre-registered criteria

The reference is the normalised P–D–e model of Phase 2 ($s=\mu$, tumour-free $e=1$, $K=1$, so burdens are fractions of $K$). Reduced models are special cases fitted in separate experiments; no earlier source file was modified (new code is additive). Six hypothetical ground truths were fixed in `configs/phase4.yaml` (sha256[:12] `cfb17cec5c80`, stamped in every Phase 4 table and figure):

**Table 6.1.** Ground-truth parameter sets (source: configs/phase4.yaml; common values $r=0.12$, $a=0.005$, $\delta_D=0.001$, $\rho=0.04$, $\eta=0.08$ except G6).

| id | role | $\beta$ | $N_0$ | $f_D$ | note |
|---|---|---|---|---|---|
| G1 | baseline | 0.04 | $1.01\times10^{-4}$ | 0.990 | NZ-only regime; first crossing of 0.1 near 220 months |
| G2 | fast | 0.02 | $1.01\times10^{-4}$ | 0.990 | weak killing, early crossing |
| G3 | slow | 0.0567 | $1.01\times10^{-4}$ | 0.990 | no crossing within 240 months |
| G4 | controlled | 0.15 | $1.0\times10^{-3}$ | 0.990 | $\beta>\beta_c$, residual decays |
| G5 | low dormant | 0.04 | $1.01\times10^{-4}$ | 0.10 | mostly proliferating start |
| G6 | stimulated | 0.06 | 0.0328 | 0.493 | two-nonzero regime ($q=0.00158$, $a=0.306$, $\rho=0.240$, $\eta=1.637$, $e_0=0.606$); Phase 3 transient fixture |

**Proof standard.** Full rank of an exact series identifiability matrix at random rational points proves *generic local* structural identifiability (finite-to-one), not global identifiability. Rank deficiency is only evidence until an explicit null vector certifies it. A full-rank numerical sensitivity matrix proves nothing structural. Practical identifiability is a statement about noisy, sparse, truncated data (profile likelihood; Raue et al. 2009).

**Decision criteria written before any result** (`decision_criteria` in the config):

- *Practically identifiable:* the 95% profile interval is finite on both sides within the bounds in at least 80% of replicates, and the median $|\ln$ recovery error$|\le\ln 1.5$ at the benchmark condition ($\sigma=0.10$, 6-month visits, 240-month follow-up truncated at detection).
- *A (descriptive comparison)*, *B (exploratory population calibration)* and *C (mechanistic calibration unsuitable)* have separate pre-registered conditions; they are answered in Section 7 using the evidence below.

Writing criteria first prevents thresholds being tuned to the outcome; one judgement (the reading of criterion B, §6.8) was made after seeing results and is flagged as such.

### 6.2 Observation model

Five levels are kept apart so that "what the model has" is never confused with "what a study sees":

1. hidden states $x(t)=(P,D,e)$;
2. total burden $N=P+D$ (a mathematical quantity);
3. hypothetical *measurement signals* $y_N=c\,N\,e^{\sigma\varepsilon}$ (likewise $y_P$, $y_D$, $y_e$), multiplicative log-normal noise, signal scale $c$ known (1) or unknown; values below a limit of detection (LOD) are *left-censored with a Tobit likelihood*, not dropped;
4. the threshold-defined recurrence event, seen only through visits as an *interval-censored* record $L<T\le R$;
5. an observed clinical outcome (time origin, detection delay, competing death; Section 7).

**Observation stops at detection.** Once $N$ reaches $N_\text{detect}$ the relapse is diagnosed and treated, so the untreated trajectory is observed only up to $T_\text{detect}$. "Long follow-up" therefore means observation up to detection (216 of 240 months for G1), not $N(240)$. Base assumptions: log-noise sd 0.10 (about 10% CV; 0.05–0.5 in sensitivity runs), LOD $10^{-6}$, visits every 6 months (1–24 in sensitivity runs), follow-up 240 months. To avoid an *inverse crime*, data are generated with Radau at rtol $10^{-12}$ and fitted with DOP853 at rtol $10^{-9}$ using analytic forward sensitivities.

### 6.3 Structural identifiability (exact)

**Method.** The Taylor coefficients of each output at $t=0$ are computed by a series recurrence (including division for $P/(h+P)$) with forward-mode automatic differentiation, modulo $p=2^{61}-1$ or in exact rationals, at 3 random rational points with generic parameters. Taylor series were verified against a $10^{-13}$ ODE solution; AD gradients against finite differences. The rank is tracked against Taylor order (flat over the last 6 of 30 orders at every point); null vectors are exact rationals.

**Result 1 – generic local identifiability.** With unknown sets A–E (up to 13 unknowns: $r,q,a,\delta_D,\beta,\rho,\eta$; plus $P_0,D_0$; plus $h,\mu$; plus $e_0$; plus signal scale $c$) and $K=1$, the rank is *full for all four observation sets* $N$, $(P,D)$, $(N,e)$, $(P,D,e)$. In principle every parameter is locally identifiable even from $N(t)$ alone.

**Result 2 – a scale gauge when $K$ and $c$ are both unknown (set F).** Deficiency is exactly 1 for every observation set. The exact rational null vector has relative components $(-0.868,+0.868,+0.868,+0.868,-0.868,+0.868)$ on $(\eta,P_0,D_0,h,c,K)$, i.e. the burden-scale gauge
$$x\to\lambda x,\quad K\to\lambda K,\quad h\to\lambda h,\quad \eta\to\eta/\lambda,\quad c\to c/\lambda.\tag{6.1}$$
This is the burden analogue of the immune-scale gauge removed in Phase 2 (M1). Fixing $K=1$ removes it.

**Result 3 – the low-burden sub-models are *not* identifiable.** In the linear regime $x'=Mx$ with $M=\begin{pmatrix} m& a\\ q&-(a+\delta_D)\end{pmatrix}$, $m=r-q-\beta e_0$, and $N$-only observation, the output has four invariants ($\lambda_\pm$ and two amplitudes) for six unknowns: symbolic rank 4 of 6 (null dimension 2). With $r$ and $\beta$ separate: rank 4 of 7 for $N$ alone, and rank 6 of 7 even for $(P,D)$, where the single deficiency is exactly that $r$ and $\beta$ enter only as $r-\beta$. The logistic term $rPN/K$ and the immune equation are what separate them, and both are second order in the burden.

**The decisive distinction.** Generic structural identifiability (full model, generic point) coexists with non-identifiability of the same parameters in the regime the data actually visit. The information separating $r$ from $\beta$ scales with tiny nonlinear terms; whether it is usable is a *practical* question (§6.4). Global identifiability was not proven.

**Result 4 – event-only observation.** One event time is one scalar. For G1: $dT/d\ln r=-563$, $dT/d\ln\beta=+218$, $dT/d\ln q=+146$, $dT/d\ln f_D=+186$, $dT/d\ln N_\text{detect}=+51$ months. By the implicit-function theorem every level set $\{T=\text{const}\}$ is a manifold of dimension $p-1=9$ in the ten inputs, so no parameter is identifiable from one event time. Numerically, changing $\beta$ by $\pm10$–$25\%$ and re-solving $r$ to keep $T=219.9$ months leaves $N(60)$ within 2–5% and $N(120)$ within 5–15% of the original curve. In the single-mode linear limit
$$T=\frac{1}{\lambda_+}\ln\frac{N_\text{detect}}{n_1u_0},\tag{6.2}$$
so scaling $N_0$ by $s$, scaling the threshold, and shifting the time origin by $\tau$ are *exactly confounded* at fixed $\lambda_+$: multiplying $N_0$ by 0.1, 0.3, 3, 10 shifts $T$ by $+52.37,+27.40,-25.06,-52.77$ months against the prediction $-\ln s/\lambda_+=+52.33,+27.36,-24.97,-52.33$ ($\lambda_+=0.044$ for G1). For a population, event data identify the push-forward law of $T$, a distribution of one scalar, not the parameter distribution.

### 6.4 Fisher spectra and profile likelihoods

**Fisher spectra.** Noise-scaled sensitivity matrices (analytic forward sensitivities, 9 log/logit parameters $r,q,a,\delta_D,\beta,\rho,\eta,N_0,f_D$; G1, 6-month visits, truncated at detection) are summarised by singular values; a singular value $s>3$ means a predicted sd below 0.33 in $\ln$. The number of such directions at $\sigma=0.10$ is **4 of 9 for $N$, 6 for $(P,D)$, 5 for $(N,e)$, 7 for $(P,D,e)$** (3, 6, 4, 6 at $\sigma=0.25$). Truncation at detection costs at most one direction because the nonlinear regime is reached only after detection.

![Figure 6.1. Fisher spectra and well-determined directions](../results/figures/fig_p4_04_identifiability_by_observed_states.png)

**Figure 6.1.** Fisher spectra and well-determined directions. *Setting:* G1 baseline, 6-month visits, LOD $10^{-7}$, linearised at the truth; SYNTHETIC. *Axes/panels:* left, singular values of the noise-scaled sensitivity matrix (log axis) against index 1–9 for the four observed-state sets, dashed line at $s=3$; centre, count of directions with $s>3$ for $\sigma_{\ln}=0.05,0.1,0.25$ (dotted line at 9 unknowns); right, the same count with and without truncation at detection. *What it shows:* the $N$-only spectrum falls steeply, from about $10^3$ at index 1 to about $10^{-5}$ at index 9, crossing $s=3$ after the fourth value, while the $(P,D,e)$ spectrum stays above 3 until index 7. The bars show 4/6/5/7 directions at $\sigma=0.05$ and 0.1, and one fewer for $N$ (3), $(N,e)$ (4) and $(P,D,e)$ (6) at 0.25. Truncation makes no visible difference (the two lines coincide). *Mathematical interpretation:* the gap between a full-rank structural matrix and a spectrum spanning eight orders of magnitude is the practical deficit; the weak directions are combinations the data barely constrain. *Biological interpretation:* total burden alone supports about four composite quantities (dominated by the growth-rate combination); observing proliferating and dormant cells and immune activity separately adds information. *Limitations:* linearisation at the truth, one ground truth, hypothetical noise.

**Information on the $r$–$\beta$ split does not grow with the initial burden.** Scaling $N_0$ from 0.1 to 100 times baseline leaves the estimator correlation at 0.9998, with Wald sd of $\ln r$ between $3\times10^2$ and $2.6\times10^3$ and of $\ln\beta$ between $5\times10^2$ and $4.6\times10^3$ (`p4_fisher_rbeta_vs_initial_burden.csv`).

![Figure 6.2. Estimator correlation heatmaps (Fisher, log-parameters)](../results/figures/fig_p4_02_parameter_correlation_heatmap.png)

**Figure 6.2.** Estimator correlation from the linearised (Gauss–Newton) covariance. *Setting:* G1, $\sigma_{\ln}=0.10$, three observation sets; $|$correlation$|>0.9$ annotated; SYNTHETIC. *Axes/panels:* each panel is a 9×9 matrix over $r,q,a,\delta_D,\beta,\rho,\eta,N_0,f_D$ (colour: $-1$ blue to $+1$ dark red); left $N$, centre $(N,e)$, right $(P,D,e)$. Panel titles list the predicted sd of each log-parameter (partly cut off in the figure; legible values include $r=2.2$, $q=6.2$, $a=13$, $\delta_D=1.3\times10^2$ for $(N,e)$ and $r=0.054$, $q=0.07$, $a=0.066$, $\beta=0.15$, $\rho=0.33$, $\eta=1.9$, $N_0=0.06$, $f_D=0.12$ for $(P,D,e)$). *What it shows:* for $N$ alone nearly all parameters except $N_0$ are locked at $|\text{corr}|\approx1.00$ (positive among $r$, $q$, $\delta_D$, $\beta$, $f_D$; opposed to $a$, $\rho$, $\eta$), with sds of order $10^2$ in $\ln$; $N_0$ is nearly uncorrelated. With $(N,e)$, strong correlations persist among $\{r,q,\rho,\eta\}$ and $\{a,\beta,f_D\}$. With $(P,D,e)$ most cells are pale; the remaining strong pairs are $r$–$\beta$ (0.90), $q$–$\delta_D$ (0.90) and $\rho$–$\eta$ (0.97). *Mathematical interpretation:* correlation $\pm1$ with huge sd marks a near-null direction of the Fisher matrix, a ridge in log-parameter space (the $r-\beta$ type redundancy of §6.3). *Biological interpretation:* growth, dormancy entry/exit and immune killing trade off against each other in burden data; separating them needs compartment- or immune-level observations. *Limitations:* local linearisation; correlations far from the truth may differ.

**Profile likelihoods.** For each parameter the likelihood is maximised over the others on a grid of fixed values; the 95% interval is where $\Delta\chi^2<3.84$, and a parameter is practically identifiable if the interval is finite on both sides within $\pm\ln 10$ (grid step 0.33 in $\ln$, so widths of tight parameters are approximate).

![Figure 6.3. Profile likelihoods for N-only and (P, D, e) data](../results/figures/fig_p4_03_profile_likelihoods.png)

**Figure 6.3.** Profile likelihoods. *Setting:* G1, $\sigma_{\ln}=0.10$, one synthetic data set per row; SYNTHETIC. *Axes/panels:* nine columns (one per parameter), top row observed $N$, bottom row observed $(P,D,e)$; $x=\ln(\theta/\hat\theta)$ over the profile window, $y=\Delta\chi^2$ (clipped at 12); dashed line at 3.84; "(open)" marks an interval not closed on at least one side. *What it shows:* with $N$ only, just $N_0$ has a closed V-shaped profile; $r$ rises steeply on the low side but stays flat above, $a$ rises only on the high side, and $q,\delta_D,\beta,\rho,\eta,f_D$ are essentially flat ($\Delta\chi^2\lesssim2$ over the whole window). With $(P,D,e)$, $r,q,a,\beta,\rho,N_0,f_D$ are sharply closed, whereas $\delta_D$ is bounded only on one side and $\eta$ is flat (open). *Mathematical interpretation:* a flat profile means a likelihood ridge: parameter changes are compensated by others without degrading the fit. *Biological interpretation:* immune suppression $\eta$ and dormant-cell loss $\delta_D$ leave no measurable trace even with all three states observed. *Limitations:* a single data set per row; coarse profile grid; the surface is ridge-like rather than multi-basin (7/8 starts for $N$ and $(P,D)$, 6/8 for $(P,D,e)$ reached fits within $\Delta\chi^2<3.84$ of the best).

**Table 6.2.** Profile-likelihood identifiability, G1, $\sigma=0.10$ (source: results/tables/p4_profile_likelihood_summary.csv; interval width factor in parentheses).

| parameter | $N$ | $P,D$ | $N,e$ | $P,D,e$ |
|---|---|---|---|---|
| $r$ | one-sided | finite (×14) | finite (×11) | **finite (×1.08)** |
| $q$ | not identifiable | **finite (×1.07)** | one-sided | **finite (×1.07)** |
| $a$ | one-sided | **finite (×1.13)** | one-sided | **finite (×1.11)** |
| $\delta_D$ | not identifiable | one-sided | not identifiable | one-sided |
| $\beta$ | not identifiable | not identifiable | one-sided | **finite (×1.59)** |
| $\rho$ | not identifiable | not identifiable | finite (×14) | **finite (×1.74)** |
| $\eta$ | not identifiable | not identifiable | not identifiable | not identifiable |
| $N_0$ | finite (×1.15) | finite | finite | finite |
| $f_D$ | not identifiable | finite | not identifiable | finite |

![Figure 6.4. Structural versus practical identifiability](../results/figures/fig_p4_11_structural_vs_practical_identifiability.png)

**Figure 6.4.** Structural versus practical identifiability. *Setting:* left, exact series analysis; right, the profile results of Table 6.2; SYNTHETIC. *Axes/panels:* columns are observed-state sets $N$, $(P,D)$, $(N,e)$, $(P,D,e)$; left rows are unknown sets A (7 parameters) to F ($E+K$, burden scale unknown); right rows are the nine parameters; green = identifiable, yellow = one-sided, red = not identifiable / deficient. *What it shows:* the left panel is green in every cell for A–E and red ("deficiency 1") for F; the right panel is largely red/yellow for $N$ (only $N_0$ green), improves with $(P,D)$ and $(N,e)$, and is mostly green for $(P,D,e)$ except $\delta_D$ (yellow) and $\eta$ (red). *Mathematical interpretation:* structural full rank is necessary but not sufficient; the right panel is the one that matters for real data. *Biological interpretation:* a model can be formally identifiable and still unlearnable from any realistic dataset. *Limitations:* local generic structural analysis; the practical panel comes from one noise realisation per column.

### 6.5 Synthetic parameter recovery

**Design.** One-factor-at-a-time variation around a base condition (G1, observed $N$, $\sigma=0.10$, 6-month visits, 240-month limit truncated at detection, $N_0$ and $f_D$ unknown): 23 conditions × 8 noise replicates = 184 constrained fits, each from 8 random starts spread over bounds of truth ÷30…×30 (log scale); estimates on a bound are flagged. Fitted: $r,q,a,\delta_D,\beta,\rho,\eta,N_0,f_D$ ($\mu,h,K$ fixed, $e_0$ known).

**Base condition.** Median $|\ln$ error$|$ per parameter: $r$ 0.77, $q$ 0.64, $a$ 0.75, $\delta_D$ 1.80, $\beta$ 1.78, $\rho$ 1.43, $\eta$ 2.19, $f_D$ 1.62, $N_0$ 0.03; only $N_0$ is within a factor 1.5 in all replicates. 62.5% of fits end on a bound. Yet the *crossing time is predicted within 10% in 87.5%* of replicates (median error 1.1%).

![Figure 6.5. Recovery: true versus estimated parameters](../results/figures/fig_p4_01_parameter_recovery_true_vs_estimated.png)

**Figure 6.5.** Parameter recovery, true versus estimated. *Setting:* six ground truths × four observed-state sets, each point one noisy synthetic data set; $\sigma_{\ln}=0.10$, 6-month visits, truncated at detection, 8-start multistart, box bounds truth ÷30…×30; SYNTHETIC. *Axes/panels:* nine log–log scatter panels (one per parameter), $x$ = true value, $y$ = estimate, black line = identity; colours: observed $N$ (blue), $(P,D)$ (green), $(N,e)$ (orange), $(P,D,e)$ (purple). *What it shows:* points for $N_0$ lie on the identity line across about three decades ($10^{-4}$ to $\approx0.03$); for $r$, $q$, $a$, $\delta_D$, $\beta$, $\rho$, $\eta$, $f_D$ the estimates form vertical streaks spanning one to three decades above a single true value, sometimes reaching the edges of the permitted box. Purple $(P,D,e)$ points cluster much closer to the truth for $r$, $q$, $a$, $\beta$, $\rho$ but stay spread for $\delta_D$ and $\eta$. *Mathematical interpretation:* the vertical streaks are the ridges of §6.4 seen in the estimates; ill-determined parameters wander along them. *Biological interpretation:* fitting total burden would give widely different mechanistic stories for the same trajectory. *Limitations:* only six truths with 8 replicates; bounds are hypothetical, and 62% of base fits touch a bound.

**Table 6.3.** Recovery summary by condition (source: results/tables/p4_recovery_summary.csv; $T$ error = relative error of the predicted first-crossing time; n.d. = no crossing).

| condition (observed $N$ unless stated) | directions $s>3$ | median $|\ln$ err$|$ | bound-hit fraction | median $T$ error |
|---|---|---|---|---|
| base $\sigma=0.10$ | 4 | 1.43 | 0.63 | 1.1% |
| $\sigma=0.001/0.01/0.05/0.25/0.5$ | 6/5/4/3/2 | 0.35/0.61/1.12/2.21/2.11 | 0/0/0.13/1.0/0.88 | 0.0%/0.2%/0.6%/1.6%/4.6% |
| schedule 1/3/12/24 months | 4 | 1.28/1.78/1.46/1.68 | 0.50/0.38/0.75/0.63 | 0.6%/0.7%/1.2%/1.5% |
| follow-up 60/120 months | 2/3 | 2.91/1.70 | 1.0/1.0 | **43.5%/26.8%** |
| observed $(P,D)$ | 6 | 0.36 | 0.63 | 0.4% |
| observed $(N,e)$ | 5 | 1.69 | 1.0 | 1.1% |
| observed $(P,D,e)$ | 7 | 0.07 | 0.75 | 0.3% |
| initial state known | 3 | 1.33 | 0.38 | 0.7% |
| scale $c$ unknown (true 2) | 4 | 1.42 | 0.38 | 25% |
| high LOD $3\times10^{-4}$ | 3 | 1.83 | 0.50 | 1.1% |
| G2 / G3 / G4 / G5 | 3/4/2/3 | 1.17/1.68/1.87/1.78 | 0.88/0.88/0.88/0.50 | 2.2%/27.5%/n.d./2.0% |
| **G6 (nonlinear transient)** | **5** | **0.09** | 1.0 | n.d. |

![Figure 6.6. Recovery and prediction error versus measurement noise](../results/figures/fig_p4_05_recovery_error_vs_noise.png)

**Figure 6.6.** Error versus noise. *Setting:* G1, observed $N$, 6-month visits, truncated at detection; SYNTHETIC. *Axes/panels:* left, median $|\ln$ recovery error$|$ of each parameter against noise sd (natural-log units, $10^{-3}$ to 0.5), dotted line at ×1.5; right, median crossing-time relative error and median $|\ln N(240)|$ error against noise. *What it shows:* only $N_0$ falls clearly with decreasing noise (from about 0.2 at $\sigma=0.5$ to about $4\times10^{-4}$ at 0.001). The other parameters stay between roughly 0.25 and 3.5, above the ×1.5 line almost everywhere; at $\sigma=0.001$ the median error over parameters is still 0.35, and the figure footnote attributes the floor to non-convergence and ridges rather than data quality. Prediction error is monotone: crossing-time error rises from about $2\times10^{-4}$ at $\sigma=0.001$ to about 0.05 at 0.5, and $|\ln N(240)|$ error from about 0.007 to about 0.4. *Mathematical interpretation:* removing noise cannot remove a structural/practical ridge, only shrink the $\sigma$-dependent part of the error. *Biological interpretation:* better assays alone would not identify the mechanisms from burden. *Limitations:* 8 replicates per point; error floors depend on the optimiser.

![Figure 6.7. Recovery versus visit frequency and follow-up](../results/figures/fig_p4_06_recovery_error_vs_observation_frequency.png)

**Figure 6.7.** Sampling frequency and follow-up. *Setting:* G1, observed $N$, $\sigma_{\ln}=0.10$; SYNTHETIC. *Axes/panels:* left, median $|\ln$ recovery error$|$ against the interval between visits (1–24 months); centre, the same against the follow-up limit (60–240 months; observation also stops at detection); right, number of observations (black, left axis) and median crossing-time relative error (orange, right axis) against the visit interval. *What it shows:* in the left panel the curves wander without trend ($\eta$ and $\delta_D$ stay between about 1.3 and 3.5; $N_0$ stays below 0.1), so denser sampling does not repair the weak directions. In the centre, errors decline modestly as follow-up grows from 60 to 240 months, but $\eta$, $\delta_D$, $\beta$, $f_D$ remain above 1. In the right panel observations fall from about 220 (1 month) to about 10 (24 months) while the prediction error only rises from about 0.006 to 0.015. *Mathematical interpretation:* the missing information is in the dynamics visited, not the sampling density. *Biological interpretation:* more frequent scans do not teach the mechanism; observing the right compartments or the nonlinear phase would. *Limitations:* short follow-up still wrecks extrapolation (Table 6.3: 43.5% $T$ error at 60 months).

![Figure 6.8. Observed versus reconstructed trajectories](../results/figures/fig_p4_07_observed_vs_reconstructed_trajectories.png)

**Figure 6.8.** Observed versus reconstructed trajectories. *Setting:* synthetic data with $\sigma_{\ln}=0.10$ and 8-start fits; red dashed line = detection threshold; SYNTHETIC. *Axes/panels:* burden (fraction of $K$) or signal on a log axis against time; left G1 observed $N$; centre G1 observed $P$, $D$ and $e$; right G6 observed $N$ over the full 240-month window. Black = hidden true $N(t)$, markers = data, orange = near-best fits ($\chi^2\le$ best + 3.84) extrapolated. *What it shows:* in the left panel the data hug the true curve and observation stops at about 216 months at detection; the near-best fits agree in-window but fan out upward beyond the stop. In the centre, $P$ starts near $5\times10^{-6}$ and rises, $D$ stays near $10^{-4}$ before rising, and $e$ (green, about 1–1.5) is observed throughout. In the right panel G6 rises to about 0.06, drops to about $2\times10^{-4}$ near month 70, then forms a second hump near $5\times10^{-3}$. *Mathematical interpretation:* a good in-window fit does not guarantee correct extrapolation or correct parameters; G6, which enters the nonlinear regime, constrains $r,\beta,\rho,\eta$. *Biological interpretation:* a non-monotone transient carries mechanistic information that smooth monotone growth does not. *Limitations:* single realisations; the fits shown are a small sample of the ridge.

**Findings from Table 6.3 and Figures 6.5–6.8.**

1. *Noise.* Even at $\sigma=0.001$ the median error is 0.35: six directions are determined and three are not. Raising noise steadily degrades the determined directions (6 → 2).
2. *Sampling frequency barely matters* (Figure 6.7).
3. *Truncation and short windows hurt prediction:* 60-month follow-up gives 43.5% $T$ error, 120 months gives 26.8%.
4. *Regime matters.* In the controlled regime G4 (decay only) two directions are determined. For G6, which enters the immune-feedback regime, $r,\beta,\rho,\eta$ are recovered from $N$ alone ($|\ln$ errors$|$ 0.07, 0.09, 0.05, 0.09) but $q$ and $\delta_D$ are not (3.4, the bound): which parameters are learnable depends on which regime the data visit.
5. *Observed states.* Adding $P$ and $D$ recovers $r,q,a,f_D$ (errors 0.02–0.36); adding $e$ recovers $\beta$ and $\rho$ (0.07, 0.11). $\delta_D$ (1.41) and $\eta$ (2.10) stay unidentified.
6. *Prediction is not recovery.* Predicting $T$ to 1–2% does not support any mechanistic parameter.
7. *Wald diagnostics* cover the well-determined directions only 8/32 times because estimates sit far from the linearisation point; profile likelihoods are the reliable tool.

### 6.6 Redundancy and model reduction

Candidate reductions were derived from the exact null space and Fisher directions and tested on data generated by the *full* model (6 replicates per cell; observed $N$ and $(P,D,e)$; truths G1, G2, G5, G6): **M7** immune frozen ($\rho=\eta=0$), **M6** M7 with $\delta_D=0$, **M5** M6 with $f_D$ fixed at the equilibrium ratio $q/(q+a)$. M9 denotes the full model.

**Accept/reject on three criteria:** (i) AIC, (ii) parameters better determined, (iii) prediction retained.

- *Burden-only, sub-threshold data (G1, G2, G5):* M7, M6, M5 are preferred by AIC (median $\Delta$AIC $-3.5$ to $-6.4$; not worse than full + 2 in 100% of replicates) and predict $T$ equally well (median error 0.3–3.0% vs 0.3–1.4%): criteria (i) and (iii) pass. Criterion (ii) *fails*: individual parameters of the reduced models stay badly determined (median max Wald sd of $\ln$ parameters 2–484, above 30 in 10 of 24 cells). The reduction removes untestable parameters, not the redundancy between $r$, $\beta$, $q$.
- *$(P,D,e)$ data and the nonlinear G6:* reductions are rejected decisively ($\Delta$AIC +29 to +30,000), because $\rho,\eta$ matter.

![Figure 6.9. Model reduction: AIC, prediction and precision](../results/figures/fig_p4_13_model_reduction.png)

**Figure 6.9.** Nested model reduction. *Setting:* data generated from the full model; M9 → M7 → M6 → M5 (full → immune frozen → $\delta_D=0$ → $f_D$ at equilibrium ratio); lines per truth/observation pair (G1, G2, G5, G6; solid = observed $N$, dashed = observed $(P,D,e)$); SYNTHETIC. *Axes/panels:* left, AIC relative to the full model (lower = reduced model preferred); centre, crossing-time relative error (log axis); right, Wald sd of $\ln(r-\beta)$ (log axis). *What it shows:* on the left most lines stay near 0 (reduction harmless or preferred), but G6 with $N$ rises to about $8\times10^3$–$9\times10^3$ and G6 with $(P,D,e)$ to about $1.9\times10^4$ at M7 and $3\times10^4$ at M5, i.e. reduction is rejected there. In the centre, crossing-time errors are mostly 0.3–3% (G5 with $(P,D,e)$ reaches about 0.1 at M5). On the right, the sd of $\ln(r-\beta)$ is of order 10–100 for several $N$-only fits at M9, falling to about 1 at M5 for some truths, and is $10^{-2}$–$10^{-1}$ for $(P,D,e)$; G6 curves have missing points where no estimate was produced. *Mathematical interpretation:* $r-\beta$ (more generally $\lambda_+$) is the well-determined combination even when $r$ and $\beta$ separately are not. *Biological interpretation:* the data constrain net growth after immune control, not the separate strengths of proliferation and killing. *Limitations:* 6 replicates per cell; truths limited to four.

**Which combination is recovered?** From $N$ alone the dominant low-burden eigenvalue $\lambda_+$ is recovered to about 3% (median $|\ln$ error$|$ 0.032 G1, 0.030 G2, 0.028 G5, 0.004 G6) together with $N_0$ (3–6%); $m=r-q-\beta$ to 5–9% (22% G6); whereas $g=r-\beta$ (38%, 30%, 42%), $\lambda_-$ (errors 0.42, 0.49, 3.4) and the individual parameters are poor. With $(P,D,e)$: $\lambda_+$ to 0.7–1.2%, $\lambda_-$ to 3–12%, $m$ to 1–7%.

**Proposed reduced parameterisation** (justified only for the sub-threshold, burden-only regime): report $(\lambda_+,N_0)$, plus $\lambda_-$ only with $(P,D)$ data, and treat $r,\beta,q,a,\delta_D,\rho,\eta,f_D$ as nuisance parameters. The full model remains the reference implementation.

### 6.7 Time-origin, threshold and initial-burden confounding

A clinical time is hypothetically $t_\text{clin}=\tau+t_\text{model}+d$ (recurrence) or $\tau+t_\text{model}+L$ (death first), with $\tau$ the unknown offset between the clinical origin and "post-treatment residual disease", $d$ a detection delay and $L$ a death lag. These are *nuisance assumptions, not estimates*. The population is the Phase 3 prior P1 run A (n = 5000, $N_\text{detect}=0.1$).

**Detection threshold.** Moving $N_\text{detect}$ over 0.01–0.5 changes $F(240)$ from 0.689 to 0.214 and the time to 10% cumulative probability from 76 to 182 months (no sample starts above any threshold).

![Figure 6.10. Detection-threshold effect on the event-time distribution](../results/figures/fig_p4_09_detection_threshold_uncertainty.png)

**Figure 6.10.** Detection-threshold uncertainty. *Setting:* prior P1, n = 5000, thresholds 0.01–0.5 (model-defined burdens, not imaging limits); SYNTHETIC. *Axes/panels:* left, cumulative first-crossing probability against model time, one curve per threshold (colour bar, log scale from $10^{-2}$); centre, time to reach 10%, 20%, 30% cumulative probability against $N_\text{detect}$; right, $F$ at 60, 120, 240 months against $N_\text{detect}$. *What it shows:* the curves fan out from the lowest threshold (top, reaching about 0.69 at 240 months) to the highest (bottom, about 0.21). In the centre panel the 10% time rises from about 76 months at 0.01 to about 182 at 0.5, and the 30% level is reached within 240 months only for thresholds below about 0.26. In the right panel $F(240)$ falls monotonically from about 0.69 to about 0.21, $F(120)$ from about 0.35 to about 0.01, and $F(60)$ from about 0.04 to 0. *Mathematical interpretation:* the threshold acts both as a horizontal shift (crossing a larger burden takes longer) and as a vertical loss (fewer cross before 240 months). *Biological interpretation:* what counts as "detectable relapse" changes the apparent late-relapse fraction by a factor of about three. *Limitations:* the threshold is a hypothetical modelling choice with no clinical anchor.

**Time origin and delay.** The shift is exact and additive: $\tau$ and $d$ enter only through $\tau+d$. The time to 10% cumulative probability moves from 125 to 128/131/137/143/149 months for $\tau+d=3/6/12/18/24$, while $F(120)$ falls from 0.085 to 0.077/0.066/0.051/0.038/0.025. Scaling the initial burden ×3 and ×10 shifts the median by $-18.76$ and $-39.51$ months against the linear-regime prediction $-\ln s/\lambda_+$ of $-18.66$ and $-39.11$ (and $+36.68$ vs $+36.62$ at ×0.1). A death-first share of 0.105 or 0.272 (lag mean 24 months) moves the time to 20% from 150 to 153/155 months and $F(240)$ from 0.424 to 0.420/0.415: a small effect.

![Figure 6.11. Time-origin offset, burden scale and equivalence map](../results/figures/fig_p4_10_time_origin_uncertainty.png)

**Figure 6.11.** Time-origin uncertainty and the equivalence map. *Setting:* prior P1 event times, hypothetical mappings; SYNTHETIC. *Axes/panels:* top row: (i) cumulative probability against clinical time for $\tau=0,3,6,12,18,24$ and $\tau+d=36$; (ii) horizontal shift of the curve at the 5/10/20/30% levels against $\tau$; (iii) curves for initial burden ×0.1, ×0.33, ×1, ×3, ×10. Bottom row: (iv) death-first mapping with shares 0, 0.105, 0.272; (v) the $\tau=12$ reference curve (black, with a light-blue band) against the best threshold-alone (0.13572) and best $N_0$-scale-alone (×0.5) reproductions; (vi) minimum sup-norm distance between curves for $\tau=3,6,12,24$ with DKW noise lines for n = 2380 (about 0.028) and n = 5000 (about 0.019). *What it shows:* in (i) the curves are translated copies; (ii) is a straight line of slope 1; (iii) shifts the curves left and upward for larger $N_0$ (about 0.51 at ×10 versus about 0.34 at ×0.1 at 240 months); (iv) the curves nearly overlap; (v) the $\tau=12$ curve is almost indistinguishable from the two alternatives; in (vi) the distances for $\tau=12$ are about 0.009 ($N_0$ scale) and 0.022 (threshold), both below the n = 2380 noise level, whereas $\tau=24$ by threshold alone (about 0.046) exceeds it and by $N_0$ alone (about 0.025) does not. *Mathematical interpretation:* the confounding is exact in the linear regime (§6.3) and only weakly broken by nonlinearity. *Biological interpretation:* a 12-month ambiguity in when disease "starts" is equivalent to being off by a factor of 2 in the residual burden. *Limitations:* one population, one death-lag assumption; the map shows what cannot be separated, not what $\tau$ is.

**Equivalence-map summary** (sup-norm distance over $[0,240]$ after the best compensating change, versus the Dvoretzky–Kiefer–Wolfowitz noise $\varepsilon=0.0278$ for n = 2380 and 0.0192 for n = 5000): $\tau=12$ months is reproduced to 0.0092 by an $N_0$ rescale ×0.5 (below $\varepsilon$ for both n) and to 0.0224 by a threshold of 0.136 (below $\varepsilon$ for n = 2380); $\tau=24$ by $N_0\times0.2$ to 0.0250 (below $\varepsilon$ for n = 2380). A factor-3 change in $N_0$ equals a shift of 19 months ($N_0/3$: 0.0134, below $\varepsilon$) or $-28$ months ($N_0\times3$: 0.0450, above); a factor-3 threshold change equals +49 months (0.0490, above). **At METABRIC sample size, an offset up to about 12 months cannot be separated from a factor-2 change in the initial burden using event times alone.**

![Figure 6.12. Clinical endpoint mapping: what is documented and what is a nuisance](../results/figures/fig_p4_12_clinical_endpoint_mapping_limitations.png)

**Figure 6.12.** Clinical endpoint mapping limitations (the one CLINICAL, descriptive figure of Phase 4; interpreted fully in Section 7). *Setting:* a schematic of the order of events plus a descriptive Kaplan–Meier summary of the cBioPortal METABRIC export (RFS composite; non-disease deaths censored; n = 2,380, 964 events); no model fitted or overlaid. *Axes/panels:* left, a schematic (not to scale) running from initial diagnosis (origin of OS: documented) through surgery, end of treatment (not recorded), model $t=0$, model threshold crossing and recorded recurrence to disease-specific death, with arrows marking the unknown windows $\tau$, $d$, $L$; right, cumulative RFS event probability against months from the undocumented RFS origin at 60, 120, 180, 240 months with confidence bars and numbers at risk (1581, 984, 521, 162). *What it shows:* the cumulative RFS event probability rises from roughly 0.26 to 0.50 while the at-risk count shrinks to 162 at 240 months. *Mathematical interpretation:* $\tau$, $d$ and $L$ shift or mix the observed curve exactly as in Figure 6.11, and the data supply only event times. *Biological interpretation:* the clinical event is a composite whose origin and delay are unrecorded. *Limitations:* descriptive only, no calibration; details in Section 7.

### 6.8 Synthetic survival validation

**Truth.** The pooled prior-P1 population (3 independent runs, 15,000 individuals; $F(60/120/240)=0.0005/0.0910/0.4363$). 200 resampled cohorts of n = 2380 (and n = 1000) were drawn.

**Right censoring.** Kaplan–Meier under administrative censoring at 60/120/240 months, with or without independent loss to follow-up, is unbiased (|bias| ≤ 0.0007, about one Monte Carlo standard error), with 95% coverage 0.94–0.97 at 120 and 240 months. At 60 months coverage is 0.65–0.69 because the truth (0.0005) yields about one event per cohort; that interval is degenerate, not meaningful.

**Interval censoring.** Turnbull with annual or 6-monthly visits is unbiased (|bias| ≤ 0.0007). Turnbull interval coverage at n = 2380 was *not computed* (bootstrap only for n = 1000, 50 cohorts).

**Naive and wrong-definition estimators.** Treating censored records as event-free is biased low (−0.017 at 120 and −0.116 at 240 months). Midpoint imputation was essentially unbiased here (not a useful contrast). Analysing persistent exceedance when the target is first crossing biases $F$ by −0.0024 (120 months) and −0.0042 (240 months): small but about 10 MC standard errors; against the persistent target the bias is +0.0008, a pure definition effect. Early-heavy, late-heavy and mostly-no-crossing populations give |bias| ≤ 0.0013.

![Figure 6.13. Synthetic event-free survival recovery and population hyper-parameters](../results/figures/fig_p4_08_synthetic_event_free_survival_recovery.png)

**Figure 6.13.** Synthetic survival validation. *Setting:* truth = pooled P1 population; n = 2380 cohorts; the hyper-parameter part uses G1 with $\ln\beta\sim N(\ln 0.04,0.3^2)$, annual visits, 50 replicates; SYNTHETIC. *Axes/panels:* (a) cumulative event probability against time for one cohort: truth (black), Kaplan–Meier with administrative 240 + random loss (blue with band), naive censored-as-event-free (orange dashed), Turnbull annual (purple steps), midpoint imputation (green dotted); (b) bias of $F(h)$ (mean over 200 cohorts, 95% MC band) at 120 and 240 months for KM admin, KM random loss, naive, Turnbull annual and 6-month, midpoint, and KM persistent-only versus first crossing; (c) profile $\Delta\chi^2$ for $m=\text{mean}\ln\beta$ (blue) and $\ln s$ (orange) against estimate − truth with a dashed 3.84 line; (d) profile of $m_\beta$ with $r$ known (blue, extremely narrow) and with $r$ also varying (purple). *What it shows:* in (a) all estimators except the naive one follow the truth up to about 0.44 at 240 months, the naive curve ending near 0.32; in (b) all bars are near zero except naive (about −0.017 and −0.116) and persistent-only (about −0.003 and −0.004); in (c) both profiles are closed parabolas around zero (widths about ±0.02 for $m$, about ±0.06 for $\ln s$), with coverage 0.92 and 0.96; in (d) the $r$-known profile is a spike at 0 while the $r$-varying profile is flat from about +0.1 to +0.4 and below 3.84 from about −0.4 to +0.6, a ridge of about one log unit. *Mathematical interpretation:* censoring and interval censoring are handled correctly by KM/Turnbull; a monotone map $T(\beta)$ makes one heterogeneity source identifiable, but a second source creates a ridge. *Biological interpretation:* survival curves recover the distribution of event times, not the mechanisms producing them. *Limitations:* a single hypothetical population; 60-month results degenerate; informative censoring not studied.

**Population hyper-parameters.** For $\ln\beta\sim N(m,s^2)$ ($m=\ln 0.04$, $s=0.3$), n = 2380, annual visits, 50 replicates, all else known: $m$ bias $-0.0007$ (RMSE 0.0078; Wald coverage 0.94, profile 0.92), $\ln s$ bias $+0.0038$ (RMSE 0.022; coverage 0.96/0.96). Because $T(\beta)$ is monotone, the event-time law determines the $\beta$ law.

**A second varying parameter destroys it.** If $\ln r$ also varies (4 unknowns), 11 of 17 grid points of the $m_\beta$ profile lie inside the 95% region ($m_\beta\approx-3.62$ to $-2.62$, about one log unit, versus $\pm0.014$ with $r$ known). The best fit has $m_\beta=-2.73$ (truth $-3.22$) while the event-time CDF is fitted to a sup-norm error of 0.030 and moves by at most 0.009 along the ridge.

**Criterion B tested with nuisances left free** (`p4_population_nuisance*.csv`). Synthetic population with true $\tau=12$ months and initial-burden scale 1, n = 2380, annual visits to 240 months, one heterogeneity source, exact model structure:

- *$\tau$ free $(m,s,\tau)$:* over 30 replicates $\hat m=-3.217$ (truth $-3.219$; sd 0.024), $\hat s=0.300$ (sd 0.027), $\hat\tau=11.6$ months (sd 4.7); the replicate correlation of $\hat m$ and $\hat\tau$ is $-0.947$. Wrongly fixing $\tau=0$ biases $\hat m$ by $+0.051$ and $\hat s$ by $-0.062$.
- *$\tau$ and initial-burden scale $s_0$ both free (8 replicates):* $\hat m=-3.230$ (sd 0.055), $\hat s=0.310$ (sd 0.045); $m$ and $s$ keep closed intervals (replicate 0: $m\in[-3.38,-3.08]$, median $\beta$ within $\pm16\%$). But the nuisances are not identified: $\tau\in[-7,+46]$ months (sd 9.2 across replicates), the $s_0$ window ×0.22–×4.5 lies fully inside the 95% region, and $\hat\tau$ and $\ln\hat s_0$ correlate 0.88.
- *A second varying parameter* turns $m_\beta$ into a one-log-unit ridge even when $\tau$ and $s_0$ are known.

Criterion B as literally written (a finite profile interval for some population combination with $\tau$ and threshold/$N_0$ free) is therefore met only in an idealised setting of one known heterogeneity source and an exactly correct model. Which sources of heterogeneity exist in real patients cannot be checked from event times; the decision based on this judgement, made after seeing these results, is reported in Section 7.

### 6.9 Conclusions of Phase 4

1. **Structural vs practical.** All parameters are generically locally identifiable ($K=1$), but the low-burden regime reduces the information to a few combinations (rank 4 of 6 or 7 from $N$), and a burden-scale gauge exists when $K$ and $c$ are both unknown.
2. **Practical identifiability** from total burden: only $N_0$ (and the combination $\lambda_+$). Observing $P$, $D$ adds $r,q,a,f_D$; adding $e$ adds $\beta,\rho$. $\delta_D$ and $\eta$ are never practically identifiable.
3. **Good fits and good $T$ predictions coexist with parameters wrong by factors of 5–10;** prediction does not validate mechanism. Short follow-up (60/120 months) is the main prediction hazard (43.5%/26.8% error).
4. **Noise and sampling** are not the bottleneck; regime and observed compartments are.
5. **Model reduction** (immune frozen, $\delta_D=0$, $f_D$ fixed) is acceptable for burden-only sub-threshold prediction but rejected for $(P,D,e)$ data and nonlinear regimes; report $(\lambda_+,N_0)$ and treat the rest as nuisance.
6. **Nuisances confound event-time curves:** $\tau+d$ is an exact shift; up to about 12 months of offset cannot be separated from a factor-2 change in $N_0$ at n = 2380; thresholds from 0.01 to 0.5 move $F(240)$ from 0.69 to 0.21.
7. **Survival estimators** (KM, Turnbull) are unbiased under the stated censoring; event times identify a distribution, not mechanisms.
8. **Remaining limitations:** local generic structural analysis only; practical results conditional on the assumed observation model, six hypothetical truths and 6–8 replicates; ridge-like likelihoods with multistart fits; Wald diagnostics poorly calibrated far from the truth; Turnbull coverage at n = 2380 not computed.

These findings feed the feasibility decision in Section 7.
