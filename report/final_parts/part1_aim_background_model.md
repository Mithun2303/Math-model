## 1. Introduction, Aim and Objectives

### 1.1 The clinical problem

A patient treated for breast cancer is often declared disease-free, yet some patients relapse many years later. For oestrogen-receptor-positive (ER+) disease the recurrence rate after year 5 does not vanish: the EBCTCG meta-analysis (Pan et al. 2017; abstract-level evidence only) reports that recurrences "occurred at a steady rate throughout the study period from 5 to 20 years" in 62,923 women who were disease-free after 5 years of endocrine therapy. Rueda et al. (2019; abstract-level) describe late-recurring ER+/HER2− subgroups with a high risk of recurrence up to 20 years after diagnosis. A common biological explanation is **tumour dormancy**: residual disease persists, clinically silent, for years to decades (Aguirre-Ghiso 2007; abstract-level), held in check by cellular quiescence (G0–G1 arrest), by limited blood supply, or by immune surveillance.

Two features make this a mathematical problem and not only a clinical one. First, relapse after a long silent interval is a *threshold-crossing* phenomenon: a tiny residual population must grow past a detectable size. Second, the clinical data that exist (Kaplan–Meier curves of relapse-free survival) record only *when* the threshold is crossed, never the hidden cell populations. Whether a mechanistic model of dormancy can say anything quantitative about such data is the question this project investigates.

### 1.2 Aim

The aim of the project is to build, analyse and stress-test a deterministic three-compartment ordinary-differential-equation model of **tumour dormancy and immune control** (proliferating cells *P*, dormant cells *D*, effective immune activity *E*), and to determine rigorously

1. what the model does mathematically (equilibria, stability, bifurcations, basins, relapse times);
2. which parameters control relapse timing, and how uncertain simulated relapse probabilities are under stated assumptions;
3. whether the model's parameters can be recovered from the kinds of data that exist (identifiability), and
4. whether a clinical comparison with a public breast-cancer cohort (METABRIC) is feasible, and at what level (descriptive, phenomenological or mechanistic).

**What is explicitly not claimed.** All parameter values are hypothetical and illustrative; every probability produced is a *model-defined simulated quantity*, not a patient risk; nothing here is a clinical prediction. Negative and limiting results are treated as first-class outputs (research_gap_and_objectives.md §4).

### 1.3 Research questions and hypotheses: final verdicts

The research questions and falsifiable hypotheses were fixed at the end of Phase 1 (research_gap_and_objectives.md §6). Table 1.1 gives each one together with its verdict after Phases 1–4. The verdicts are summaries; the evidence is developed in the later parts of this report.

**Table 1.1.** Research questions, hypotheses and final verdicts (source: research_gap_and_objectives.md §6; PHASE_2_SUMMARY.md; PHASE_3_SUMMARY.md; PHASE_4_SUMMARY.md).

| RQ | Hypothesis as pre-stated | Verdict | Basis |
|---|---|---|---|
| RQ1 What bounds immune control of residual disease? | H1: the tumour-free state is locally stable iff $r < ks/\mu + q\delta_D/(a+\delta_D)$ | **Supported (proved, local).** Equivalent to $\beta < r_\mathrm{eff}$, with $k_c = 0.28333$ at baseline. Qualification: *local* stability does not guarantee control from every initial state: in the bistable cases of Phase 3 the tumour-free basin is bounded. | Phase 2 derivation (Sotomayor conditions); Phase 3 basins |
| RQ2 Does dormancy alone prevent escape? | H2: if $\delta_D=0$ the boundary reduces to $r<ks/\mu$, so switching alone does not change tumour-free stability | **Supported algebraically** ($r_\mathrm{eff}=r-q\delta_D/(a+\delta_D)$ equals $r$ at $\delta_D=0$). The sources read do not report a separate numerical sweep over $(q,a)$. | Phase 1 algebra; Phase 2 §0 |
| RQ3 How is relapse time controlled? | H3: $T\approx\ln(N_\mathrm{det}/N_0)/\lambda_+$ | **Not supported in this crude form.** At $k=0.5k_c$ the estimate is 237.5 months vs 375.5 simulated. The gap (138 months) is explained by eigenmode projection (+44.0), logistic saturation (+11.5) and immune feedback during growth (+82.5; order-dependent split). Refined estimators A1–A3 were built; at $N_\mathrm{det}=0.1$ the best closed form is A1 (10% median error). | Phase 2 §6 |
| RQ4 Can relapse likelihood be estimated? | H4: the Monte Carlo cumulative event curve depends strongly on the prior (5-year probability differing by more than a factor of 2 between reasonable priors) | **Supported, with a caveat.** Across seven priors P(5 y) ranges 0–0.036, P(10 y) 0.034–0.232, P(20 y) 0.317–0.468; the prior spread (0.15 at 20 years) is 5.5× the Monte Carlo interval width, and the threshold choice moves P(20 y) from 0.689 to 0.214. The pooled 5-year probability under the primary prior is very small (8/15,000). | Phase 3 §5–8 |
| RQ5 Can the model reproduce group-level recurrence patterns? | H5: a prior exists such that the simulated curve lies in the METABRIC Kaplan–Meier band at 5 and 10 years (not evidence of mechanism) | **Tested descriptively in Phase 5 (pre-registered, no fitting).** Reference-scenario simulated incidence is below the METABRIC 95% CI at 5 and 10 years in all or nearly all of 1,764 scenarios; agreement at 20 years occurs in 11.1% of scenarios. A prior in the METABRIC band at both 5 and 10 years was not found (0% inside the CI at 5 y). See Section 12. | Phase 4 §9 |
| RQ6 Is the P/D split identifiable from N alone? | H6: no | **Supported.** From total burden alone, only the initial burden $N_0$ is practically identifiable; the low-burden linear regime has symbolic rank 4 of 6; 62.5% of base-condition fits ended on a parameter bound while the crossing time was still predicted within 10% in 87.5% of replicates. | Phase 4 §2–4 |

### 1.4 Objectives and their status

**Table 1.2.** Objectives O1–O7 with success criteria and final status (source: research_gap_and_objectives.md §5; Phase summaries).

| # | Objective | Phase | Status |
|---|---|---|---|
| O1 | Verify and extend the equilibrium/stability analysis; locate and classify a second boundary or interior saddle | 2 | **Done.** Transcritical point verified (Sotomayor), at most three nonzero equilibria (cubic), saddle-node folds and two kinds of bistability found; no Hopf found in scanned ranges. |
| O2 | Define operational dormancy and relapse; show sensitivity of conclusions to the thresholds | 2–3 | **Done.** Event = first upward crossing of $N_\mathrm{det}$; thresholds 0.01–0.5 and a dormancy-duration definition were varied. |
| O3 | Quantify parameter influence (local and global) | 3 | **Done.** Local elasticities and flip margins; Spearman/PRCC, Morris and Sobol (45,056 runs). |
| O4 | Relapse-time distributions at 5/10/20 years with at least two prior scenarios | 3–4 | **Done** (seven prior scenarios, three independent runs pooled for the primary prior). |
| O5 | Structural identifiability by synthetic recovery | 4 | **Done** (exact structural analysis, profile likelihoods, 184 constrained fits). |
| O6 | Group-level comparison with METABRIC RFS at 5, 10 (cautiously 15) years | 5 | **Partly done.** Dataset audit and endpoint feasibility were completed (Phases 1 and 4); decision: descriptive comparison only. No model-versus-data overlay was produced. |
| O7 | Ablation: logistic vs P/D vs P/D/E | 4–5 | **Only partly addressed.** Phase 4 tested reductions of the full model (immune frozen, $\delta_D=0$, fixed dormant fraction), not the planned logistic-versus-P/D-versus-P/D/E comparison on a common clinical target. |

### 1.5 How to read this report

The report follows the order in which the work was done, and each step is justified by the one before it. Section 2 reviews the literature and states the research gap. Section 3 (this part) defines the model, its parameters, its assumptions and the project workflow. The following parts present, in order, the mathematical analysis of the model (equilibria, bifurcations, normalisation and relapse-time approximations; Phase 2), sensitivity, basins and Monte Carlo simulation (Phase 3), identifiability and synthetic validation (Phase 4), and the clinical dataset investigation and conclusions. Every figure carries a caption with its setting, axes, what it shows, and a mathematical and a biological interpretation, followed by its limitations. Throughout, three kinds of statement are kept apart: **proved** (analytic derivation or exact computation), **numerical observation** (simulation or scan within stated ranges, which is not a proof of absence or generality) and **hypothetical** (anything depending on the assumed parameter values or priors). Figures are in `results/figures/`, tables in `results/tables/`, and the figure and table captions name their source files.

---

## 2. Literature Review and Research Gap

### 2.1 How the review was done and what was actually read

Full texts were not equally accessible, so every source carries an explicit verification level (literature_review.md §0). **Full text (FT)** means equations, methods and results were read; **abstract (ABS)** means only the abstract and bibliographic metadata were read, and anything beyond them is unverified. Of the five required papers only Wilkie & Hahnfeldt 2013 (*Cancer Research*) was read in full text, from an open author manuscript. Wilkie & Hahnfeldt 2013 (*Interface Focus*), Yáñez Feliú et al. 2026, Blath et al. 2025 and Rueda et al. 2019 are abstract-level (PubMed/PMC blocked automated access; Nature is paywalled; the 2026 article body was not served). Consequently **equation-level statements are made only for Wilkie & Hahnfeldt 2013 (*Cancer Research*)**. The search used Europe PMC title/topic queries and the reference list of that paper; it is bounded and not a systematic review.

### 2.2 Biological background

*Dormancy is an umbrella term.* Aguirre-Ghiso (2007, ABS) summarises three mechanisms, **angiogenic dormancy, cellular dormancy (G0–G1 arrest) and immunosurveillance**, with latency "from years even to decades". The P–D–E model combines cellular state switching (P↔D) with immune control (E); it has no angiogenic mechanism.

*Late relapse is quantitatively substantial.* Pan et al. (2017, ABS): 5–20-year distant-recurrence risk in ER+ disease depends strongly on tumour size and nodal status (13% for T1N0, 34% for T1N4–9, 41% for T2N4–9). The consequence for modelling is that the empirical hazard after year 5 is non-vanishing, which a model with absorbing "cure" cannot reproduce unless its parameter distribution contains slow-escape members.

*Relapse dynamics differ by subgroup.* Rueda et al. (2019, ABS) fit a multistate model (locoregional recurrence, distant recurrence, breast-cancer death, other death) in 3,240 patients (1,980 with molecular data), reporting four late-recurring ER+/HER2− subgroups (about 26% of such tumours; mean recurrence risk 47–62% up to 20 years) and a triple-negative subgroup that rarely recurs after five years. This motivates a model able to generate slow, heterogeneous escape times. The multistate model is richer than the single composite endpoint in the public cBioPortal export (Section 13 of this report).

### 2.3 The closest mathematical precedent: Wilkie & Hahnfeldt 2013 (*Cancer Research*), full text

This paper models cancer cells *C*, immune cells *I* and a dynamic environmental carrying capacity $K_C$:

$$\frac{dC}{dt}=\mu\alpha\,(1+\Psi(I,C))\,C\left(1-\left(\frac{C}{K_C}\right)^{\alpha}\right),\qquad \frac{dI}{dt}=\lambda (I+rC)\left(1-\frac{I}{K_I}\right),$$

with a saturating predation term $\Psi(I,C)=-\theta\big(I^\beta/(\phi C^\beta+I^\beta)+\varepsilon\log_{10}(1+I)\big)$, $\varepsilon=0.01$ assumed, and $dK_C/dt = pC - qK_C C^{2/3}$ (from Hahnfeldt et al. 1999). They analyse equilibria and Jacobians analytically: four equilibria, with an interior point that is a **saddle** when both boundary equilibria are stable, and two transcritical bifurcations in the immune-efficacy parameter θ. Parameters were estimated by MCMC from murine fibrosarcoma growth data under the assumption of no immune predation. Findings relevant here: (i) dormancy is "a transient period of growth that necessarily ends in either tumor elimination or tumor escape" (the dormant interior state is a saddle); (ii) four parameter sets fit untreated growth equally well but predict very different outcomes under the same immunotherapy, a direct warning about non-identifiability from growth curves alone; (iii) their operational dormancy definition (at least 30 days within the first 150 days with $(1/C)\,dC/dt<0.04$ and diameter 1–30 mm) shows that any dormancy definition is operational. Differences from our model: their immune recruitment is tumour-triggered with $I(0)=0$ in the treatment simulations, so no immune-surveilled tumour-free state exists, whereas ours has a constant replenishment $s$ giving the tumour-free level $E^*=s/\mu$; they have a dynamic carrying capacity and no explicit quiescent compartment; they use mouse data and make no time-to-event comparison with patients.

### 2.4 Other modelling literature (abstract-level unless stated)

**Table 2.1.** Supporting models and what their abstracts establish (source: literature_review.md §2.2–2.5).

| Source | What the abstract establishes | Relevance to P–D–E |
|---|---|---|
| Wilkie & Hahnfeldt 2013 (*Interface Focus*) | Mouse data: longer dormancy gives more resistance to T-cell killing; decline in immune recruitment is a stronger determinant of escape than increased resistance; immune-resistant subpopulations and immunoediting analysed | Supports time-varying $k$ or $s/\rho$ as a future extension; escape can be immune failure or tumour adaptation, which constant parameters cannot distinguish |
| Yáñez Feliú et al. 2026 | Evolutionary game with environmental conditioning factor; history dependence; large variation in time to outgrowth near thresholds | Motivates an environmental signal driving reactivation; the constant reactivation rate $a$ is a strong simplification |
| Blath et al. 2025 | Stochastic agent-based model of spontaneous and treatment-induced switching into dormancy; a small dormant population can cause therapy failure | Supports P↔D switching and a stochastic extension; concerns therapy response, not post-treatment latency |
| Kuznetsov et al. 1994 | Cytotoxic T-cell/tumour model with a dormant state, bifurcations computed, parameters estimated for mouse B-lymphoma | Classic bifurcation-analysed precedent; equation-level comparison not done |
| Page & Uhr 2005 | Immune-induced cell-cycle arrest; parameters controlling dormant tumour size and proliferating fraction | Nearest precedent for a proliferating-versus-arrested split in an immune model |
| Wilkie 2013 (review) | Quiescent-state transition is a recognised augmentation of predator–prey dormancy models | The quiescent compartment itself is *not* novel |
| de Pillis et al. 2005 | NK and CD8⁺ T-cell ODE model validated against mouse and human data | A data-validated route to immune parameter ranges; ours are unvalidated |
| Enderling et al. 2009 | Cell-based model; increased death can paradoxically accelerate long-term growth | Space and stem-cell effects are outside a well-mixed ODE |
| Mehdizadeh et al. 2021 | A dormant-like population that is *more* immunogenic in murine 4T1 | Counter-evidence for the assumption that immune killing acts only on P |
| Saudemont et al. 2005 | Persistent dormant leukaemic cells resist CTL lysis but are killed by NK cells | Supports cell-state-dependent immune susceptibility |
| Baffoe & Bouchnita 2026 | Dormancy as a continuous ERK/p38 phenotype axis | Direct critique of the binary P/D simplification |

### 2.5 Clinical and data sources

Curtis et al. (2012) and Pereira et al. (2016) describe the METABRIC discovery/validation cohorts with long-term follow-up (abstract-level); Pan et al. (2017) and Rueda et al. (2019) provide the late-hazard targets; Tolaney et al. (2021, STEEP 2.0) justifies insisting that "RFS", "DFS" and "OS" are not interchangeable endpoints.

### 2.6 Evaluation of the model against the literature

**Table 2.2.** Assessment of the model assumptions (source: research_gap_and_objectives.md §2).

| Question | Assessment |
|---|---|
| Is the P/D separation defensible? | Defensible as a simplification, not as a measured partition. Cellular dormancy and state switching are established (Aguirre-Ghiso 2007; Blath et al. 2025; Page & Uhr 2005), but binary dormancy is contested (Baffoe & Bouchnita 2026). |
| Is the immune interaction reasonable? | Reasonable as an aggregate with two caveats: mass-action killing has no saturation (Wilkie & Hahnfeldt use a saturating ratio-dependent term), and the constant source $s$ gives a surveilled tumour-free state. |
| Should killing act only on P? | Contested (Mehdizadeh et al. 2021; Saudemont et al. 2005); the assumption is one of several literature-consistent options. |
| What do q and a represent? | First-order switching rates with mean dwell times $1/q$ and $1/(a+\delta_D)$; population-averaged mathematical rates, not measured quantities. |
| Which mechanisms are omitted? | Environmental feedback on carrying capacity or awakening, immunoediting, angiogenic dormancy, space, treatment, stochastic extinction of small populations, time-varying parameters, delays. |
| Which parameters have support? | None quantitatively for human post-treatment disease. |

### 2.7 The parameter evidence audit

Each of the 11 model parameters was audited for sourced support (parameter_evidence.md). **The result is that all 11 parameters are hypothetical for human post-treatment residual disease.** The only numerical values actually retrieved are from the murine fibrosarcoma model of Wilkie & Hahnfeldt (2013): an intrinsic growth rate of about 0.115 per day (≈ 3.5 per month), about **30× faster** than the baseline $r=0.12$ per month, and with a different growth law and predation form, so not transferable. One qualitative plausibility check exists: with $a=0.005$ and $\delta_D=0.001$ the mean dormant dwell time $1/(a+\delta_D)\approx167$ months ($\approx13.9$ years) is consistent in order of magnitude with latencies of "years even to decades". The same audit derived, from the equations alone, four structural observations that shaped later phases: (i) a scale redundancy between $k$, $s$ and the unit of $E$; (ii) at low burden only the product $ks/\mu$ matters, so $\rho$, $h$, $\eta$ are invisible at first order; (iii) with only total burden observed, the linearised tumour subsystem cannot determine its four matrix entries; (iv) Kaplan–Meier event times carry only a distribution of threshold-crossing times.

### 2.8 Research gap and claimed contribution

Based on the sources read (bounded, non-systematic):

1. Dormancy–immune ODE models are mostly analysed for **qualitative regimes** (elimination, dormancy, escape) and parameter-window bifurcations. No abstract-level evidence was found that they are converted into **model-based time-to-event distributions** and compared with **right-censored recurrence data**.
2. Clinical late-recurrence work characterises hazards statistically; it does not propose a mechanistic generative model of the dormant compartment.
3. **Identifiability and endpoint-definition mismatch** (competing risks, time origin, composite events) is seldom confronted when linking mechanistic models to survival data.
4. Dormancy is "a stable state" in classic models but "transient" in Wilkie & Hahnfeldt; consequences for relapse-time *distributions* are not examined in the abstracts read.
5. Binary versus continuous dormancy, and the immune susceptibility of dormant cells, are unresolved.

The contribution that can legitimately be claimed is therefore: an interpretable P–D–E framework with a closed-form tumour-free stability boundary and verified bifurcation structure; a protocol turning a deterministic dormancy ODE into relapse-time distributions with explicit uncertainty and censoring; a transparent assessment of what group-level clinical data can and cannot constrain, with the mismatches (time origin, composite event, competing risks) quantified; and negative or limiting results reported as outputs. **Not claimed:** clinical validity, patient-level prediction, identification of dormant-cell counts, or novelty of the quiescent compartment itself.

### 2.9 Caveats of the review

Four of the five required papers were abstract-level; abstract-level claims should be re-checked against full texts before being quoted beyond the abstract. Full texts of Rueda 2019, Wilkie & Hahnfeldt 2013 (*Interface Focus*), Yáñez Feliú 2026, Blath 2025, Kuznetsov 1994 and de Pillis 2005 remain unread.

---

## 3. The Mathematical Model

### 3.1 State variables and equations

The model has three state variables, with time in months and tumour burdens normalised to the carrying capacity $K$ (so $K=1$ means "fraction of carrying capacity"):

- $P$: proliferating tumour burden;
- $D$: dormant tumour burden, with $N=P+D$ the total burden;
- $E$: effective immune activity (arbitrary effector units).

$$\begin{aligned}
\frac{dP}{dt}&=rP\left(1-\frac{N}{K}\right)-qP+aD-kEP\\
\frac{dD}{dt}&=qP-(a+\delta_D)D\\
\frac{dE}{dt}&=s+\rho\,\frac{P}{h+P}\,E-\mu E-\eta N E
\end{aligned}\tag{3.1}$$

Term by term: $rP(1-N/K)$ is logistic growth of proliferating cells limited by the *total* burden; $-qP$ moves cells to dormancy and $+aD$ returns them (reactivation); $-kEP$ is mass-action immune killing of proliferating cells only; $qP-(a+\delta_D)D$ gives dormant cells an inflow from $P$ and an outflow by reactivation $a$ and loss $\delta_D$; in the immune equation, $s$ is a constant replenishment (immune surveillance), $\rho P/(h+P)\,E$ is tumour-stimulated expansion saturating at half-saturation burden $h$, $-\mu E$ is turnover, and $-\eta N E$ is burden-associated suppression.

### 3.2 Assumptions

The following assumptions are built into the baseline (research_gap_and_objectives.md §2 and §8; Phase 1 summary). They are recorded as limitations, not defects to be fixed in the baseline.

1. Binary dormancy: cells are either proliferating or dormant; no continuous dormancy depth.
2. Immune killing is mass-action ($kEP$) with no saturation, acting only on $P$; dormant cells are immune-invisible.
3. The carrying capacity $K$ is constant; no environmental feedback.
4. All parameters are constant in time; no treatment, delays, immunoediting or time-varying immune efficacy.
5. The model is deterministic and well mixed (no space, no stochastic extinction of small populations).
6. The immune system has a constant source $s$, so the tumour-free state is an immune-surveilled equilibrium $E^*=s/\mu$.
7. Relapse is defined operationally by a detection threshold (design choice, not an imaging limit).
8. Initial state represents hypothetical post-treatment residual disease, mostly dormant.

### 3.3 Parameters

**Table 3.1.** Baseline parameters (source: configs/baseline.yaml; parameter_evidence.md). All values are hypothetical illustrative choices, not fitted estimates.

| Symbol | Meaning | Units | Baseline |
|---|---|---|---|
| $K$ | carrying capacity (normalised) | burden | 1.0 |
| $r$ | proliferation rate | month⁻¹ | 0.12 |
| $q$ | proliferating → dormant rate | month⁻¹ | 0.04 |
| $a$ | dormant → proliferating (reactivation) rate | month⁻¹ | 0.005 |
| $\delta_D$ | dormant-cell loss rate | month⁻¹ | 0.001 |
| $k$ | immune killing effectiveness | (E·month)⁻¹ | 0.1 |
| $s$ | baseline immune replenishment | E·month⁻¹ | 0.02 |
| $\mu$ | immune turnover | month⁻¹ | 0.05 |
| $\rho$ | tumour-stimulated immune expansion | month⁻¹ | 0.04 |
| $h$ | half-saturation of immune stimulation | burden | 0.02 |
| $\eta$ | burden-associated immune suppression | (burden·month)⁻¹ | 0.08 |

Initial conditions: $P_0=10^{-6}$, $D_0=10^{-4}$ (so the residual disease is about 99% dormant) and $E_0=s/\mu=0.4$ (the tumour-free immune level). The baseline $r=0.12$ gives a doubling time $\ln 2/r\approx5.8$ months, and $1/(a+\delta_D)\approx167$ months is the mean dormant dwell time.

### 3.4 Derived quantities and the tumour-free stability condition

The following derived constants are used throughout (Phase 2 report §0).

**Table 3.2.** Derived constants (source: report/phase_2_mathematical_analysis.md §0; parameter_evidence.md §4).

| Symbol | Definition | Baseline |
|---|---|---|
| $c$ | $q/(a+\delta_D)$ (equilibrium $D/P$ ratio) | 6.667 |
| $\varphi$ | $1+c=N^*/P^*$ at every nonzero equilibrium | 7.667 |
| $r_\mathrm{eff}$ | $r-q\delta_D/(a+\delta_D)$ | 0.11333 |
| $E_f$ | $s/\mu$ (tumour-free immune level) | 0.4 |
| $k_c$ | $r_\mathrm{eff}\,\mu/s$ | 0.28333 |
| $\beta$ | $ks/\mu$ (effective killing) | 0.04 at $k=0.1$ |
| dormant fraction | $q/(q+a+\delta_D)$ at nonzero equilibria | 0.8696 |

**The tumour-free equilibrium** is $(P,D,E)=(0,0,s/\mu)$: setting $P=D=0$ in (3.1) leaves $\dot E=s-\mu E$. At that point the $E$-column of the Jacobian vanishes in the $P$ and $D$ rows, so the Jacobian is block lower-triangular, with eigenvalues $-\mu$ and those of the $2\times2$ block

$$M=\begin{pmatrix} r-q-\beta & a\\ q & -(a+\delta_D)\end{pmatrix},\qquad \beta=\frac{ks}{\mu}.$$

Since $\operatorname{tr}M<0$ whenever $\det M>0$, the state is locally asymptotically stable exactly when $\det M>0$, i.e. (proved, Phase 2 §3.2)

$$r<\frac{ks}{\mu}+\frac{q\,\delta_D}{a+\delta_D}\quad\Longleftrightarrow\quad \beta>r_\mathrm{eff}\quad\Longleftrightarrow\quad k>k_c=\frac{r_\mathrm{eff}\,\mu}{s}.\tag{3.2}$$

This is hypothesis H1 and it is a statement about *local* stability only. At baseline $k_c=0.28333$ with the other two eigenvalues $-0.05$ and $-0.0393$ at $k=k_c$. If $\delta_D=0$ the condition reduces to $r<ks/\mu$ (H2). The role of dormancy is visible in (3.2): the term $q\delta_D/(a+\delta_D)$ is the growth lost to dormant-cell death, so dormancy helps control only through $\delta_D>0$.

### 3.5 Event definition, censoring and solver settings

*Relapse* is the first **upward** crossing of the burden $N=P+D$ through a model-defined detection threshold $N_\mathrm{detect}=0.1$ (fraction of $K$). Follow-up has a horizon of 240 months (20 years); a trajectory that never crosses within the horizon is **right-censored** at 240 months, and no relapse time is invented for it. The event is detected as a terminal event of the integrator (direction $+1$), so the time is the exact crossing time and not a grid value. A subtlety recorded for later phases: if $N_0\ge N_\mathrm{detect}$ the legacy solver returns an event at time 0; later phases exclude such samples as having no defined event.

Integration uses SciPy's `solve_ivp` with the explicit Runge–Kutta method DOP853, relative tolerance $10^{-9}$ and absolute tolerance $10^{-12}$ (configs/baseline.yaml). Later phases verify results against other solvers (Radau, LSODA) and tolerances from $10^{-6}$ to $10^{-12}$.

### 3.6 Baseline behaviour (Figure 3.1)

Phase 1 ran four scenarios with the killing strength set to fractions of the critical value $k_c$ (outputs in `results/tables/baseline_scenarios.csv`). Figure 3.1 shows the result: one scenario produces a relapse within 240 months; the three stronger-killing scenarios do not.

![Figure 3.1. Baseline trajectories at four immune killing strengths](../results/figures/baseline_trajectories.png)

**Figure 3.1.** Baseline trajectories at four immune killing strengths. *Setting:* baseline parameters of Table 3.1 except $k$; initial state $P_0=10^{-6}$, $D_0=10^{-4}$, $E_0=0.4$; integration to 240 months or to the first crossing of $N_\mathrm{detect}=0.1$. *Axes/panels:* four panels (weak $0.25k_c$, $k=0.071$; moderate $0.5k_c$, $k=0.142$; near boundary $0.98k_c$, $k=0.278$; strong $2k_c$, $k=0.567$); horizontal axis time in months (0–250), vertical axis burden on a logarithmic scale from $10^{-8}$ to $10^{0}$; blue $P$, orange $D$, black $N=P+D$, dashed red line $N_\mathrm{detect}$. *What it shows:* in the weak-killing panel $N$ grows from about $10^{-4}$ and crosses the threshold at 172 months, where the curve ends. At $0.5k_c$ the burden rises to about $2\times10^{-2}$ by 240 months without crossing. At $0.98k_c$ the burden stays almost flat near $10^{-4}$, and at $2k_c$ it decays slowly to about $3\times10^{-5}$. *Mathematical interpretation:* below $k_c$ the tumour-free state is unstable and growth is exponential at the rate of the dominant eigenvalue after a short transient in which $P$ rises quickly from $10^{-6}$ towards the dormant-dominated ratio; the growth rate shrinks as $k\to k_c$ and becomes negative above it, in agreement with (3.2). *Biological interpretation:* the same residual disease can stay silent for the whole follow-up (or decline) or relapse after about 14 years, depending only on immune killing. *Limitations:* a four-point illustration at hypothetical parameters; "no event" at $0.5k_c$ is a censoring at 240 months, not control (Phase 2 later shows that this trajectory does cross 0.1 at 375.5 months in a longer run); it does not show equilibria or basins.

The case at $0.5k_c$ is the origin of the project's first analytic problem: the crude estimate $\ln(N_\mathrm{det}/N_0)/\lambda_+$ gives about 240 months, while the longer simulation gives 375.5 months (hypothesis H3, Table 1.1).

### 3.7 Project workflow: how each phase led to the next

**Table 3.3.** Project phases, what each produced and what it implied for the next (source: PHASE_1_SUMMARY.md–PHASE_4_SUMMARY.md).

| Phase | Question | Main outputs | Led to |
|---|---|---|---|
| 1. Literature, evidence and dataset feasibility | Is the model defensible, and what do parameters and data support? | Literature review with verification levels; parameter audit (all hypothetical); structural observations; METABRIC and TCGA audits; **conditional GO** for a group-level exploratory comparison only | A failed first check of H3 and an observed transcritical-like behaviour at $k_c$, both needing proof; the $E$-gauge redundancy |
| 2. Formal mathematical analysis | What exactly does the model do? | Exact equilibrium reduction (at most three nonzero equilibria); transcritical bifurcation verified; saddle-nodes and two kinds of bistability; no Hopf found; normalised model (M1); relapse-time estimators A0–A3 and the $k_\mathrm{det}$ boundary | A scale-free parameterisation and the knowledge that relapse depends on a few combinations, so global sensitivity and uncertainty analysis became meaningful |
| 3. Sensitivity, basins and uncertainty | Which inputs matter, and how uncertain is simulated relapse? | Local/global sensitivity; basins of attraction; Monte Carlo event probabilities under seven priors; threshold robustness | The finding that assumptions (priors, threshold) dominate sampling error, so any clinical calibration must first answer what can be identified |
| 4. Identifiability, synthetic validation and endpoint mapping | Can parameters be recovered, and can data constrain them? | Exact structural and practical identifiability; synthetic recovery; model reduction; survival validation; time-origin mapping; METABRIC endpoint audit; pre-registered feasibility decision | Decision: descriptive comparison (A) is feasible; restricted calibration (B) not supported; mechanistic calibration (C) unsuitable |
| 5 (planned) | Descriptive METABRIC Kaplan–Meier analysis | A script exists in the repository; no outputs or summary were found | Not part of the conclusions of this report |

The logic of the sequence is: Phase 1 showed that no parameter has empirical support and that a naive relapse-time formula fails; Phase 2 therefore established the exact mathematical structure and the correct approximations; Phase 3 used that structure to ask which parameters matter and how far simulated probabilities depend on assumptions; Phase 4 asked whether the parameters could be recovered at all, and whether the data available could constrain the model. The pre-registered decision criteria in Phase 4 were written before the results were seen (configs/phase4.yaml, config hash `cfb17cec5c80`); Phase 3 configuration hash is `c7e70720ebfc`, and every Phase 3 and Phase 4 table and figure is stamped with it.

### 3.8 Repository and code structure

The project is a Python package with experiment scripts, configurations, tests and written reports (README.md).

**Table 3.4.** Repository layout (source: README.md; PHASE_n_SUMMARY.md).

| Location | Content |
|---|---|
| `src/dormancy/` | `model.py` (right-hand side, analytic Jacobian, second-derivative form), `parameters.py`, `solver.py` (first-upward-crossing relapse event with right-censoring), `equilibria.py` (exact cubic reduction, enumeration), `stability.py` (Routh–Hurwitz, labels), `continuation.py`, `bifurcation.py`, `relapse.py` (estimators A0–A3), `normalized.py`; Phase 3: `uncertainty.py`, `sensitivity.py`, `basins.py`; Phase 4: `identifiability.py`, `observation.py`, `synthetic_validation.py`, `endpoint_mapping.py`; `km_descriptive.py` |
| `configs/` | `baseline.yaml`, `normalized.yaml`, `phase3.yaml`, `phase4.yaml` (priors, seeds, pre-registered criteria) |
| `experiments/` | scripts that regenerate every table and figure (for example `baseline.py`, `bifurcation_analysis.py`, `global_sensitivity.py`, `monte_carlo_relapse.py`, `identifiability_analysis.py`, `synthetic_recovery.py`, `metabric_endpoint_investigation.py`) |
| `tests/` | the automated test suite |
| `data/` | raw METABRIC and TCGA-BRCA clinical exports, data dictionary, feasibility assessment |
| `results/tables/`, `results/figures/` | all numerical outputs and figures |
| `report/`, `PHASE_n_SUMMARY.md` | detailed phase reports and summaries |

Design rules followed throughout: each phase adds code and does not modify earlier modules, configurations or results (the baseline configuration is unchanged; `baseline_scenarios.csv` regenerates byte-identical); random sampling is seeded and reproducible (parallel runs equal serial runs exactly); and Phase 3 and 4 outputs carry the hash of their configuration file.

### 3.9 Tests

The automated test suite grew with the project and no earlier test was changed:

**Table 3.5.** Test counts by phase (source: phase summaries).

| After phase | Tests passing | New tests in that phase | Examples of what is tested |
|---|---|---|---|
| 1 | 16 | 16 | model, solver, equilibria and stability code (test list not itemised in the Phase 1 summary; the Phase 1 stability boundary was checked against eigenvalues) |
| 2 | 71 | 55 | continuation against the exact branch, cubic vs multistart equilibria, Sotomayor conditions, normalised-model equivalence over gauges, relapse estimators |
| 3 | 111 | 40 | LHS stratification, copula correlation, dense runner vs verified solver, start-above-threshold exclusion, Sobol on the Ishigami function, basin classification |
| 4 | 138 | 27 | Taylor series vs a $10^{-13}$ ODE solution, forward sensitivities vs finite differences, exact rank, Tobit and censoring, noise-free recovery within 1%, Kaplan–Meier worked example, Turnbull |

All tests passed at the end of each phase (0 failures). Passing tests establish that the code implements the stated mathematics; they do not make the biological assumptions or parameter values true.

### 3.10 Summary of this part

The model (3.1) is a hypothetical three-compartment system whose local tumour-free stability is governed by the single condition (3.2). The literature supports each ingredient qualitatively but offers no numerical parameter for human post-treatment disease, so every conclusion in this report is a conclusion about the *structure* of the model and about what can or cannot be inferred from it, not about patients. The next part develops the exact mathematical analysis that Phase 1 called for.
