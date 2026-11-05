# Literature Review: Tumour Dormancy, Relapse and Mathematical Modelling

*Phase 1 deliverable. Written 2026-11-04. Citation keys refer to `report/references.bib`.*

## 0. How this review was done, and what was actually read

Full texts are not equally accessible. To avoid presenting abstract-level knowledge as full-paper knowledge, every source carries an explicit **verification level**:

| Level | Meaning |
|---|---|
| **FT** | Full text read (equations, methods, results). |
| **ABS** | Abstract and bibliographic metadata read only. Anything beyond the abstract is marked *not verified*. |

| Source | Key | Level | How retrieved |
|---|---|---|---|
| Wilkie & Hahnfeldt 2013, *Cancer Res* | `Wilkie2013CancerRes` | **FT** | NCBI E-utilities (NIH author manuscript, PMC3955200) |
| Wilkie & Hahnfeldt 2013, *Interface Focus* | `Wilkie2013Interface` | ABS | PMC full text blocked / not open access via API |
| Yáñez Feliú et al. 2026, *npj Syst Biol Appl* | `YanezFeliu2026` | ABS | Publisher page returned abstract and metadata but not the article body |
| Blath et al. 2025, *J Theor Biol* | `Blath2025` | ABS | Europe PMC |
| Rueda et al. 2019, *Nature* | `Rueda2019` | ABS (+ reference list) | Publisher page (paywalled) |
| All other papers | see bib | ABS | Europe PMC metadata |

PubMed/PMC pages block automated fetching, and several publishers withhold the body text, so **equation-level statements below are made only for Wilkie & Hahnfeldt 2013 (Cancer Res)**. For the others I report what the abstract states and list what remains to be checked in the full text. Direct quotations are short and attributed.

---

## 1. Biological background

**Dormancy is an umbrella term for distinct mechanisms.** Aguirre-Ghiso's review (`Aguirre2007`, ABS) frames cancer dormancy as "a stage in cancer progression in which residual disease is present but remains asymptomatic", with latency periods "that range from years even to decades", and summarises experimental and clinical evidence for three mechanisms: **angiogenic dormancy, cellular dormancy (G0–G1 arrest) and immunosurveillance**. This matches the taxonomy in our project plan (cellular, tumour-mass, immune-mediated). Our P–D–E model combines cellular state-switching (P↔D) with immune control (E), and has no angiogenic mechanism.

**Late relapse is real and quantitatively substantial in ER-positive breast cancer.** The EBCTCG meta-analysis (`Pan2017`, ABS; 88 trials, 62,923 women with ER-positive disease who were disease-free after 5 years of endocrine therapy) reports that "breast-cancer recurrences occurred at a steady rate throughout the study period from 5 to 20 years", and that the 5–20-year distant-recurrence risk depended strongly on tumour-size/nodal status (e.g. 13% for T1N0, 34% for T1N4–9, 41% for T2N4–9). The abstract is the only part read. Two implications for us: (i) the empirical hazard after year 5 is **non-vanishing**, which a model with an absorbing "cure" (tumour elimination) cannot reproduce unless its parameter distribution contains slow-escape members; (ii) group-level risk differs by baseline burden-related covariates.

**Breast-cancer relapse dynamics differ by molecular subgroup.** `Rueda2019` (ABS) models locoregional recurrence, distant recurrence, breast-cancer death and death from other causes as stages with competing risks in 3,240 patients (1,980 with molecular data), and reports four late-recurring ER+/HER2− integrative subgroups (about 26% of such tumours) with a "high risk of recurrence (mean 47–62%) up to 20 years after diagnosis", plus a triple-negative subgroup in which cancer "rarely recurs after five years". This is the clinical motivation for a model that can generate slow, heterogeneous escape times. The study's multistate model is richer than the single composite endpoint available in the public cBioPortal export of the same cohort (see `data/dataset_feasibility.md`).

---

## 2. Mathematical models of dormancy and immune control

### 2.1 Wilkie & Hahnfeldt 2013, *Cancer Research* (`Wilkie2013CancerRes`), **FT**

- **Citation:** Wilkie KP, Hahnfeldt P. *Cancer Res* 2013;73(12):3534–3544. doi:10.1158/0008-5472.CAN-12-4590.
- **Objective:** Show how cancer-cell *sensitivity to microenvironmental growth-regulating signals* shapes tumour fate, and explain immune-induced dormancy as a *transient* state.
- **Mechanisms:** Immune predation on cancer cells; immune recruitment stimulated by cancer cells; an environmental carrying capacity that can grow with the tumour (angiogenesis-like feedback).
- **Model (as written in the paper):**
  - Cancer: dC/dt = μ α (1 + Ψ(I,C)) C (1 − (C/K_C)^α)
  - Immune: dI/dt = λ (I + r C)(1 − I/K_I)
  - Predation: Ψ(I,C) = −θ ( I^β / (φ C^β + I^β) + ε log10(1+I) ), with ε = 0.01 *assumed*.
  - Dynamic capacity (from Hahnfeldt et al. 1999, `Hahnfeldt1999`): dK_C/dt = p C − q K_C C^{2/3}.
- **Analysis:** Analytic equilibria and Jacobian stability. Four equilibria: (0,0) unstable; (K_I,0) stable if 1+Ψ(K_I,0) < 0, else saddle; (K_I,K_C) stable if 1+Ψ(K_I,K_C) > 0, else saddle; and an interior point (K_I, C⁻) with 1+Ψ(K_I,C⁻)=0, a **saddle** when both boundary points are stable. θ is a bifurcation parameter with **two transcritical bifurcations** at θ_min and θ_max, bounding the immune-efficacy range that "guarantees the existence of all possible tumor outcomes". Parameters were estimated by MCMC (10 runs × 20,000 iterations) from Tanooka et al. subcutaneous 3-methylcholanthrene fibrosarcoma growth in mice, *under the assumption of no immune predation*. Simulation in MAPLE.
- **Data:** Murine tumour-volume data only. The authors state that immune size and growth-rate values "are chosen to approximately match those of the tumor since methods to measure immune system checkpoint blockades and other size-determining factors remain elusive".
- **Findings:** (i) Dormancy is predicted to be "a transient period of growth that necessarily ends in either tumor elimination or tumor escape". The dormant interior equilibrium is a saddle, so duration is set by the trajectory's proximity to its separatrix. (ii) Four parameter sets fit untreated growth equally well but predict very different outcomes under the same immunotherapy: a direct warning about **non-identifiability from growth curves alone**. (iii) Disrupting the environmental signal can accelerate regrowth ("accelerated repopulation"). (iv) They contrast their result with predator–prey-type models where dormancy is a stable steady state or limit cycle.
- **Dormancy definition used:** ≥30 days within the first 150 days where (1/C)dC/dt < 0.04 and tumour diameter is 1–30 mm. This is an *operational* definition, as ours must be.
- **Limitations (authors):** The model "does not explicitly address" immunoediting or immune polarisation, and immune parameters were not measured. *Inferred by us:* single mouse tumour type; no explicit quiescent compartment; no time-to-event comparison with patients.
- **Relevance to P–D–E:** Closest conceptual precedent for (a) analytic equilibrium/Jacobian analysis of dormancy, (b) transcritical structure at an immunity threshold (our k_c boundary is also a transcritical-type crossing of the tumour-free state, see `parameter_evidence.md` §5), and (c) using operational dormancy criteria. **Difference:** their immune recruitment is logistic and tumour-triggered with I(0)=0 in the treatment simulations, so the tumour-free state is not an equilibrium with immune surveillance; ours has constant baseline replenishment `s`, giving E* = s/μ.
- **Suggested improvements from the paper:** a dynamic carrying capacity; heterogeneity in sensitivity to environmental signals; immunoediting/resistance.

### 2.2 Wilkie & Hahnfeldt 2013, *Interface Focus* (`Wilkie2013Interface`), ABS

- **Citation:** *Interface Focus* 2013;3(4):20130010. doi:10.1098/rsfs.2013.0010.
- **Verified (abstract):** Experimental work in mice suggests that in immune-induced dormancy "the longer a cancer remains dormant in a host, the more resistant the cancer cells become to cytotoxic T-cell-mediated killing". Using a data-driven approach the authors compare *decaying efficacy of immune predation* with *declining immune recruitment* as escape mechanisms, concluding that "decline in recruitment is a stronger determinant of escape than increased resistance to predation". A mechanistic approach adds an immune-resistant subpopulation and analyses *immunoediting*, in which selective killing of immune-sensitive cells makes the tumour more homogeneous and resistant.
- **Not verified:** equations, parameter values, datasets, figures.
- **Relevance:** Supports a future extension in which `k` (killing) or `s`/`ρ` (immune supply) changes over time, and shows escape can be driven by either *immune failure* or *tumour adaptation*. Our constant-parameter baseline cannot distinguish them.

### 2.3 Yáñez Feliú, Rossato, Valleriani & Cipitria 2026 (`YanezFeliu2026`), ABS

- **Citation:** *npj Syst Biol Appl* (online 13 Oct 2026). doi:10.1038/s41540-026-00827-8.
- **Verified (abstract):** An evolutionary-game-theory framework with environmental feedback. Cancer cells produce a conditioning factor that is cleared by the microenvironment and "alters cancer–host competition once it approaches a characteristic concentration". The coupling yields distinct long-term regimes; "prior niche conditioning can also create history dependence, so that the same cancer population can follow different fates depending on whether the microenvironment is already primed"; repeated small waves of arriving cells can cross a threshold cumulatively; and "closeness to either threshold generates large variation in the time to outgrowth".
- **Not verified:** the equations, analysis details, parameter values and any data.
- **Relevance:** Directly motivates the optional M(t) extension in the project plan (§19) and an *a priori* expectation that outgrowth-time variance is large near thresholds. For our baseline it is a reminder that a constant awakening rate `a` is a strong simplification.

### 2.4 Blath, Kraut, Paul & Tóbiás 2025 (`Blath2025`), ABS

- **Citation:** *J Theor Biol* 2025;597:111995. doi:10.1016/j.jtbi.2024.111995.
- **Verified (abstract):** An idealised stochastic agent-based model of cells "switching into and out of a dormant state both spontaneously as well as in response to treatment", with multi-drug protocols. Analysis uses simulations of the many-particle limit. Result: "even a small initial population (of explicitly quantifiable size) of dormant cells can lead to therapy failure under classical single-drug treatments that would successfully eradicate the tumour in the absence of dormancy", with outcomes depending on the dormancy mechanism; several multi-drug regimes are quantified.
- **Not verified:** the exact rate laws, switching mechanisms and numerical values.
- **Relevance:** Supports the P↔D switching structure and the claim that dormancy changes therapy outcomes, and motivates a stochastic extension. It concerns *short-term* dormancy and *therapy response*, whereas our model targets *post-treatment* long-latency relapse.

### 2.5 Supporting modelling literature (ABS)

| Paper | What the abstract establishes | Relevance to P–D–E |
|---|---|---|
| **Kuznetsov et al. 1994**, *Bull Math Biol* 56:295–321 (`Kuznetsov1994`) | Model of the cytotoxic T-lymphocyte response to an immunogenic tumour; reproduces immunostimulation, "sneaking through", and a tumour "dormant state"; parameters estimated against BCL1 B-lymphoma kinetics in mouse spleen; local and global bifurcations computed; predicts for many parameters a recurrent profile "with a 3- to 4-month cycle". | The classic precedent for a bifurcation-analysed tumour–effector model with parameters estimated from data. **Equation-level comparison with our immune equation is not done**, since only the abstract was read. Our saturating tumour-stimulated recruitment term resembles models of this family, to be confirmed against the full text. |
| **Page & Uhr 2005**, *Leuk Lymphoma* 46:313–327 (`Page2005`) | Preliminary models of tumour–antibody interaction in BCL1 lymphoma; immune-induced "cell cycle arrest and apoptosis"; analyses parameters controlling "the size of the dormant tumors, and the fraction of proliferating cells". | The nearest precedent we found for a *proliferating vs arrested* split inside an immune-control model. |
| **Wilkie 2013**, *Adv Exp Med Biol* 734:201–234 (`Wilkie2013Review`) | Review of predator–prey-type cancer–immune dormancy models, including augmentations such as "the cancer cell transition to a quiescent state" and "time delay in T-cell activation". | Establishes that a quiescent compartment in an immune model is a recognised augmentation, so our novelty is not the compartment itself. |
| **de Pillis et al. 2005**, *Cancer Res* 65:7950–7958 (`deP2005`) | NK and CD8⁺ T-cell tumour-surveillance ODE model fitted/validated against published mouse and human data (chromium-release lysis, tumour growth) with sensitivity analysis. | Provides a *data-validated* route for immune-parameter ranges. Our `k`, `s`, `μ`, `ρ`, `h`, `η` are unvalidated. |
| **Enderling et al. 2009**, *Cancer Res* 69:8814–8821 (`Enderling2009`) | Cell-based stem/non-stem model; dormancy as balanced proliferation and death; increased cell death produces short-term reduction but "paradoxical accelerated long-term growth" via liberation of cancer stem cells. | Warns that raising a death term can backfire; also indicates that spatial/space constraints matter. Both are outside our well-mixed ODE. |
| **Mehdizadeh et al. 2021**, *Cancers* 13:245 (`Mehdizadeh2021`) | Deterministic model with two subpopulations (high type-I-IFN signalling / immunogenic / slower-growing vs the opposite) in murine 4T1 triple-negative breast cancer; reproduces dormancy, elimination and escape; dormant-phenotype cells can suppress growth of an aggressive counterpart. | **Counter-evidence for our assumption 6**: here the dormant-like population is *more* immunogenic, not immune-invisible. |
| **Saudemont et al. 2005**, *Blood* 105:2428–2435 (`Saudemont2005`) | In a mouse AML model, long-term persistent (dormant) leukaemic cells resist CTL lysis (overexpressing B7-H1) but are killed by NK cells. | Supports a *cell-state-dependent immune susceptibility*. Dormant cells may escape one effector arm but not another. |
| **Baffoe & Bouchnita 2026**, *Front Immunol* 17:1861769 (`Baffoe2026`) | Dormancy as a continuous ERK/p38 phenotype axis rather than a binary state; regimes of elimination, equilibrium and escape; reduced model with critical thresholds in mean phenotype. | Direct critique of the binary P/D simplification and a route to a continuous-phenotype extension. |

---

## 3. Clinical and data literature

| Paper | What the abstract establishes | Use |
|---|---|---|
| `Rueda2019` | See §1. Cohort = METABRIC-type long-term-follow-up breast cohort; the study's `meta_study.txt` entry in cBioPortal lists Rueda et al. 2019 among the study's citations. | Clinical motivation; model-based competing-risk framework. |
| `Curtis2012`, `Pereira2016` | METABRIC discovery/validation cohorts (997 + 995 tumours) with "long-term clinical follow-up"; targeted sequencing of 2,433 primary tumours with copy-number, expression and follow-up data. | Provenance of the dataset audited in Phase 1. |
| `Pan2017` | See §1. | Late-hazard target for group-level comparison. |
| `Tolaney2021` (STEEP 2.0) | Standardised endpoint definitions for adjuvant breast-cancer trials; simulations show that including second non-breast primary cancers in invasive-disease-free survival can reduce power and "mask differences in recurrence rates". | Justifies insisting that "RFS", "DFS" and "OS" are not interchangeable. The abstract does not itself list the endpoint definitions. |

---

## 4. Comparison of models

| Feature | **Our P–D–E** | Wilkie & Hahnfeldt 2013 (CR) | Kuznetsov 1994 | Page & Uhr 2005 | Blath 2025 | Yáñez Feliú 2026 | Baffoe 2026 |
|---|---|---|---|---|---|---|---|
| Quiescent/dormant cell state | **Yes (binary D)** | No | No (abstract) | Arrested vs proliferating (abstract) | Yes (stochastic) | Not stated in abstract | Continuous ERK/p38 axis |
| Immune compartment | Yes (E, with surveillance source `s`) | Yes (I, logistic recruitment) | Yes (effector cells) | Antibody/immune response | No (treatment) | "Host" competition | Yes |
| Environment feedback | No | **Yes** (dynamic K_C) | No | No | No | **Yes** (conditioning factor) | Stress/growth factors |
| Analysis | Equilibria, Jacobian, simulation | Equilibria, Jacobian, bifurcation, MCMC | Local/global bifurcation, fit | Parameter analysis | Many-particle simulation | Regime analysis (abstract) | Reduced-model thresholds |
| Stochastic/uncertainty | Planned (Monte Carlo) | MCMC parameter sets | No | No | **Yes** | Time-to-outgrowth variation | No |
| Fitted to data | **No (illustrative)** | Mouse growth data (immune parameters assumed) | Mouse lymphoma | Not stated | Not stated | Not stated | Not stated |
| Clinical time-to-event link | Planned | No | No | No | No | No | No |

---

## 5. Research gaps (relative to the verified literature)

1. **Few models link dormancy–immune ODEs to time-to-event distributions.** Among the sources read, none (by abstract) maps a dormancy ODE onto a clinical recurrence-free curve. This is where the project is positioned, with the caveat that the search was bounded and is not a systematic review.
2. **Dormancy is often "a stable state" in classic models and "transient" in Wilkie & Hahnfeldt.** The consequences for *relapse-time distributions* (as opposed to qualitative regimes) are not examined in the abstracts we read.
3. **Identifiability of the dormant compartment from total-burden or clinical data** is acknowledged implicitly (Wilkie & Hahnfeldt's four equally fitting parameter sets) but not treated systematically.
4. **Competing risks and the time origin** are central in clinical data (`Rueda2019`, `Tolaney2021`) yet are rarely part of mechanistic-model validation.
5. **Binary versus continuous dormancy** (`Baffoe2026`) is unresolved for population-level models.
6. **Immune susceptibility of dormant cells is contested** (`Mehdizadeh2021`, `Saudemont2005` vs the common assumption of immune evasion).

## 6. Caveats of this review

- Of the five required papers, only one was read in full text (four were abstract-level). Abstract-level claims should be re-checked against full texts before appearing in the final report.
- Equation-level comparison with Kuznetsov 1994, d'Onofrio 2005 and Interface Focus 2013 is deferred.
- The search for additional papers used Europe PMC title/topic queries and the reference list of Wilkie & Hahnfeldt 2013; it is not exhaustive and was limited to sources whose metadata could be verified in this session.
