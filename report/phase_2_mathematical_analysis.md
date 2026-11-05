# Phase 2: Formal Bifurcation Analysis, Parameter Normalisation and Relapse-Time Approximation


Reproduce: `.venv/bin/python experiments/bifurcation_analysis.py`, `experiments/relapse_time_approximation.py`, `experiments/normalized_model_comparison.py`, then `.venv/bin/python -m pytest` (71 tests). Tables: `results/tables/`, figures: `results/figures/`.

---

## 0. Model, notation, baseline

State (P, D, E): proliferating tumour burden, dormant tumour burden (both fractions of K), effective immune activity. Time in months.

```
dP/dt = r P (1 − N/K) − q P + a D − k E P            N = P + D
dD/dt = q P − (a + δ_D) D
dE/dt = s + ρ P/(h+P) E − μ E − η N E
```

Baseline: r = 0.12, q = 0.04, a = 0.005, δ_D = 0.001, k = 0.1, s = 0.02, μ = 0.05, ρ = 0.04, h = 0.02, η = 0.08, K = 1.

Derived constants used throughout:

| symbol | definition | baseline |
|---|---|---|
| c | q/(a+δ_D) | 6.667 |
| φ | 1 + c = N*/P* on every nonzero equilibrium | 7.667 |
| r_eff | r − qδ_D/(a+δ_D) | 0.11333 |
| E_f | s/μ (tumour-free immune level) | 0.4 |
| k_c | r_eff μ / s | 0.28333 |
| β | k s / μ (effective killing, §5) | 0.04 at k = 0.1 |

---

## 1. Equilibria: exact reduction and exhaustive enumeration

**Tumour-free equilibrium.** P = D = 0 gives E = s/μ. Conversely, dD = 0 with P = 0 forces D = 0, so there is no equilibrium with P = 0, D > 0.

**Nonzero equilibria.** From dD = 0: D = cP, N = φP. Substituting into dP = 0 and dividing by P > 0, using −q + a c = −qδ_D/(a+δ_D):

```
r_eff − (r φ / K) P = k E ,        E = s / B(P) ,   B(P) = μ + η φ P − ρ P/(h+P)        (1)
```

Because k appears linearly, every nonzero equilibrium lies on the **explicit graph**

```
k = k(P) = A(P) B(P) / s ,      A(P) = r_eff − (r φ/K) P ,      D = cP,  E = s/B(P)       (2)
```

(`equilibria.branch_k`). Consequences, all exact:

* **Admissible range.** E > 0 ⇒ k(P) > 0 ⇒ A(P) > 0, so P < r_eff K/(r φ): nonzero equilibria are bounded, and none exist if r_eff ≤ 0.
* **At most three nonzero equilibria.** Multiplying (1) by (h+P) gives the cubic `A(P)·[(μ+ηφP)(h+P) − ρP] − k s (h+P) = 0`. Hence the number of nonzero equilibria is ≤ 3 for every parameter set; the Phase-1 grid search and the multistart check below can be compared against this hard bound (`equilibria.cubic_coefficients`).
* **Dormant fraction is an identity.** D*/N* = c/(1+c) = q/(q+a+δ_D) = 0.8696 at *every* nonzero equilibrium for *every* k (Fig. D1). It is a property of the dormancy switching rates, not an outcome of immune control.
* **Constant term.** At P = 0 the cubic's constant term is h(r_eff μ − k s), which vanishes exactly at k = k_c: the nonzero branch passes through the tumour-free state at k_c (proof that the two branches intersect).

**Numerical verification** (`enumerate_equilibria`: cubic roots → rejection of complex / P ≤ 0 / E ≤ 0 → Newton polish on the full 3-D system → residual check). Independent `multistart_equilibria` (300 random log-uniform starts per k, no knowledge of the cubic) agrees in the number and location of equilibria at all 8 tested k (`equilibrium_enumeration_crosscheck.csv`); residuals ≤ 6×10⁻¹⁸. Rejected roots are returned with a reason (e.g. the P < 0 root for k > k_c).

---

## 2. Continuation (`src/dormancy/continuation.py`)

Pseudo-arclength predictor–corrector on F(x, k) = 0, x = (P, D, E), scaled variables, analytic Jacobian and ∂F/∂k = (−EP, 0, 0).

* **Tangent** = null vector of the 3×4 matrix [F_x F_k] (SVD); at a branch point the null space is 2-D and the direction closest to the previous tangent is chosen, so the continuation passes straight through intersections.
* **Corrector** Newton on [F = 0; t·(u − u_pred) = 0]. Convergence is judged on the Newton *step* (≤ 10⁻¹⁴), not on the residual: near the transcritical point P ~ 10⁻⁸, F_P ∝ P is tiny even for points that are not on the branch (an early version with a residual-only test mislocated k_c by 9×10⁻⁵).
* **Adaptive step.** ×1.6 after ≤ 3 Newton iterations with tangent angle < ≈ 2.5°; ×0.6 after ≥ 6; halved on failure or on tangent angle > ≈ 11°; additionally reduced when linear extrapolation of a monitored function predicts a sign change within the next step.
* **Monitored functions** (sign changes located by bisection on the corrector curve): det J (zero real eigenvalue: fold or branch point), a₁a₂ − a₃ (Hopf, with a₂ > 0 required for a complex pair; with a₂ < 0 recorded as a "neutral saddle" and *not* a bifurcation), and P (leaving the orthant). At a branch point the corrector is singular (every k solves F_P = 0 on P = 0) so the bisection stops at |P| < 2×10⁻⁷ and extrapolates quadratically in P to P = 0 (error O(P³)).
* **Every accepted point** stores the full Jacobian eigenvalues, a stability label from **all** eigenvalues (stable node/focus, saddle n_u = 1/2, unstable node/focus, marginal), max |F|, and a `physical` flag (P, D, E ≥ 0). Nonphysical points (the P < 0 continuation) are tracked but never reported as equilibria.
* **Multi-seed tracing** (`trace_all_branches`): exact admissible equilibria at 25 values of k across [0, 3k_c] seed forward/backward continuations; a seed already on a traced branch is skipped, so disconnected branches would be found.

**Results at baseline** (`bifurcation_branches_baseline.csv`, `bifurcation_events_baseline.csv`): 217 points; max |F| = 8.7×10⁻¹⁸; max |k − k(P)| = 2.5×10⁻¹⁶ against the exact curve (2); one nonzero branch; the only bifurcation is the branch point at k = 0.283333333333 (analytic 0.283333333333). The "neutral saddle" events at k = 0.0839, 0.1850 on the tumour-free branch are points where two Jacobian eigenvalues sum to zero (−μ + λ₊ = 0 etc.), real, not a Hopf. **Convergence in the step cap** (`continuation_convergence.csv`): ds_max = 0.2 … 0.01 gives identical k(P) accuracy (≈ 2×10⁻¹⁶) and branch-point error ≤ 7×10⁻¹⁰.

---

## 3. The transcritical bifurcation at k = k_c

### 3.1 Definitions used

A **transcritical bifurcation** at (x₀, k₀) of ẋ = F(x, k) requires, with a *simple* zero eigenvalue of J = F_x(x₀,k₀) with right/left null vectors w, ℓ (ℓ·w = 1) and the other eigenvalues off the imaginary axis (Sotomayor / Kuznetsov):

1. F_k(x₀,k₀)·ℓ = 0 (the equilibrium persists for k ≠ k₀: here x₀ is an equilibrium for all k);
2. **transversality** α = ℓ·[F_{xk} w] ≠ 0 (the eigenvalue crosses zero with nonzero speed);
3. **nondegeneracy** β = ½ ℓ·D²F(w,w) ≠ 0.

Then, on the 1-D centre manifold, u̇ = α(k−k₀) u + β u² + O(u³, (k−k₀)u²): two branches u = 0 and u* = −α(k−k₀)/β intersect at k₀ and **exchange stability**, since the eigenvalues are λ(0) = α(k−k₀) and λ(u*) = −α(k−k₀).

### 3.2 Derivation for the model

At x₀ = (0, 0, E_f) the Jacobian is block lower-triangular, because the E-column vanishes in the P and D rows (the only entries there, −kP, are zero at P = 0), so the (P, D) block decouples:

```
J = [[ r − q − kE_f,  a,        0 ],
     [ q,            −(a+δ_D),  0 ],
     [ (ρ/h − η)E_f, −ηE_f,    −μ ]]
```

so eigenvalues = {−μ} ∪ eig(M), M = [[r − q − β, a],[q, −(a+δ_D)]], β = kE_f. det M = −(a+δ_D)(r − q − β) − aq = 0 ⇔ β = r_eff, i.e. **k_c = r_eff μ/s = 0.28333**, and tr M = −qa/(a+δ_D) − (a+δ_D) < 0 there, so the other two eigenvalues are −μ = −0.05 and −0.03933 (simple zero eigenvalue, others in Re < 0; numerically verified).

Null vectors at k_c (M = [[−qa/(a+δ_D), a],[q, −(a+δ_D)]]):

```
w = (1, c, w_E),   w_E = E_f [ρ/h − η(1+c)] / μ = 11.093      (third row of J w = 0, w_P = 1)
ℓ = (ℓ_P, ℓ_D, 0) = ( 1/(1 + qa/(a+δ_D)²),  ℓ_P·a/(a+δ_D),  0 ),  ℓ_P = 0.15254, ℓ_D = 0.12712
```

Conditions (computed with the analytic Jacobian, analytic F_k and analytic D²F, `bifurcation.transcritical_report`):

| condition | value |
|---|---|
| ℓ·F_k at (x₀,k_c) | 0 (F_k = (−EP,0,0) = 0 at P = 0) |
| α = ℓ·(J_k w) = −E_f ℓ_P w_P | **−0.061017** ≠ 0 |
| β = ½ ℓ·D²F(w,w) | **−0.61980** ≠ 0 |
| dλ/dk (finite-difference of the Jacobian eigenvalue) | −0.0610169 (= α to 7 digits) |
| simple zero eigenvalue, others Re < 0 | yes |

`D²F` is checked against symbolic differentiation (sympy) in the tests.

### 3.3 Branch behaviour near k_c (derived and compared)

From (2): k′(0) = [−rφμ/K + r_eff(ηφ − ρ/h)]/s = **−10.158** (baseline). Hence near k_c the nonzero branch is

```
P*(k) = (k − k_c)/k′(0) + O((k−k_c)²) = −0.098447 (k − k_c) + …           (3)
```

The centre-manifold prediction −α w_P/β = **−0.0984467294** equals 1/k′(0) = −0.0984467294 to 16 digits; algebraically β = −α w_P k′(0), so **nondegeneracy (β ≠ 0) ⇔ k′(0) ≠ 0** (an exact identity checked numerically). Numerical branch vs (3) (`transcritical_asymptotics.csv`): relative error 17.9% at ε = 0.1, 1.7×10⁻³ at ε = 10⁻³, 1.7×10⁻⁴ at ε = 10⁻⁴ (first-order convergence, as expected).

### 3.4 Stability exchange

`transcritical_asymptotics.csv`, 12 values of k = k_c(1 ± ε), ε ∈ [10⁻⁴, 10⁻¹]: for **every** sampled k exactly one of the two branches is stable.

* k < k_c: tumour-free is a saddle (max Re λ = +α(k−k_c) > 0), the nonzero equilibrium (P > 0, admissible) is stable.
* k > k_c: tumour-free is stable; the algebraic continuation of the nonzero root has **P < 0** and is unstable (max Re λ ≈ −α(k−k_c) > 0).
* The ratios λ_TF/[α(k−k_c)] and λ_nz/[−α(k−k_c)] tend to 1 (1.0002 and 0.9997 at ε = 10⁻⁴).

### 3.5 What may and may not be claimed

* **Established (mathematically):** at k = k_c the system has a simple zero eigenvalue, the tumour-free branch and the nonzero branch intersect, the eigenvalue crosses zero transversally (α ≠ 0), the quadratic coefficient is nonzero (β ≠ 0), and stability is exchanged. These are the Sotomayor conditions, so **the point is a transcritical bifurcation of the extended system (P ∈ (−h, ∞))**, verified analytically (§3.2), by continuation (k_c located to ≤ 10⁻⁹) and by root enumeration.
* **Needs the qualification:** the exchanging branch is physical (P > 0) only on one side. On the nonzero branch D = cP, so it passes through the origin of the (P, D) plane and P changes sign with k − k_c whenever k′(0) ≠ 0; an exchange with admissible branches on *both* sides is therefore impossible. In the orthant the picture is a **forward ("supercritical") bifurcation of the nonzero state out of the tumour-free state at k_c**: below k_c a stable low-burden state exists, above k_c only the tumour-free state is admissible. The label "transcritical" is therefore the right local normal form, but its biological reading is "branching-off at k_c", not "two admissible states swapping".
* **Direction is not universal** (§4.2): the above holds for k′(0) < 0 (baseline). For k′(0) > 0 the same point is still transcritical but backward (subcritical): the nonzero branch lives at k > k_c.

---

## 4. Additional bifurcations

### 4.1 Method

Because (2) is a graph over P, a *complete* scan for every k at once is possible: sample P over the admissible range (A > 0, B > 0), compute the Jacobian, its characteristic polynomial λ³ + a₁λ² + a₂λ + a₃, and evaluate (i) dk/dP (sign changes = saddle-nodes; exact roots of a cubic, `equilibria.fold_P`), (ii) Routh–Hurwitz stability, (iii) H = a₁a₂ − a₃ (a sign change with a₂ > 0 is a purely imaginary pair, i.e. a Hopf candidate). `bifurcation.scan_branch`.

Scans (`parameter_scan_summary.csv`): 8,000 Latin-hypercube samples each, log-uniform.
* **Scan A:** a ∈ [10⁻⁴,1], q ∈ [10⁻³,1], ρ ∈ [10⁻³,2], η ∈ [10⁻³,5]; others at baseline. 7,278 valid (r_eff > 0).
* **Scan B:** A plus r ∈ [0.02,0.5], δ_D ∈ [10⁻⁴,0.1], μ ∈ [0.01,0.5], h ∈ [10⁻³,0.5]. 6,777 valid.
* **k** needs no sampling: the branch graph (2) spans every k ≥ 0 at once.

| finding | A | B |
|---|---|---|
| forward at k_c / backward at k_c | 5,602 / 1,676 | 5,382 / 1,395 |
| ≥ 1 saddle-node | 3,078 | 2,638 |
| two saddle-nodes | 850 | 773 |
| tumour-free **and** a nonzero state both stable at some k | 2,580 | 2,212 |
| ≥ 2 stable *nonzero* equilibria at one k (window always below k_c) | 1,349 | 1,206 |
| **Hopf candidate** | **0** | **0** |
| "neutral saddle" (a₁a₂ = a₃, a₂ < 0; not a bifurcation) | 672 | 481 |

### 4.2 Saddle-node bifurcations and bistability (verified)

* **Backward transcritical bifurcation:** k′(0) > 0 ⇔ η > η* = [ρ/h + rφμ/(K r_eff)]/φ = **0.3138** (baseline others), or ρ < ρ* = h(ηφ − rφμ/(K r_eff)) = **0.00415**. At η = η*, β = 0 (degenerate; a codimension-2 point). Example η = 0.6 (`bifurcation_*_eta0p6.csv`, Fig. B2): continuation locates a fold at **k_SN = 0.860834** (P = 0.058486, N = 0.4484), matching the exact root of the cubic to 10⁻⁸; between k_c = 0.2833 and k_SN the tumour-free state (stable) coexists with a stable high-burden state and an unstable saddle between them.
* **Simulation check** (`bistability_simulations.csv`, η = 0.6): at k = 0.572 and 0.775, an initial condition with small residual disease (P₀ = 10⁻⁶, D₀ = 10⁻⁴, E₀ = E_f) converges to the tumour-free state; an initial condition with large burden (0.3, 0.3, E = 0.1) converges to the nonzero equilibrium (P* = 0.0436 and 0.0524 respectively). At k = 0.5 k_c both converge to the same nonzero state (one stable equilibrium).
* **Fold pairs with a forward transcritical point (two distinct kinds of coexistence).** For smaller η (or larger ρ) the nonzero branch is S-shaped *without* being backward at k_c (k′(0) < 0). Where the folds lie relative to k_c decides what coexists (verified with `enumerate_equilibria` and eigenvalue labels):
  * **η = 0.12** (baseline others): folds at k = 0.22933 and 0.23025, **both below k_c**. In (0.22933, 0.23025) there are three nonzero equilibria (stable, saddle, stable) while the tumour-free state is a saddle: **bistability between two *nonzero* states**, no tumour-free coexistence.
  * **η = 0.164**: folds at k = 0.2561 (below k_c) and 0.2840 (above k_c). On (0.2561, 0.2833) two stable nonzero states (a focus at P = 0.0026 and a node at P = 0.065, at k = 0.27) coexist with a saddle (P = 0.027) and an *unstable* tumour-free state. Only on (0.2833, 0.2840) do the stable tumour-free state and one stable nonzero state coexist (at k = 0.28367: saddle P = 0.0447, stable node P = 0.0506).
  * First fold appears at η ≈ 0.116; tumour-free + nonzero bistability onsets at η ≈ 0.164 (and ρ ≲ 0.005 at baseline η); at η = 0.4 (backward) the tumour-free + nonzero window is (0.283, ≈ 0.59). The baseline (η = 0.08, ρ = 0.04) has neither kind of coexistence, about a factor 2 below the onset in η (Fig. F1, F2).
* **Two stable nonzero equilibria at one k (nonzero–nonzero multistability)** occurs in about 18% of the sampled sets (1,349 + 1,206 of 14,055; in *every* one the coexistence window lies below k_c, where the tumour-free state is a saddle, `n_multi_nonzero_with_unstable_TF_only`). Example (Scan A, row `multi_stable_nonzero`): a = 0.306, q = 1.6×10⁻³, ρ = 0.240, η = 1.637, k = 0.15 has three nonzero roots, P = 0.0027 (stable focus), 0.111 (saddle) and 0.977 (stable node), plus an unstable tumour-free state. 30 random initial conditions simulated for 3×10⁵ months end only at P = 0.0027 or 0.977. Interpretation: a low-burden immune-controlled state versus a near-carrying-capacity state; strong immune stimulation (ρ/μ = 4.8 here) is the mechanism.

### 4.3 Hopf bifurcation and oscillations: not found

* Along the exact branch graph in both scans: **no sign change of a₁a₂ − a₃ with a₂ > 0** among 14,055 valid parameter sets. Stable *foci* (complex eigenvalues with Re < 0, i.e. damped oscillation) occur on the stable branches in 354 of the parameter sets examined for this (every 4th draw of each scan, first 50 stable-branch points each); none crosses the imaginary axis.
* In the continuation of the *nonphysical* P < 0 continuation at η = 0.6, one Hopf candidate is detected (k = 0.2396, P = −0.0046, a₂ > 0, eigenvalues ±0.0061 i): it is on a part of the branch with P < 0, **not an admissible state**, and is not reported as a model bifurcation.
* **Direct simulation search** (`oscillation_search*.csv`): 263 random (parameter set, k ∈ [0, 2k_c], initial condition) draws integrated to 1.5×10⁴ months. 15 were flagged by a crude amplitude criterion; re-integrating each to 3×10⁶ months (Radau) shows all converge to an equilibrium with residual ≤ 10⁻¹⁷ (extrema counted at relative amplitude 10⁻¹⁰ are round-off). **0 persistent oscillations.**
* **Conclusion:** "no Hopf bifurcation and no sustained oscillation in the scanned ranges". This is *not* a proof of absence: ranges are bounded and sampled, and no global argument (Dulac, Lyapunov) excludes limit cycles. Global boundedness of trajectories is also not proven (for ρ > μ the immune equation is not obviously dissipative; none of the sampled trajectories diverged).

---

## 5. M1: normalising the immune-effector scale

**Redundancy.** The gauge map (E, s, k) → (cE, cs, k/c) leaves P and D dynamics unchanged. Only β = ks/μ enters the tumour equations at low burden; E has no natural unit.

**Transformation.** Let e = E/(s/μ) = μE/s (tumour-free level e = 1) and β = ks/μ. Then E = (s/μ)e, kE = βe, and de/dt = (μ/s)dE/dt:

```
dP/dt = r P (1 − N/K) − q P + a D − β e P
dD/dt = q P − (a + δ_D) D
de/dt = μ (1 − e) + ρ P/(h+P) e − η N e
```

(derived: (μ/s)(s + ρ g E − μE − ηNE) = μ + (ρ g − μ − ηN) e = μ(1−e) + ρ g e − ηNe). Parameters k and s disappear (10 parameters: K, r, q, a, δ_D, β, μ, ρ, h, η). Units: β [1/month], e dimensionless. Consequences, all exact:

* β_c = r_eff = r − qδ_D/(a+δ_D) = 0.11333 (k_c = β_c μ/s); the ratio β/β_c = k/k_c.
* The Jacobians are conjugate by S = diag(1,1,μ/s). Hence eigenvalues, bifurcation points and branch burdens are identical.
* In the original model, E = (s/μ) e; everything about E is only meaningful once a normalisation is declared.

**Implementation (backward compatible).** `src/dormancy/normalized.py`: `NormalizedParams` (separate dataclass), an independent `rhs_normalized`, `jacobian_normalized`, `simulate_normalized`, the maps `to_params(s)`, `from_params`, state maps, and `gauge_transform(p, c)`; `configs/normalized.yaml` (β = 0.04 ⇔ baseline k = 0.1, s = 0.02, μ = 0.05). `Params.beta`, `Params.r_eff`, `Params.phi` were added as read-only properties. **No existing parameter, config, result file or default changed.** `baseline_scenarios.csv` regenerates byte-identical.

**Verification** (`normalized_model_equivalence.csv`, Fig. H, `tests/test_normalized.py`): the original model at s ∈ {0.002, 0.02, 0.2, 2} (k rescaled so that β is fixed; a 10³ range) at k/k_c ∈ {0.25, 0.5, 0.75}: relapse-time difference to the normalised model ≤ 8.4×10⁻⁷ months (solver tolerance), max |ΔP, ΔD| at the event ≤ 7.6×10⁻¹² and E = (s/μ)e to 10⁻⁸ relative. RHS agreement to 10⁻¹³; a negative control (β changed by 7%) changes T by > 1 month, so the test is not vacuous. The two shipped configs agree on β.

**What M1 does *not* resolve.** It fixes the unit of E, not identifiability of ρ, h, η from burden data, and it does not remove the dependence on μ (turnover and source are tied in e: source = μ).

---

## 6. Relapse-time approximation

Setting: x₀ = (P₀, D₀, E_f), N₀ = 1.01×10⁻⁴ (P₀ = 10⁻⁶, D₀ = 10⁻⁴), threshold N_det, k < k_c. Numerical truth = first upward crossing of N = N_det by DOP853 (rtol 10⁻⁹), horizon 10⁵ months.

### 6.1 Estimators

Let λ₊ be the positive real eigenvalue of J at the tumour-free state, w its right eigenvector (w_P = 1), ℓ the left eigenvector (ℓ·w = 1), n₁ = w_P + w_D, u₀ = ℓ·(x₀ − x_TF).

* **A0 (Phase 1):** T = ln(N_det/N₀)/λ₊. Assumes the whole initial burden sits on the growing mode.
* **A1 (eigen-projected):** the linear solution is x − x_TF = Σ cᵢ vᵢ e^{λᵢ t}; the growing mode's amplitude in N is n₁ u₀, hence T = ln(N_det/(n₁u₀))/λ₊ (`T_projected`). **A1b** solves N_lin(t) = N_det exactly (both modes).
* **A2 (second-order unstable manifold).** The 1-D unstable manifold x = x_TF + u w + u² h₂, with u̇ = λ₊ u + γ u²:
  γ = ½ ℓ·D²F(w,w), (2λ₊I − J) h₂ = ½D²F(w,w) − γ w, N(u) = n₁u + n₂u², n₂ = h₂,P + h₂,D.
  u̇ = λu + γu² is Bernoulli; with v = 1/u: v(t) = (v₀ + γ/λ)e^{−λt} − γ/λ, so
  **T = (1/λ₊) ln[(1/u₀ + γ/λ₊)/(1/u_det + γ/λ₊)]**, u_det the smallest positive root of N(u) = N_det. If the bracket is non-positive (u saturates at u* = −λ/γ before reaching u_det), the estimator returns ∞ ("no relapse predicted"). The immune component is slaved correctly because x includes the E-direction of w and h₂.
* **A3 (matched; SEMI-analytic, not a closed form).** T = T_A2(N_match) + τ, with N_match = 10⁻³ ≪ λ/|γ| and the **tail τ obtained by numerical integration** from the second-order-manifold point to N_det. The closed-form part depends on the initial condition (through u₀); τ depends only on the parameters. At baseline, k = 0.5 k_c, N_det = 0.1, the closed-form part is 58.4 / 77.7 / 123.5 months for initial dormant fractions 0 / 0.5 / 0.99, while **τ = 252.046 months for all three** (tail share of T: 81% / 76% / 67%; column `A3_tail_fraction_of_T`). The agreement of A3 with a numerical first-crossing time at N_det = 0.1 is therefore *largely numerics agreeing with numerics*. Its value is structural: **T = T_lin(x₀) + τ(parameters)**, so initial-condition dependence enters only through ln u₀ (Fig. G3), and A3 shows that the second-order closed form is accurate up to N_match, not that a closed form exists at N_det = 0.1.

### 6.2 Why the original estimate failed (baseline, k = 0.5 k_c, N_det = 0.1)

λ₊ = 0.02904/month, other eigenvalues −0.01171 and −0.05. Numerical T = 375.5 months, A0 = 237.5 (gap 138.0 months):

| contribution | months | evidence |
|---|---|---|
| Initial condition is 99% dormant: only part of the burden is on the growing mode (A1 − A0) | **+44.0** | A1 = 281.5; ℓ = (0.860, 0.123, 0) |
| Logistic saturation and other no-feedback nonlinearity (numerical at ρ = η = 0 − A1) | **+11.5** | ρ = η = 0: numerical 293.0, A1 281.6, A2 293.4 |
| Immune feedback during growth (baseline numerical − numerical at ρ = η = 0) | **+82.5** | E rises to 1.40 E_f at detection (stimulation ρ > suppression η), slowing growth |
| **Total gap** | **138.0** | 375.5 − 237.5 |

(The split is order-dependent; the nonlinear terms do not add exactly. The robust statement is that the projection explains about a third of the gap and the immune response about 60%.)

Variant check (`relapse_time_attribution.csv`): ρ = 0 only: numerical 287.2, A2 287.2; η = 0 only: numerical 459.9 (immune boost uncancelled, E reaches 1.59 E_f). So the largest omitted mechanism is **the second-order immune response**, followed by the initial-condition projection; logistic saturation is minor (≈ 11 months). Quadratic-closure breakdown: P_det/h = N_det/(φ h) = 0.65 at N_det = 0.1, so ρP/(h+P) is outside its Taylor radius, which is why A2 returns ∞ at N_det = 0.1 but succeeds at 10⁻² (P_det/h = 0.065).

### 6.3 Results (`relapse_time_approximation.csv`, `relapse_time_accuracy_summary.csv`, Fig. G)

At k = 0.5 k_c (months; numerical / A0 / A1 / A2 / A3):

| N_det | D₀ fraction 0 | 0.5 | 0.99 |
|---|---|---|---|
| 10⁻³ | 59.8 / 78.9 / 57.9 / 58.4 / 58.4 | 77.8 / 78.9 / 77.2 / 77.7 / 77.7 | 122.9 / 78.9 / 123.0 / 123.5 / 123.5 |
| 10⁻² | 143.5 / 158.2 / 137.2 / 143.8 / 143.4 | 162.7 / 158.2 / 156.5 / 163.1 / 162.7 | 208.4 / 158.2 / 202.2 / 208.9 / 208.5 |
| 0.1 | 310.5 / 237.5 / 216.5 / ∞ / 310.4 | 329.7 / 237.5 / 235.8 / ∞ / 329.8 | 375.5 / 237.5 / 281.5 / ∞ / 375.6 |

Accuracy over k/k_c ∈ [0.05, 0.9 k_det/k_c], three compositions (relative error of T):

| N_det | A0 median / max | A1 | A2 | A3 (semi-analytic; tail share of T: 0% / median 56% / median 76%) |
|---|---|---|---|---|
| 10⁻³ | 20% / 49% | 1.4% / 3.2% | 0.6% / 2.4% | 0.6% / 2.4% |
| 10⁻² | 8.7% / 33% | 2.9% / 13% | 0.19% / 3.4% | 0.054% / 0.14% |
| 0.1 | 24% / 44% | 10% / 42% | 1.7% / 11% (9 of 21 predict no event) | 0.007% / 0.036% |

**Honest summary at N_det = 0.1:** the best *closed-form* estimate is A1 (10% median, 42% max error); A2 is a good closed form only for N_det ≲ 10⁻² (≤ 3.4%); the 0.007% of A3 includes a numerically integrated tail covering 76% (median) of T.

The Phase-1 estimate A0 is also **independent of the initial composition**: at N_det = 10⁻³ it overestimates the all-proliferating start by 32% and underestimates the 99%-dormant start by 36%; A1 removes this because it projects x₀ onto the growing mode.

### 6.4 The validity boundary is k_det, not k_c

For k in (k_det, k_c) the stable nonzero equilibrium has N* < N_det. k_det solves N*(k) = N_det, i.e. k_det = k(P = N_det/φ) from (2): **k_det/k_c = 0.9954, 0.9568, 0.7548** for N_det = 10⁻³, 10⁻², 0.1.

*Is "no event beyond k_det" exact?* That needs the approach to N* to be monotone from below: the nonzero equilibrium is a stable *focus* over much of the range (Fig. C1, weak damping), so overshoot above N* is conceivable. Fine scan (`relapse_time_kdet_boundary.csv`): k/k_det ∈ {0.98, … , 1.06} (12 values) × D₀ fractions {0, 0.5, 0.99} × N_det ∈ {10⁻³, 10⁻², 0.1}, horizon 10⁵ months. **0 of 54 runs with k > 1.001 k_det produced an event**; the maximum N along these trajectories stays below N* (e.g. 0.088 at k = 1.02 k_det for N_det = 0.1); the event time grows as k → k_det⁻ (N_det = 0.1, fD = 0.99: 2,113 / 2,560 / 3,002 / 4,024 months at k/k_det = 0.98 / 0.99 / 0.995 / 0.999, and ≈ 16,300 at k = k_det). So at baseline, and for these initial conditions (N₀ ≤ 0.09, E₀ = E_f), the no-event boundary coincides with k_det and T diverges there. This is verified only for these parameters and initial conditions; for other parameter sets (strongly focus-like equilibria, or initial states with depleted E) overshoot above N* could place the true boundary above k_det and has not been excluded.

For N_det = 0.1 and k ∈ (0.755, 1) k_c, A0 and A1 return finite times (Fig. G1, shaded region) although the simulation never crosses; A2 returns ∞ there and A3 returns ∞ because its tail never reaches N_det. The Phase 1 suggestion that long transients near k_c produce delayed relapse therefore does not hold for N_det = 0.1 at k > 0.755 k_c at baseline: those trajectories settle at N* < N_det.

### 6.5 Conditions of validity

| estimator | needs | fails when |
|---|---|---|
| A0 | N₀ already on the growing mode; N_det ≪ λ/|γ|; no immune feedback | dormant-dominated starts; k near k_det; (24% median error at N_det = 0.1) |
| A1 / A1b | N_det ≪ n₁λ₊/|γ| (here ≈ 0.046) and N_det ≪ φ h (0.153); k ≤ 0.9 k_det | N_det ≳ 0.02 (errors ≥ 13%); near k_det |
| A2 | N_det ≲ 0.3·min(n₁λ₊/|γ|, φh) ≈ 0.014 (here λ₊/|γ| = 0.0215, n₁ = 2.14, φh = 0.153); empirically fine at 10⁻², not at 0.1 | N_det = 0.1 (P_det/h = 0.65): returns ∞ |
| A3 | N_match ≪ λ/|γ|, φh; numerical tail | semi-analytic (tail ≈ 2/3 of T at N_det = 0.1); ∞ beyond k_det, consistent with no event |

---

## 7. Numerical verification summary

* 71 tests pass (16 Phase 1 + 55 new, `pytest`): continuation vs exact branch (≤ 10⁻¹²), residuals ≤ 10⁻¹², cubic vs multistart vs Phase-1 grid agreement, rejection of negative roots, fold location to 10⁻⁸, stability classification (including a hidden positive eigenvalue), Routh–Hurwitz vs eigenvalues on 500 random matrices, synthetic Hopf function test, critical-parameter formulas, Sotomayor conditions, eigenvalue crossing speed, stability exchange, degenerate η*, ρ*, normalised-model equivalence over gauges, A0–A3 formulas and validity, k_det divergence, step-size and solver-tolerance convergence.
* Phase 1 outputs regenerate unchanged (`baseline_scenarios.csv` byte-identical).

---

## 8. Limitations and unresolved questions

1. **Hypothetical parameters.** None of the numbers (k_c, η*, k_det, λ₊, …) has an empirical basis; the *structure* (formulas (2), the cubic, Sotomayor conditions, the scaling, the failure modes of the estimates) is what transfers.
2. **No proof of global behaviour.** Boundedness, absence of limit cycles and of chaos are checked only by scans/simulation in bounded ranges. "No Hopf found" ≠ "no Hopf exists".
3. **Scan design.** Log-uniform sampling, 8,000 draws, parameters outside the listed ranges untested; random-initial-condition oscillation search is small (263). Bistability window sizes depend on the grid (n = 401; ≈ 3×10⁻⁴ resolution at a fold).
4. **The transcritical label** is justified in the extended phase space; admissible states exist on one side only (§3.5). A fully rigorous centre-manifold reduction at higher order was not carried out; only the quadratic normal-form coefficients (which suffice for the Sotomayor conditions) were computed.
5. **Event locator accuracy at the branch point** is ~10⁻⁹–10⁻⁶ in k for steep branches (η = 0.6: error 1.5×10⁻⁶ in one run); the analytic value is exact and is used for conclusions.
6. **Relapse approximations** were tested at baseline structural parameters only for the k/N_det/composition sweep, and the "no event beyond k_det" statement (§6.4) only there; A2/A3 validity thresholds (e.g. N_det ≲ 0.02) are for this parameter set and scale with λ₊/|γ| and φh. A3 is semi-analytic.
7. **N₀ held fixed** (1.01×10⁻⁴); dependence on N₀ enters only through ln u₀ in A1–A3 (not scanned).
8. **E-gauge only partially resolved:** the unit of E is fixed, but ρ, h, η remain unidentifiable from low-burden data (Phase 1, §3), and in the normalised model the source equals μ (not a free parameter).
9. **Phase-1 documents not edited.** Phase 1 text listing "stable nonzero equilibrium at every k < k_c" and "no unstable interior dormant state" is true *at baseline*; in §4 we show the unstable interior state (saddle) and bistability appear for η ≳ 0.164 or ρ ≲ 0.005.

---

## 9. Figures

| file | content |
|---|---|
| `fig_A_equilibrium_branches.png` | A. tumour-free and nonzero branches (burden, E) vs k, with the P < 0 continuation marked |
| `fig_B_stability_bifurcation_diagram.png` | B. stability-classified diagrams, baseline (B1) and η = 0.6 (B2: fold, bistable window) |
| `fig_C_max_real_eigenvalue.png` | C. largest Re λ vs k; C2 zoom shows λ_TF ≈ α(k−k_c), λ_nz ≈ −α(k−k_c) |
| `fig_D_dormant_fraction.png` | D. equilibrium dormant fraction (identity) and transient D/N at detection |
| `fig_E_phase_portraits.png` | E. (N, E) projections at 0.9, 1.0, 1.1 k_c |
| `fig_F_additional_bifurcations.png` | F. saddle-node and bistable window vs η; regimes in (η, ρ); branch curves |
| `fig_G_relapse_time_approximation.png` | G. analytical vs simulated relapse time; error by method; initial composition |
| `fig_H_original_vs_normalized.png` | H. original vs normalised model across gauges |
