## 4. Phase 2: Equilibria, Bifurcations and Relapse-Time Analysis

**Purpose of this phase.** Phase 1 left the model frozen and produced three open questions: (i) is the observed "stable nonzero state that shrinks to the tumour-free state as $k\to k_c^-$" a genuine bifurcation, and are there other boundaries (saddle-nodes, Hopf points) that the baseline parameters hide? (ii) the immune-effector scale $E$ is a gauge (only $\beta = ks/\mu$ matters at low burden), so can the model be rewritten without it? (iii) the crude relapse-time estimate $T\approx\ln(N_{det}/N_0)/\lambda_+$ failed a first check (≈240 vs ≈375 months), so *why* does it fail and what replaces it? Phase 2 answers these with exact algebra, continuation, simulation checks and a test suite (71 tests: 16 from Phase 1 unchanged + 55 new). **Everything in this section concerns a hypothetical model with illustrative parameters; no biological or clinical validation is claimed.** Statements are labelled *proved/exact* (analytic derivation or exact computation), *numerical* (computed, verified within stated tolerance) or *assumption*.

**The model (state $P$, $D$, $e$ or $E$; time in months).**

$$
\begin{aligned}
\dot P &= rP\left(1-\tfrac{N}{K}\right) - qP + aD - kEP,\qquad N=P+D\\
\dot D &= qP - (a+\delta_D)D\\
\dot E &= s + \rho\,\tfrac{P}{h+P}E - \mu E - \eta N E
\end{aligned}
\tag{4.0}
$$

Baseline (hypothetical): $r=0.12,\ q=0.04,\ a=0.005,\ \delta_D=0.001,\ k=0.1,\ s=0.02,\ \mu=0.05,\ \rho=0.04,\ h=0.02,\ \eta=0.08,\ K=1$. Derived constants used throughout:

**Table 4.1.** Derived constants of the baseline model (source: report/phase_2_mathematical_analysis.md §0).

| symbol | definition | baseline value |
|---|---|---|
| $c$ | $q/(a+\delta_D)$ | 6.667 |
| $\varphi$ | $1+c = N^*/P^*$ on every nonzero equilibrium | 7.667 |
| $r_\mathrm{eff}$ | $r - q\delta_D/(a+\delta_D)$ | 0.11333 |
| $E_f$ | $s/\mu$ (tumour-free immune level) | 0.4 |
| $k_c$ | $r_\mathrm{eff}\mu/s$ | 0.28333 |
| $\beta$ | $ks/\mu$ (effective killing) | 0.04 at $k=0.1$ |

### 4.1 Equilibria: exact reduction and exhaustive enumeration

**Why.** Every later statement (stability, bifurcation, relapse) depends on knowing *all* equilibria, not only the one the simulator happened to find. The model is three-dimensional, so a brute-force search could miss branches; an exact reduction gives a guarantee.

**Step 1 – tumour-free state (proved).** $P=D=0$ gives $E=s/\mu$. Conversely $\dot D=0$ with $P=0$ forces $D=0$, so there is no equilibrium with $P=0,\ D>0$.

**Step 2 – nonzero equilibria (proved).** From $\dot D=0$: $D=cP$ and $N=\varphi P$. Substituting in $\dot P=0$, dividing by $P>0$ and using $-q+ac=-q\delta_D/(a+\delta_D)$:

$$
r_\mathrm{eff} - \frac{r\varphi}{K}P = kE,\qquad E=\frac{s}{B(P)},\quad B(P)=\mu+\eta\varphi P-\rho\frac{P}{h+P}.
\tag{4.1}
$$

Because $k$ appears linearly, every nonzero equilibrium lies on an **explicit graph over $P$**:

$$
k=k(P)=\frac{A(P)B(P)}{s},\quad A(P)=r_\mathrm{eff}-\frac{r\varphi}{K}P,\quad D=cP,\quad E=\frac{s}{B(P)}.
\tag{4.2}
$$

Four exact consequences:

1. **Bounded admissible range.** $E>0\Rightarrow k(P)>0\Rightarrow A(P)>0$, so $P<r_\mathrm{eff}K/(r\varphi)$; no nonzero equilibrium exists if $r_\mathrm{eff}\le 0$.
2. **At most three nonzero equilibria for any parameter set.** Multiplying (4.1) by $(h+P)$ gives the cubic $A(P)\,[(\mu+\eta\varphi P)(h+P)-\rho P]-ks(h+P)=0$.
3. **The dormant fraction is an identity, not an outcome:** $D^*/N^* = c/(1+c) = q/(q+a+\delta_D) = 0.8696$ at *every* nonzero equilibrium for every $k$. It is set by the switching rates, not by immune control.
4. **Branch crossing.** At $P=0$ the cubic's constant term is $h(r_\mathrm{eff}\mu-ks)$, which vanishes exactly at $k=k_c$: the nonzero branch passes through the tumour-free state at $k_c$.

**Step 3 – numerical verification (numerical).** `enumerate_equilibria` takes the cubic roots, rejects complex roots and roots with $P\le0$ or $E\le0$ (returning the reason), polishes with Newton on the full 3-D system and checks the residual. An *independent* multistart search (300 random log-uniform starts per $k$, no knowledge of the cubic) agrees in the number and location of equilibria at all 8 tested $k$ values (`equilibrium_enumeration_crosscheck.csv`); residuals ≤ $6\times10^{-18}$.

![Figure 4.1. Equilibrium branches versus immune killing k](../results/figures/fig_A_equilibrium_branches.png)

**Figure 4.1.** Equilibrium branches versus immune killing $k$. *Setting:* baseline parameters (footer of the figure), exact branch graph (4.2) plus continuation. *Axes/panels:* left, equilibrium burden $N^*=P^*+D^*$ (fraction of $K$) against $k\in[0,0.45]$; right, equilibrium immune level $E^*$ against $k$. Blue = tumour-free branch, orange = admissible nonzero branch ($P>0$), grey dotted = algebraic continuation with $P<0$ (non-physical); dash-dotted vertical line at $k_c=0.2833$. *What it shows:* the left panel shows the nonzero burden falling monotonically from above 0.7 at $k\approx0.09$ to 0 at $k_c$, with a steep drop near $k\approx0.18$–$0.2$ and a long shallow tail towards $k_c$; for $k>k_c$ the only admissible state is tumour-free (the continuation goes negative). In the right panel the tumour-free $E^*$ is the constant 0.4 while the nonzero branch starts near 0.22 at $k=0$, rises to a peak of about 0.48 near $k\approx0.2$ and returns to 0.4 exactly at $k_c$. *Mathematical interpretation:* the two branches meet at $k_c$ (branch point, proved via the cubic's constant term); the single nonzero branch is a graph over $P$ with no fold at baseline. The rise of $E^*$ above $E_f$ means that tumour-driven stimulation ($\rho$) outweighs burden-driven suppression ($\eta$) at moderate burden. *Biological interpretation:* weak immune killing allows a persistent, immune-boosted, high-burden state; as killing approaches the critical value the persistent state shrinks continuously to nothing. *Limitations:* hypothetical parameters; mean-field ODE; stability is not shown here (see Figure 4.3).

### 4.2 Continuation and verification of the branch structure

**Why.** The exact graph (4.2) gives the equilibria but not their stability or the exact location of bifurcation points when branches intersect. A pseudo-arclength continuation of $F(x,k)=0$, $x=(P,D,E)$, traces branches through folds and branch points and gives an independent check on the algebra.

**Method (numerical).** Predictor–corrector with analytic Jacobian and $\partial F/\partial k=(-EP,0,0)$, scaled variables. The tangent is the SVD null vector of the $3\times4$ matrix $[F_x\ F_k]$; at a branch point the null space is 2-D and the direction closest to the previous tangent is taken, so the continuation passes straight through intersections. Convergence is judged on the Newton *step* ($\le10^{-14}$), not on the residual: an early version with a residual-only test mislocated $k_c$ by $9\times10^{-5}$ because $F_P\propto P$ is tiny near the origin. Steps adapt (×1.6 after ≤ 3 iterations, ×0.6 after ≥ 6, halved on failure or tangent angle > ≈11°). Monitored functions, with sign changes located by bisection: $\det J$ (fold or branch point), $a_1a_2-a_3$ (Hopf, requiring $a_2>0$ for a complex pair; with $a_2<0$ the event is recorded as a "neutral saddle" and **not** as a bifurcation) and $P$ (leaving the orthant). Every accepted point stores all Jacobian eigenvalues, a stability label from all eigenvalues, the residual and a `physical` flag. Multi-seed tracing seeds forward/backward continuations from exact equilibria at 25 values of $k\in[0,3k_c]$, so disconnected branches would be found.

**Results at baseline** (`bifurcation_branches_baseline.csv`, `bifurcation_events_baseline.csv`): 217 points; max $|F|=8.7\times10^{-18}$; max $|k-k(P)|=2.5\times10^{-16}$ against the exact curve; one nonzero branch; the only bifurcation is the branch point at $k=0.283333333333$ (analytic value identical to 12 digits). The "neutral saddle" events at $k=0.0839$ and $0.1850$ on the tumour-free branch are points where two Jacobian eigenvalues sum to zero (e.g. $-\mu+\lambda_+=0$), real eigenvalues, not a Hopf. Reducing the step cap from 0.2 to 0.01 (`continuation_convergence.csv`) leaves $k(P)$ accuracy unchanged ($\approx2\times10^{-16}$) and the branch-point error ≤ $7\times10^{-10}$.

![Figure 4.3. Stability-classified bifurcation diagrams](../results/figures/fig_B_stability_bifurcation_diagram.png)

**Figure 4.3.** Stability-classified bifurcation diagrams. *Setting:* B1 baseline; B2 identical except $\eta=0.6$. Continuation points. *Axes/panels:* equilibrium burden $N^*$ against $k$; solid = stable, dashed = unstable, dotted = non-physical ($P<0$); blue = tumour-free, orange = nonzero; circle = fold, square = transcritical point. *What it shows:* **B1 (baseline):** the tumour-free state is unstable (dashed) for $k<k_c$ and stable (solid) beyond; the nonzero branch is entirely stable and ends at the transcritical point $k_c$; no fold, no other change of stability. **B2 ($\eta=0.6$):** the nonzero branch bends backwards: from the transcritical point (at the same $k_c$) it runs to *larger* $k$ as an unstable dashed branch, reaches a fold near $k\approx0.86$ ($N\approx0.45$) and continues back as a stable upper branch. The shaded window $k_c<k<k_{SN}$ is marked "bistable window (tumour-free and nonzero both stable)". *Mathematical interpretation:* the fold is at $k_{SN}=0.860834$ ($P=0.058486$, $N=0.4484$), matching the exact cubic root to $10^{-8}$; the backward direction corresponds to $k'(0)>0$ (Section 4.3). Inside the window the unstable middle branch is a saddle that separates the two basins. *Biological interpretation:* in B2 a patient-like system with enough killing to make the tumour-free state locally stable can still sustain a high-burden state if it starts above the saddle: local stability of "cure" does not imply global control. *Limitations:* B2 is a deliberately changed parameter set ($\eta=0.6$ is 7.5× baseline); baseline itself is monostable.

### 4.3 The transcritical bifurcation at $k=k_c$

**Why.** Phase 1 only *observed* that the stable nonzero state vanishes at $k_c$. To call this a transcritical bifurcation one must prove a set of conditions, and to know what may be claimed one must check which side is physical.

**Definition used (Sotomayor/Kuznetsov).** At $(x_0,k_0)$ with a *simple* zero eigenvalue of $J=F_x$, right/left null vectors $w,\ell$ ($\ell\cdot w=1$) and other eigenvalues off the imaginary axis: (1) $\ell\cdot F_k=0$; (2) transversality $\alpha=\ell\cdot[F_{xk}w]\ne0$; (3) nondegeneracy $\beta_\mathrm{nd}=\tfrac12\ell\cdot D^2F(w,w)\ne0$. The centre-manifold equation is then $\dot u=\alpha(k-k_0)u+\beta_\mathrm{nd}u^2+\dots$, whose two branches $u=0$ and $u^*=-\alpha(k-k_0)/\beta_\mathrm{nd}$ intersect at $k_0$ and exchange stability.

**Derivation (proved).** At $x_0=(0,0,E_f)$ the Jacobian is block lower-triangular (the $E$-column vanishes in the $P,D$ rows), so its eigenvalues are $\{-\mu\}\cup\mathrm{eig}(M)$ with $M=\begin{pmatrix}r-q-\beta&a\\q&-(a+\delta_D)\end{pmatrix}$, $\beta=kE_f$. $\det M=0\iff\beta=r_\mathrm{eff}$, i.e. $k_c=r_\mathrm{eff}\mu/s=0.28333$; there $\mathrm{tr}\,M<0$, so the other eigenvalues are $-\mu=-0.05$ and $-0.03933$. Null vectors: $w=(1,c,w_E)$, $w_E=11.093$; $\ell=(0.15254,0.12712,0)$.

**Table 4.2.** Sotomayor conditions at the baseline transcritical point (source: results/tables/transcritical_verification.csv; report §3.2).

| condition | value |
|---|---|
| $\ell\cdot F_k$ | 0 ($F_k=(-EP,0,0)=0$ at $P=0$) |
| transversality $\alpha=-E_f\ell_Pw_P$ | $-0.061017\ne0$ |
| nondegeneracy $\beta_\mathrm{nd}$ | $-0.61980\ne0$ |
| $d\lambda/dk$ (finite differences of the Jacobian eigenvalue) | $-0.0610169$ (= $\alpha$ to 7 digits) |
| zero eigenvalue simple, others $\mathrm{Re}<0$ | yes: $\{-0.05,\,0,\,-0.03933\}$ |
| $k_c$: analytic / continuation | 0.2833333333333333 / 0.2833333333333385 |

$D^2F$ was checked against symbolic differentiation (sympy) in the tests.

**Branch slope (proved + numerical).** From (4.2), $k'(0)=[-r\varphi\mu/K+r_\mathrm{eff}(\eta\varphi-\rho/h)]/s=-10.158$, hence near $k_c$: $P^*(k)=-0.098447\,(k-k_c)+\dots$. The centre-manifold prediction $-\alpha w_P/\beta_\mathrm{nd}=-0.0984467294$ equals $1/k'(0)$ to 16 digits; algebraically $\beta_\mathrm{nd}=-\alpha w_Pk'(0)$, so *nondegeneracy $\iff k'(0)\ne0$*. The numerical branch agrees with this asymptote with relative error 17.9% at $\varepsilon=0.1$, $1.7\times10^{-3}$ at $\varepsilon=10^{-3}$ and $1.7\times10^{-4}$ at $\varepsilon=10^{-4}$ (first-order convergence, as expected; `transcritical_asymptotics.csv`).

**Stability exchange (numerical).** For all 12 sampled $k=k_c(1\pm\varepsilon)$, $\varepsilon\in[10^{-4},10^{-1}]$, exactly one of the two branches is stable. For $k<k_c$ the tumour-free state is a saddle ($\max\mathrm{Re}\lambda=\alpha(k-k_c)>0$) and the nonzero state is stable; for $k>k_c$ the tumour-free state is stable and the algebraic continuation has $P<0$. The ratios $\lambda_{TF}/[\alpha(k-k_c)]$ and $\lambda_{nz}/[-\alpha(k-k_c)]$ tend to 1 (1.0002 and 0.9997 at $\varepsilon=10^{-4}$).

![Figure 4.4. Largest real part of the Jacobian eigenvalues](../results/figures/fig_C_max_real_eigenvalue.png)

**Figure 4.4.** Largest real part of the Jacobian eigenvalues versus $k$. *Setting:* baseline parameters. *Axes/panels:* C1: $\max\mathrm{Re}\,\lambda$ (1/month) against $k\in[0,0.45]$; C2: zoom on $k\in[0.255,0.312]$ (note the tenfold smaller vertical scale, ±0.002). Blue = tumour-free, orange = admissible nonzero branch, grey dotted = $P<0$ continuation. *What it shows:* C1: the tumour-free branch starts at about $+0.083$/month at $k=0$, decreases monotonically and crosses zero at $k_c$; the nonzero branch is slightly negative (about $-0.02$ at $k=0$), becoming nearly flat and close to zero between $k\approx0.18$ and $k_c$, where it touches zero. C2: the two curves cross linearly at $k_c$ with opposite slopes: $\lambda_{TF}\approx\alpha(k-k_c)$ and $\lambda_{nz}\approx-\alpha(k-k_c)$. *Mathematical interpretation:* this is the stability-exchange signature of a transcritical point with transversality $\alpha=-0.061$; the weak negative real part of the nonzero branch over a long range reflects an eigenvalue close to zero (slow relaxation, possibly oscillatory since the equilibrium is a focus over much of the range). *Biological interpretation:* near the critical killing strength the system takes very long to settle, which is why transient escape and delayed relapse are possible even when a stable state exists. *Limitations:* only the real part of the leading eigenvalue is shown, not the frequency; baseline only.

**What may and may not be claimed (proved).** The point is a transcritical bifurcation of the *extended* system ($P\in(-h,\infty)$). But on the nonzero branch $D=cP$, so it passes through the origin of the $(P,D)$ plane and $P$ changes sign with $k-k_c$ whenever $k'(0)\ne0$: an exchange with admissible branches on both sides is impossible. In the orthant the picture is a **forward ("supercritical") branching-off of a stable low-burden state from the tumour-free state**: below $k_c$ a stable low-burden state exists; above $k_c$ only the tumour-free state is admissible. The direction is not universal: $k'(0)>0$ (backward bifurcation) for $\eta>\eta^*=[\rho/h+r\varphi\mu/(Kr_\mathrm{eff})]/\varphi=0.3138$ or $\rho<\rho^*=0.00415$ (other parameters baseline); at $\eta=\eta^*$ the point is degenerate ($\beta_\mathrm{nd}=0$, codimension 2).

**Dormant fraction along trajectories.** The identity of Section 4.1 concerns equilibria only; Figure 4.2 checks it numerically and contrasts it with the transient composition at the moment of relapse.

![Figure 4.2. Dormant fraction: equilibrium identity and transient composition](../results/figures/fig_D_dormant_fraction.png)

**Figure 4.2.** Dormant fraction: equilibrium identity versus transient composition at detection. *Setting:* baseline; D2 uses $P_0=10^{-6},D_0=10^{-4},E_0=s/\mu$ (initial $D_0/N_0=0.990$) and $N_{det}=0.1$. *Axes/panels:* D1: equilibrium $D^*/N^*$ against $k$; D2: $D/N$ at first crossing of $N_{det}$ against $k$ (points missing where no crossing occurs within $2\times10^4$ months). *What it shows:* D1: the computed nonzero equilibria (orange) lie exactly on the dashed identity line $q/(q+a+\delta_D)=0.8696$ for all $k<k_c$. D2: the dormant fraction at relapse rises from about 0.33 at $k\approx0.006$ to 0.87 at $k\approx0.21$, always below the initial 0.99, and approaches the equilibrium value 0.87 from below as $k$ grows (relapse then happens only after a long, near-equilibrium transient). *Mathematical interpretation:* D1 is an exact algebraic identity (proved), D2 is numerical. The composition at detection is therefore *not* a property of the final state but of the transient: weak killing means fast proliferation and a proliferating-rich tumour at detection. *Biological interpretation:* a dormant-dominated residual population (99% dormant) can relapse as a mostly proliferating tumour if killing is weak; the dormant fraction observed at relapse should not be read as the dormant fraction of the persistent state. *Limitations:* D2 depends on the arbitrary threshold $N_{det}=0.1$ and initial condition.

**Phase portraits across $k_c$.** Figure 4.5 shows the flow itself rather than its linearisation.

![Figure 4.5. Phase portraits across the transcritical boundary](../results/figures/fig_E_phase_portraits.png)

**Figure 4.5.** Projections of the 3-D flow onto $(N,E)$ for $k=0.9k_c,\ k_c,\ 1.1k_c$. *Setting:* baseline, trajectories integrated to 3000 months from a grid of initial states (dots), orange curves starting at $D/N=0$, blue at $D/N=0.87$. *Axes/panels:* burden $N$ (0–0.07) against $E$ (0.32–0.50). Left ($k=0.255$): open square = tumour-free saddle, filled orange circle = stable nonzero equilibrium at $N\approx0.026$, $E\approx0.43$; middle ($k=k_c$): tumour-free state marginal (open square), no nonzero state; right ($k=0.3117$): filled blue square = stable tumour-free state. *What it shows:* in all three panels trajectories first move quickly in $E$ (vertical lines to the quasi-steady immune level near $E\approx0.41$–$0.44$) and then drift slowly along a common curve in $N$. At $0.9k_c$ they converge to the nonzero point; at $k_c$ they creep towards the origin along the slow curve; at $1.1k_c$ they converge to the origin. Crossings in the projection are not intersections in 3-D. *Mathematical interpretation:* a fast–slow structure: $E$ is slaved to $N$ on a slow manifold, and the single stable equilibrium moves continuously onto the origin at $k_c$. No closed orbits or spirals are visible, consistent with the absence of Hopf points (Section 4.5). *Biological interpretation:* immune effectors equilibrate quickly relative to tumour burden; the tumour burden then decides the outcome slowly, over years. *Limitations:* a 2-D projection of a 3-D flow; finite set of initial conditions.

### 4.4 Additional bifurcations: saddle-nodes and the two kinds of bistability

**Why.** The baseline is structurally simple, but a model whose conclusions might change under modest parameter changes must be mapped. Because (4.2) is a graph over $P$, a *complete* scan over all $k$ at once is possible.

**Method (numerical, exact per sample).** For each parameter draw, $P$ is sampled over the admissible range; the Jacobian, its characteristic polynomial $\lambda^3+a_1\lambda^2+a_2\lambda+a_3$ and three diagnostics are computed: sign changes of $dk/dP$ (saddle-nodes; exact roots of a cubic), Routh–Hurwitz stability and $H=a_1a_2-a_3$ (a sign change with $a_2>0$ is a purely imaginary pair, i.e. a Hopf candidate). Scan A: 8,000 Latin-hypercube samples (log-uniform) of $a\in[10^{-4},1]$, $q\in[10^{-3},1]$, $\rho\in[10^{-3},2]$, $\eta\in[10^{-3},5]$ (7,278 valid, $r_\mathrm{eff}>0$). Scan B adds $r\in[0.02,0.5]$, $\delta_D\in[10^{-4},0.1]$, $\mu\in[0.01,0.5]$, $h\in[10^{-3},0.5]$ (6,777 valid).

**Table 4.3.** Outcomes of the two parameter scans (source: results/tables/parameter_scan_summary.csv; report §4.1).

| finding | Scan A | Scan B |
|---|---|---|
| valid parameter sets | 7,278 | 6,777 |
| forward / backward at $k_c$ | 5,602 / 1,676 | 5,382 / 1,395 |
| at least one saddle-node | 3,078 | 2,638 |
| two saddle-nodes | 850 | 773 |
| tumour-free **and** a nonzero state both stable at some $k$ | 2,580 | 2,212 |
| $\ge2$ stable *nonzero* states at one $k$ (window always below $k_c$) | 1,349 | 1,206 |
| Hopf candidates | **0** | **0** |
| "neutral saddle" (not a bifurcation) | 672 | 481 |

**Kind (a): tumour-free + nonzero bistability (backward transcritical).** For $\eta>\eta^*$ or $\rho<\rho^*$ the nonzero branch bends to the right of $k_c$. Example $\eta=0.6$ (Figure 4.3, B2): fold at $k_{SN}=0.860834$; for $k\in(0.2833,0.8608)$ the stable tumour-free state coexists with a stable high-burden state and an unstable saddle. *Simulation check* (`bistability_simulations.csv`, $\eta=0.6$): at $k=0.572$ and $0.775$ an initial condition with small residual disease ($P_0=10^{-6},D_0=10^{-4},E_0=E_f$) converges to the tumour-free state, whereas $(0.3,0.3,E=0.1)$ converges to the nonzero state ($P^*=0.0436$ and $0.0524$); at $k=0.5k_c$ both go to the same nonzero state.

**Kind (b): two stable nonzero states while the tumour-free state is unstable (fold pair below $k_c$).** Verified with root enumeration and eigenvalue labels: $\eta=0.12$ has folds at $k=0.22933$ and $0.23025$, both below $k_c$; in that window three nonzero equilibria (stable, saddle, stable) coexist while the tumour-free state is a saddle. At $\eta=0.164$ folds lie at $0.2561$ (below $k_c$) and $0.2840$ (above): on $(0.2561,0.2833)$ a stable focus ($P=0.0026$) and a stable node ($P=0.065$ at $k=0.27$) coexist with a saddle ($P=0.027$); only on $(0.2833,0.2840)$ do tumour-free and one nonzero state coexist. The first fold appears at $\eta\approx0.116$ and kind-(a) bistability onsets at $\eta\approx0.164$ (or $\rho\lesssim0.005$ at baseline $\eta$). The baseline is about a factor 2 in $\eta$ below the onset. Kind (b) occurs in about 18% of the 14,055 valid sampled sets (1,349 + 1,206), always with the window below $k_c$. Example (Scan A): $a=0.306,\ q=1.6\times10^{-3},\ \rho=0.240,\ \eta=1.637,\ k=0.15$ has stable equilibria at $P=0.0027$ (focus) and $P=0.977$ (node) with a saddle at $P=0.111$; 30 random initial conditions integrated for $3\times10^5$ months end only at one of the two stable states. The mechanism is strong immune stimulation ($\rho/\mu=4.8$).

![Figure 4.6. Saddle-node windows and regimes in the (eta, rho) plane](../results/figures/fig_F_additional_bifurcations.png)

**Figure 4.6.** Saddle-node and coexistence windows. *Setting:* all other parameters at baseline; the footer notes that no Hopf point was found in the scanned ranges, so none is drawn. *Axes/panels:* **F1** critical killing $k$ against immune suppression $\eta\in[0,1.5]$: orange line = largest saddle-node $k_{SN}(\eta)$, blue horizontal line = $k_c$ (independent of $\eta$), lavender shading = tumour-free and nonzero both stable, green sliver = two stable nonzero states; vertical dash-dot line at $\eta^*=0.314$ ("backward for $\eta>\eta^*$"). **F2** regime map in $(\eta,\rho)$: grey = no coexistence of stable states, green = two stable nonzero states ($k<k_c$), light lavender = tumour-free + one nonzero, dark blue = both kinds, black line = $k'(0)=0$ i.e. $\rho^*(\eta)$; star = baseline. **F3** nonzero branch $N^*$ against $k$ for $\eta=0,0.2,0.4,0.6,1.0$ (curves only). *What it shows:* F1: $k_{SN}$ rises roughly linearly with $\eta$ (above 2 at $\eta=1.5$), so the coexistence window widens quickly once $\eta$ exceeds about 0.16; the green window below $k_c$ is very thin. F2: the baseline star lies in the grey region, close to the green/blue boundary; coexistence regions occupy most of the plane at large $\eta$ and small $\rho$. F3: the branch is monotone for $\eta=0$, becomes S-shaped (two folds) for $\eta=0.2$ and bends increasingly far to the right (backward) for $\eta\ge0.4$. *Mathematical interpretation:* folds are roots of the cubic discriminant structure; $k'(0)=0$ is the codimension-2 line separating forward and backward branching. *Biological interpretation:* strong burden-driven immune suppression (large $\eta$) or weak tumour-driven immune stimulation (small $\rho$) makes control fragile and history-dependent; strong stimulation can create a second, high-burden attractor. *Limitations:* two-parameter slices of a 10-parameter model; the window size depends on the scan grid ($n=401$; resolution ≈ $3\times10^{-4}$ at a fold).

### 4.5 Hopf bifurcation and oscillations: searched for, not found

**Why.** Immune–tumour models can oscillate. A Hopf bifurcation would produce sustained cycles and radically change the meaning of "relapse time".

**Evidence (numerical; absence is not proved).**

* Along the exact branch graph in both scans: **no sign change of $a_1a_2-a_3$ with $a_2>0$** among 14,055 valid parameter sets. Stable *foci* (damped oscillation) occur in 354 of the examined sets (every 4th draw, first 50 stable-branch points each), but none crosses the imaginary axis.
* One Hopf candidate appears in the continuation of the *non-physical* $P<0$ continuation at $\eta=0.6$ ($k=0.2396$, $P=-0.0046$, eigenvalues $\pm0.0061i$). It is not an admissible state and is not reported as a model bifurcation.
* **Direct simulation search** (`oscillation_search*.csv`): 263 random (parameter set, $k\in[0,2k_c]$, initial condition) draws integrated to $1.5\times10^4$ months; 15 were flagged by a crude amplitude criterion; re-integrating each to $3\times10^6$ months (Radau) showed all converge to an equilibrium (residual ≤ $10^{-17}$; extrema at relative amplitude $10^{-10}$ are round-off). **0 persistent oscillations.**

**Conclusion.** "No Hopf bifurcation and no sustained oscillation in the scanned ranges." This is *not* a proof of absence: ranges are bounded and sampled, no Dulac/Lyapunov argument excludes limit cycles, and global boundedness for $\rho>\mu$ is not proved (no sampled trajectory diverged). The phase portraits (Figure 4.5) show no spiralling or closed orbits in the baseline regime.

### 4.6 M1: normalising the immune-effector scale

**Why.** In the original model the map $(E,s,k)\to(cE,cs,k/c)$ leaves $P$ and $D$ unchanged: $E$ has no natural unit and $k,s$ cannot be separated (a structural redundancy found in Phase 1). Priors and sensitivity analyses in later phases must not be placed on gauge-dependent quantities.

**Transformation (proved).** Let $e=\mu E/s$ (tumour-free level $e=1$) and $\beta=ks/\mu$. Then $kE=\beta e$ and $\dot e=(\mu/s)\dot E$ give

$$
\begin{aligned}
\dot P &= rP\left(1-\tfrac{N}{K}\right)-qP+aD-\beta eP\\
\dot D &= qP-(a+\delta_D)D\\
\dot e &= \mu(1-e)+\rho\tfrac{P}{h+P}e-\eta Ne
\end{aligned}
\tag{4.3}
$$

Parameters $k$ and $s$ disappear (10 parameters left: $K,r,q,a,\delta_D,\beta,\mu,\rho,h,\eta$). Consequences: $\beta_c=r_\mathrm{eff}=0.11333$ and $\beta/\beta_c=k/k_c$; the Jacobians are conjugate by $S=\mathrm{diag}(1,1,\mu/s)$, so eigenvalues, bifurcation points and branch burdens are identical. Baseline corresponds to $\beta=0.04$.

**Implementation.** A separate dataclass and independent right-hand side (`normalized.py`, `configs/normalized.yaml`); read-only properties `beta`, `r_eff`, `phi` were added to `Params`. No existing parameter, config or Phase 1 output was changed; `baseline_scenarios.csv` regenerates byte-identical.

**Verification (numerical).** The original model with $s\in\{0.002,0.02,0.2,2\}$ (a $10^3$ range, $k$ rescaled to keep $\beta$ fixed) at $k/k_c\in\{0.25,0.5,0.75\}$ gives relapse times that differ from the normalised model by ≤ $8.4\times10^{-7}$ months (solver tolerance), state differences at the event ≤ $7.6\times10^{-12}$ and $E=(s/\mu)e$ to $10^{-8}$ relative; RHS agreement $10^{-13}$. A negative control (changing $\beta$ by 7%) shifts $T$ by more than 1 month, so the test is not vacuous.

![Figure 4.7. Original versus normalised model across gauges](../results/figures/fig_H_original_vs_normalized.png)

**Figure 4.7.** Original versus normalised model across immune-scale gauges. *Setting:* baseline structure, $(k,s)\to(k/c,cs)$ with $\beta=ks/\mu$ fixed ($\beta=0.0567$ i.e. $k/k_c=0.5$ shown in H1/H2), $s\in\{0.002,0.02,0.2,2\}$ with $k=1.42,0.142,0.0142,0.00142$. *Axes/panels:* **H1** burden $N$ (log) against time (0–380 months), four original-model curves overlaid by the dashed normalised curve; **H2** immune level $E$ (log, arbitrary units) against time, four curves at levels differing by factors of ten, with the dashed $e=\mu E/s$ matching each; **H3** absolute relapse-time difference (original − normalised, months, log scale $10^{-13}$–$10^{-6}$) against $\beta/\beta_c=k/k_c\in\{0.25,0.5,0.75\}$ for the four values of $s$. *What it shows:* the four burden curves are visually indistinguishable from the normalised curve, rising from $10^{-4}$ to the threshold 0.1 at ≈375 months; the $E$ trajectories are rescaled copies; differences in H3 range from about $10^{-13}$ to $9\times10^{-7}$ months, growing with $k/k_c$ and largest for the smallest $s$. *Mathematical interpretation:* confirms the exact equivalence (4.3); the small differences are integration error, largest where the stiffness/scale mismatch is greatest. *Biological interpretation:* the absolute size of the immune "effector" count is arbitrary; only the dimensionless killing strength $\beta$ and the relative immune level $e$ carry meaning. *Limitations:* M1 fixes the unit of $E$ only; $\rho,h,\eta$ remain unidentifiable from low-burden data and the immune source is tied to $\mu$ (Phase 4 returns to this).

### 4.7 Relapse-time approximation: why the first estimate failed and what replaces it

**Setting.** $x_0=(P_0,D_0,E_f)$ with $N_0=1.01\times10^{-4}$ ($P_0=10^{-6},D_0=10^{-4}$), threshold $N_{det}$, $k<k_c$. The numerical truth is the first upward crossing of $N=N_{det}$ computed with DOP853 (rtol $10^{-9}$), horizon $10^5$ months. Let $\lambda_+$ be the positive real eigenvalue of $J$ at the tumour-free state, $w$ its right eigenvector ($w_P=1$), $\ell$ the left eigenvector, $n_1=w_P+w_D$ and $u_0=\ell\cdot(x_0-x_{TF})$.

**Estimators.**

* **A0 (Phase 1):** $T=\ln(N_{det}/N_0)/\lambda_+$ — assumes all of the initial burden is on the growing mode.
* **A1 (eigen-projected):** $T=\ln\!\big(N_{det}/(n_1u_0)\big)/\lambda_+$. (A1b solves the two-mode linear solution exactly.)
* **A2 (second-order unstable manifold):** on $x=x_{TF}+uw+u^2h_2$, $\dot u=\lambda_+u+\gamma u^2$ is a Bernoulli equation; with $v=1/u$,
$$T=\frac{1}{\lambda_+}\ln\frac{1/u_0+\gamma/\lambda_+}{1/u_{det}+\gamma/\lambda_+},\tag{4.4}$$
where $u_{det}$ is the smallest positive root of $N(u)=N_{det}$ and $\gamma=\tfrac12\ell\cdot D^2F(w,w)$. If the bracket is non-positive, $u$ saturates before reaching $u_{det}$ and the estimator returns $\infty$ ("no relapse predicted").
* **A3 (matched; semi-analytic, not a closed form):** $T=T_{A2}(N_{match})+\tau$, $N_{match}=10^{-3}$, with the tail $\tau$ computed by numerical integration. At baseline, $k=0.5k_c$, $N_{det}=0.1$ the closed-form part is 58.4 / 77.7 / 123.5 months for initial dormant fractions 0 / 0.5 / 0.99, while $\tau=252.046$ months for all three; the tail is 81% / 76% / 67% of $T$. So A3's small error at $N_{det}=0.1$ is largely numerics agreeing with numerics. Its value is structural: $T=T_{lin}(x_0)+\tau(\text{parameters})$, i.e. the initial condition enters only through $\ln u_0$.

**Why the Phase 1 estimate failed (baseline, $k=0.5k_c$, $N_{det}=0.1$, 99% dormant start).** $\lambda_+=0.02904$/month; other eigenvalues $-0.01171$ and $-0.05$. Numerical $T=375.5$ months versus A0 $=237.5$ (gap 138.0).

**Table 4.4.** Attribution of the gap between the numerical relapse time and the Phase 1 estimate (source: results/tables/relapse_time_attribution.csv; report §6.2).

| contribution | months | evidence |
|---|---|---|
| Only part of the 99%-dormant start sits on the growing mode (A1 − A0) | **+44.0** | A1 = 281.5; $\ell=(0.860,0.123,0)$ |
| Logistic saturation / other no-feedback nonlinearity (numerical at $\rho=\eta=0$ minus A1) | **+11.5** | $\rho=\eta=0$: numerical 293.0, A1 281.6, A2 293.4 |
| Immune feedback during growth (baseline numerical minus numerical at $\rho=\eta=0$) | **+82.5** | $E$ rises to $1.40E_f$ at detection and slows growth |
| **Total gap** | **138.0** | 375.5 − 237.5 |

The split is order-dependent; the robust statement is that the initial-state projection explains about a third of the gap and the immune response about 60%. Variants: $\rho=0$ only: numerical 287.2 = A2 287.2; $\eta=0$ only: numerical 459.9 ($E$ reaches $1.59E_f$). Quadratic closure also breaks down because $P_{det}/h=N_{det}/(\varphi h)=0.65$ at $N_{det}=0.1$ (the term $\rho P/(h+P)$ lies outside its Taylor radius), which is why A2 returns $\infty$ at $N_{det}=0.1$ but works at $10^{-2}$ ($P_{det}/h=0.065$).

**Table 4.5.** Relapse times in months at $k=0.5k_c$ (numerical / A0 / A1 / A2 / A3) for three initial dormant fractions (source: results/tables/relapse_time_approximation.csv).

| $N_{det}$ | $D_0$ fraction 0 | 0.5 | 0.99 |
|---|---|---|---|
| $10^{-3}$ | 59.8 / 78.9 / 57.9 / 58.4 / 58.4 | 77.8 / 78.9 / 77.2 / 77.7 / 77.7 | 122.9 / 78.9 / 123.0 / 123.5 / 123.5 |
| $10^{-2}$ | 143.5 / 158.2 / 137.2 / 143.8 / 143.4 | 162.7 / 158.2 / 156.5 / 163.1 / 162.7 | 208.4 / 158.2 / 202.2 / 208.9 / 208.5 |
| 0.1 | 310.5 / 237.5 / 216.5 / ∞ / 310.4 | 329.7 / 237.5 / 235.8 / ∞ / 329.8 | 375.5 / 237.5 / 281.5 / ∞ / 375.6 |

**Table 4.6.** Relative error of each estimator over $k/k_c\in[0.05,\,0.9\,k_{det}/k_c]$ and three compositions, median / max (source: results/tables/relapse_time_accuracy_summary.csv).

| $N_{det}$ | A0 | A1 | A2 | A3 (semi-analytic) |
|---|---|---|---|---|
| $10^{-3}$ | 20% / 49% | 1.4% / 3.2% | 0.6% / 2.4% | 0.6% / 2.4% |
| $10^{-2}$ | 8.7% / 33% | 2.9% / 13% | 0.19% / 3.4% | 0.054% / 0.14% |
| 0.1 | 24% / 44% | 10% / 42% | 1.7% / 11% (9 of 21 predict no event) | 0.007% / 0.036% |

**Honest summary at $N_{det}=0.1$:** the best *closed-form* estimate is A1 (10% median, 42% maximum error); A2 is a good closed form only for $N_{det}\lesssim10^{-2}$ (≤ 3.4%); A3's 0.007% includes a numerically integrated tail covering 76% (median) of $T$. A0 is also independent of the initial composition: at $N_{det}=10^{-3}$ it overestimates the all-proliferating start by 32% and underestimates the 99%-dormant start by 36%.

![Figure 4.8. Analytical estimates versus simulated relapse time](../results/figures/fig_G_relapse_time_approximation.png)

**Figure 4.8.** Analytical relapse-time estimates against simulation. *Setting:* baseline, E$_0=s/\mu$, numerical = DOP853 first upward crossing (rtol $10^{-9}$); G1 at $N_{det}=0.1$, initial dormant fraction 0.99. *Axes/panels:* **G1** relapse time (months, log scale, $10^2$–$>10^3$) against $k/k_c$: black dots = numerical, grey dashed = A0, blue dash-dot = A1, red dotted = A2, green = A3; pink region for $k>k_{det}$ ($k/k_c>0.755$) labelled "$N^*<N_{det}$, no event in simulation; A0/A1 still return finite $T$". **G2** relative error (estimate − numeric)/numeric in % against $k/k_c$, line style = threshold (solid $10^{-3}$, dashed $10^{-2}$, dotted 0.1). **G3** relapse time against initial dormant fraction $D_0/(P_0+D_0)$ at $k=0.5k_c$, fixed $N_0=1.01\times10^{-4}$. *What it shows:* G1: numerical $T$ grows from about 120 months at $k/k_c=0.05$ to over $10^3$ near 0.7–0.75 and diverges towards $k_{det}$; A3 follows the numerical points exactly, A1 and A0 fall progressively below, A2 stops (returns $\infty$) once beyond its validity. G2: A0 errors are −26% to −50% (it underestimates), A1 stays within about ±2% at $10^{-3}$ and drifts to −37% at $N_{det}=0.1$, A2 is within about 3% for $10^{-3}$/$10^{-2}$ and up to +8% at 0.1 where it exists, A3 ≈ 0. G3: numerical $T$ rises from about 310 to 375 months as the start becomes more dormant; A0 is flat at 237.5 (ignores composition); A1 rises from about 216 to 281 but stays below. *Mathematical interpretation:* the error budget of Table 4.4 and the validity conditions below; divergence of $T$ near $k_{det}$ is the approach to the saddle-node-like slow-down of reaching an equilibrium below the threshold. *Biological interpretation:* a residual population that is mostly dormant takes substantially longer to relapse than one that is proliferating, and immune stimulation during growth delays relapse further. *Limitations:* baseline only; N$_0$ not varied; A3 is semi-analytic.

**Conditions of validity (summary).** A0: $N_0$ already on the growing mode, $N_{det}\ll\lambda/|\gamma|$, no immune feedback. A1: $N_{det}\ll n_1\lambda_+/|\gamma|\ (\approx0.046)$ and $\ll\varphi h\ (0.153)$, $k\le0.9k_{det}$. A2: $N_{det}\lesssim0.3\min(n_1\lambda_+/|\gamma|,\varphi h)\approx0.014$; empirically fine at $10^{-2}$, not at 0.1. A3: requires $N_{match}\ll\lambda/|\gamma|,\varphi h$ and a numerical tail. These thresholds were measured at baseline and scale with $\lambda_+/|\gamma|$ and $\varphi h$.

### 4.8 The validity boundary of relapse is $k_{det}$, not $k_c$

**Why.** For $k\in(k_{det},k_c)$ the stable nonzero equilibrium has $N^*<N_{det}$: the burden approaches a level below the detection threshold and never crosses it. The boundary $k_{det}$ solves $N^*(k)=N_{det}$, i.e. $k_{det}=k(P=N_{det}/\varphi)$ from (4.2).

**Values (exact).** $k_{det}/k_c=0.9954,\ 0.9568,\ 0.7548$ for $N_{det}=10^{-3},10^{-2},0.1$.

**Is "no event beyond $k_{det}$" exact? (numerical test).** It requires the approach to $N^*$ to be monotone from below, but the equilibrium is a stable *focus* over much of the range (weak damping, Figure 4.4), so overshoot is conceivable. A fine scan (`relapse_time_kdet_boundary.csv`): $k/k_{det}\in\{0.98,\dots,1.06\}$ (12 values) × initial dormant fractions $\{0,0.5,0.99\}$ × $N_{det}\in\{10^{-3},10^{-2},0.1\}$, horizon $10^5$ months. **0 of 54 runs with $k>1.001k_{det}$ produced an event**; the maximum $N$ along those trajectories stays below $N^*$ (e.g. 0.088 at $k=1.02k_{det}$ for $N_{det}=0.1$). Event time increases as $k\to k_{det}^-$: for $N_{det}=0.1$ and 99% dormant, 2,113 / 2,560 / 3,002 / 4,024 months at $k/k_{det}=0.98/0.99/0.995/0.999$ and ≈16,300 months at $k_{det}$. (For $N_{det}=10^{-3}$, no dormant, the scan gives 5,073 / 9,451 / 16,118 months at 0.98 / 0.99 / 0.995.)

**Consequence.** At baseline and for these initial conditions ($N_0\le0.09$, $E_0=E_f$) the no-event boundary coincides with $k_{det}$ and $T$ diverges there. For $N_{det}=0.1$ and $k\in(0.755,1)k_c$, A0 and A1 return finite times although the simulation never crosses (shaded region of Figure 4.8). This corrects the Phase 1 suggestion that long transients near $k_c$ produce delayed relapse: at $N_{det}=0.1$ and $k>0.755k_c$ the trajectories settle at $N^*<N_{det}$. The statement is verified only at baseline; a strongly focus-like equilibrium or an immune-depleted start could place the true boundary above $k_{det}$ and this was not excluded.

### 4.9 Summary of Phase 2 conclusions and implications for Phase 3

**Table 4.7.** Critical parameter values at baseline (hypothetical; source: PHASE_2_SUMMARY.md).

| quantity | value |
|---|---|
| $k_c$ (transcritical) | 0.28333 (continuation 0.283333333333) |
| $\beta_c=r_\mathrm{eff}$ | 0.11333 |
| $\eta^*$ ($k'(0)=0$) / $\rho^*$ | 0.3138 / 0.00415 |
| first fold appears ($\eta$) / bistability onset ($\eta$) | ≈ 0.116 / ≈ 0.164 |
| $\eta=0.6$: fold $k_{SN}$ | 0.86083 |
| $k_{det}/k_c$ at $N_{det}=10^{-3},10^{-2},0.1$ | 0.9954, 0.9568, 0.7548 |
| $k'(0)$, $dP^*/dk$ at $k_c$ | −10.158, −0.098447 |

**Conclusions (labelled).**

1. *Proved:* all nonzero equilibria lie on the explicit graph (4.2), satisfy a cubic (at most three), and have dormant fraction $q/(q+a+\delta_D)=0.8696$.
2. *Proved + numerical:* the baseline critical point $k_c=0.28333$ is a transcritical bifurcation (simple zero eigenvalue, $\alpha=-0.0610$, $\beta_\mathrm{nd}=-0.620$, centre-manifold slope equal to the exact slope to 16 digits); physically it is a forward branching of a stable low-burden state out of the tumour-free state.
3. *Numerical:* the direction depends on parameters; saddle-nodes and two kinds of bistability exist outside the baseline (about 18% of sampled valid sets have two stable nonzero states). The baseline is about a factor 2 in $\eta$ from the bistability onset.
4. *Numerical, not a proof:* no Hopf bifurcation and no oscillation found (0 candidates in 14,055 sets; 0 of 263 long simulations).
5. *Proved:* M1 removes the E-gauge exactly; equivalence verified over a 1000-fold range of the gauge.
6. *Numerical:* the Phase 1 relapse estimate fails because of the dormant-dominated start (+44 months), logistic saturation (+11.5) and especially immune feedback (+82.5); A1 is the best closed form at $N_{det}=0.1$ and A3 is only semi-analytic.
7. *Numerical (baseline only):* the relapse boundary is $k_{det}$ ($k_{det}/k_c=0.7548$ at $N_{det}=0.1$), not $k_c$.

**Implications for Phase 3.** (i) Sensitivity and Monte Carlo should be done in the normalised model with priors on $\beta$, not on $k$ or $s$. (ii) The Monte Carlo output must be the probability of an event before a horizon *and* the event time, with censoring by $k_{det}$ built in. (iii) The numerical first-crossing time is the production event time; A1 and A3 serve as interpretable checks. (iv) The prior must decide explicitly whether the multistable regimes are in its support, since in them the outcome depends on the initial burden. (v) Natural scale-free coordinates are $\beta/\beta_c$, $k_{det}$, $\varphi$, $\lambda_+$ and $N_{det}/(\varphi h)$.

**Limitations of Phase 2.** All parameter values are hypothetical (only the structure transfers); no global proof of boundedness or of the absence of limit cycles; scans are bounded and sampled; the transcritical label applies to the extended system and only quadratic normal-form coefficients were computed; event-location accuracy at steep branch points is $10^{-9}$–$10^{-6}$ in $k$ (analytic values are used for conclusions); relapse-estimator validity thresholds and the $k_{det}$ statement were tested at baseline only and $N_0$ was not varied; M1 fixes the unit of $E$ but does not make $\rho,h,\eta$ identifiable.
