"""Phase 2 experiment: continuation-based bifurcation analysis in k, transcritical verification,
search for saddle-node / Hopf / bistability / oscillations, and figures A-F.

Run:  .venv/bin/python experiments/bifurcation_analysis.py
Outputs: results/tables/*.csv, results/figures/fig_{A..F}_*.png
Everything is a mathematical statement about the hypothetical model, not a clinical result.
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).parent))
from _style import *  # noqa: E402,F403
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from dormancy import Params, critical_k, enumerate_equilibria, load_config, rhs, simulate, jacobian  # noqa: E402
from dormancy.bifurcation import (attractor_summary, coexisting_stable_nonzero, degenerate_eta,  # noqa: E402
                                  degenerate_rho, extended_root_near_zero, scan_branch, stable_overlap_windows,
                                  transcritical_report)
from dormancy.continuation import continue_branch, trace_all_branches  # noqa: E402
from dormancy.equilibria import branch_k, branch_k_prime, fold_P, multistart_equilibria, residual  # noqa: E402
from dormancy.stability import max_real_eigenvalue, stability_label  # noqa: E402

cfg = load_config(ROOT / "configs" / "baseline.yaml")
p0: Params = cfg["params"]
kc = critical_k(p0)
KMAX = 3.0 * kc
t_start = time.time()


def log(msg):
    print(f"[{time.time() - t_start:6.1f}s] {msg}", flush=True)


def branch_frame(branches, p, tag):
    rows = []
    for bi, b in enumerate(branches):
        for i in range(len(b.lam)):
            x = b.x[i]
            N = x[0] + x[1]
            rows.append({"tag": tag, "branch": bi, "k": b.lam[i], "P": x[0], "D": x[1], "E": x[2], "N": N,
                         "dormant_fraction": x[1] / N if abs(N) > 1e-14 else np.nan,
                         "max_re_eig": b.max_real[i], "eig1": b.eigs[i, 0], "eig2": b.eigs[i, 1], "eig3": b.eigs[i, 2],
                         "stability": b.labels[i], "physical": bool(b.physical[i]), "residual_max_abs_F": b.residual[i]})
    return pd.DataFrame(rows)


def event_frame(branches, tag):
    rows = []
    for bi, b in enumerate(branches):
        for e in b.events:
            rows.append({"tag": tag, "branch": bi, "kind": e.kind, "k": e.lam, "P": e.x[0], "D": e.x[1], "E": e.x[2],
                         "eigs": np.round(e.eigs, 8).tolist(), "detail": e.detail})
    return pd.DataFrame(rows)


# =========================================================================== 1. baseline continuation
log("1. continuation at baseline parameters")
br0 = trace_all_branches(p0, (0.0, KMAX), n_seed=25)
df0 = branch_frame(br0, p0, "baseline")
ev0 = event_frame(br0, "baseline")
df0.to_csv(TAB / "bifurcation_branches_baseline.csv", index=False)
ev0.to_csv(TAB / "bifurcation_events_baseline.csv", index=False)
nz_phys = df0[(df0.P > 1e-12) & df0.physical]
# verification of every continuation point against the exact explicit branch k(P)
nz_all = df0[df0.P.abs() > 1e-12]
k_err = np.abs(nz_all.k.values - branch_k(nz_all.P.values, p0))
log(f"   {len(df0)} points, max |F| = {df0.residual_max_abs_F.max():.2e}, "
    f"max |k - k(P)| = {k_err.max():.2e}, physical nonzero stable: "
    f"{sorted(set(nz_phys.stability))}")
print(ev0[ev0.kind != 'neutral saddle'][["branch", "kind", "k", "P"]].to_string())

# independent exhaustive check: cubic enumeration vs multistart Newton on a k grid
rows = []
for f in [0.0, 0.1, 0.5, 0.9, 0.99, 0.999, 1.001, 1.5]:
    p = p0.with_(k=f * kc)
    a = enumerate_equilibria(p)
    b = multistart_equilibria(p, 300, seed=2)
    rows.append({"k_over_kc": f, "n_cubic": len(a), "n_multistart": len(b),
                 "max_residual": max([residual(e, p) for e in a] + [0.0]),
                 "agree": len(a) == len(b) and all(np.allclose(x.state, y.state, rtol=1e-6, atol=1e-10) for x, y in zip(a, b))})
pd.DataFrame(rows).to_csv(TAB / "equilibrium_enumeration_crosscheck.csv", index=False)
log("   enumeration cross-check: " + str(all(r["agree"] for r in rows)))

# =========================================================================== 2. transcritical verification
log("2. transcritical verification")
rep = transcritical_report(p0)
tc_rows = [
    ("k_c analytic = (r - q dD/(a+dD)) mu / s", rep.k_c),
    ("k_c from continuation (branch point, tumour-free branch)",
     min([e.lam for b in br0 for e in b.events if e.kind.startswith("branch point")], key=lambda v: abs(v - kc))),
    ("eigenvalues of J at (x_TF, k_c)", str(np.round(rep.eigenvalues_at_kc, 8).tolist())),
    ("simple zero eigenvalue, other two Re<0", rep.simple_zero_eigenvalue),
    ("l.F_k at (x_TF,k_c)  [must be 0]", rep.F_k_projection),
    ("alpha = l.(D_x F_k) w  (transversality)", rep.alpha),
    ("beta = (1/2) l.D2F(w,w)  (nondegeneracy)", rep.beta),
    ("d(lambda_0)/dk by finite difference", np.nan),
    ("k'(0) exact branch slope dk/dP at P=0", rep.dk_dP_exact),
    ("dP/dk centre-manifold = -alpha w_P / beta", rep.dP_dk_normal_form),
    ("dP/dk exact = 1/k'(0)", rep.dP_dk_exact),
    ("direction", rep.direction),
    ("transcritical (Sotomayor conditions hold)", rep.is_transcritical),
]


def lam0(k):
    pk = p0.with_(k=k)
    return float(np.max(np.real(np.linalg.eigvals(jacobian(np.array([0, 0, pk.E_free]), pk)))))


h_ = 1e-6
tc_rows[7] = (tc_rows[7][0], (lam0(kc + h_) - lam0(kc - h_)) / (2 * h_))
pd.DataFrame(tc_rows, columns=["quantity", "value"]).to_csv(TAB / "transcritical_verification.csv", index=False)

# local asymptotics + exchange of stability (nonzero branch continued to P<0 for k>k_c)
rows = []
for eps in [1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 1e-4]:
    for sgn in (-1, +1):
        k = kc * (1 + sgn * eps)
        p = p0.with_(k=k)
        res = extended_root_near_zero(p)
        x, e, phys = res
        tf_e = np.linalg.eigvals(jacobian(np.array([0, 0, p.E_free]), p))
        P_lin = (k - kc) / rep.dk_dP_exact
        rows.append({"k_over_kc": 1 + sgn * eps, "k": k, "P_nonzero_root": x[0], "P_linear_prediction": P_lin,
                     "rel_error_linear": abs(x[0] - P_lin) / abs(x[0]), "physical": phys,
                     "TF_max_re_eig": tf_e.real.max(), "nonzero_max_re_eig": e.real.max(),
                     "TF_stable": bool(tf_e.real.max() < 0), "nonzero_stable": bool(e.real.max() < 0),
                     "lambda_TF_over_alpha_kappa": tf_e.real.max() / (rep.alpha * (k - kc)),
                     "lambda_nz_over_minus_alpha_kappa": e.real.max() / (-rep.alpha * (k - kc))})
asym = pd.DataFrame(rows)
asym.to_csv(TAB / "transcritical_asymptotics.csv", index=False)
exch = all(r["TF_stable"] != r["nonzero_stable"] for r in rows)
log(f"   exchange of stability at every sampled k != k_c: {exch}; "
    f"max rel. error of linear branch at eps=1e-3: {asym[asym.k_over_kc.sub(1).abs().round(6) == 0.001].rel_error_linear.max():.2e}")

# =========================================================================== 3. convergence of continuation
log("3. continuation convergence study")
rows = []
nz_seed = enumerate_equilibria(p0.with_(k=0.5 * kc))[0]
for ds_max in [0.2, 0.1, 0.05, 0.02, 0.01]:
    b = continue_branch(p0.with_(k=0.5 * kc), nz_seed.state, "k", +1, (0.0, KMAX), ds0=min(0.02, ds_max), ds_max=ds_max,
                        lam0=0.5 * kc)
    pos = np.abs(b.x[:, 0]) > 1e-12
    bp = [e.lam for e in b.events if e.kind.startswith("branch point")]
    rows.append({"ds_max": ds_max, "n_points": len(b.lam), "max_residual": b.residual.max(),
                 "max_abs_k_minus_exact_k_of_P": float(np.max(np.abs(b.lam[pos] - branch_k(b.x[pos, 0], p0)))),
                 "branch_point_k": bp[0] if bp else np.nan, "branch_point_error": abs(bp[0] - kc) if bp else np.nan})
pd.DataFrame(rows).to_csv(TAB / "continuation_convergence.csv", index=False)
print(pd.DataFrame(rows).to_string())

# =========================================================================== 4. search for other bifurcations
log("4. parameter scans (saddle-node / Hopf / bistability)")
rng = np.random.default_rng(20261104)


def loguni(n, lo, hi):
    u = (rng.permutation(n) + rng.random(n)) / n
    return np.exp(np.log(lo) + u * (np.log(hi) - np.log(lo)))


scan_rows = []
examples = {"hopf": [], "bistable": [], "multi_stable_nonzero": []}
for label, N, ranges in [
        ("scan_A: a,q,rho,eta varied", 8000, dict(a=(1e-4, 1.0), q=(1e-3, 1.0), rho=(1e-3, 2.0), eta=(1e-3, 5.0))),
        ("scan_B: + r,delta_D,mu,h varied", 8000, dict(a=(1e-4, 1.0), q=(1e-3, 1.0), rho=(1e-3, 2.0), eta=(1e-3, 5.0),
                                                       r=(0.02, 0.5), delta_D=(1e-4, 0.1), mu=(0.01, 0.5), h=(1e-3, 0.5)))]:
    draws = {k: loguni(N, *v) for k, v in ranges.items()}
    st = dict(n_valid=0, n_forward=0, n_backward=0, n_with_fold=0, n_two_folds=0, n_bistable=0,
              n_hopf_candidate=0, n_neutral_saddle=0, n_multiple_stable_nonzero_same_k=0, n_multi_nonzero_with_unstable_TF_only=0,
              n_complex_eigs_on_stable_branch=0)
    for i in range(N):
        p = Params(**{k: float(v[i]) for k, v in draws.items()})
        if p.r_eff <= 0:
            continue
        sc = scan_branch(p, n=401)
        if sc.P.size == 0:
            continue
        st["n_valid"] += 1
        st["n_forward" if sc.direction.startswith("forward") else "n_backward"] += 1
        st["n_with_fold"] += bool(sc.folds)
        st["n_two_folds"] += len(sc.folds) >= 2
        st["n_bistable"] += sc.bistable
        st["n_hopf_candidate"] += bool(sc.hopf_candidates)
        st["n_neutral_saddle"] += bool(sc.neutral_saddles)
        if sc.hopf_candidates and len(examples["hopf"]) < 5:
            examples["hopf"].append({k: float(v[i]) for k, v in draws.items()})
        if sc.bistable and len(examples["bistable"]) < 5 and label.startswith("scan_A"):
            examples["bistable"].append({k: float(v[i]) for k, v in draws.items()})
        # two stable nonzero equilibria at one k (nonzero-nonzero multistability), verified with the exact cubic
        found = coexisting_stable_nonzero(p, sc)
        if found:
            st["n_multiple_stable_nonzero_same_k"] += 1
            st["n_multi_nonzero_with_unstable_TF_only"] += int(all(hi <= sc.k_c * (1 + 1e-9) for lo, hi in stable_overlap_windows(sc)))
            if len(examples["multi_stable_nonzero"]) < 5:
                examples["multi_stable_nonzero"].append({k: float(v[i]) for k, v in draws.items()})
        if i % 4 == 0 and sc.stable_k_ranges:       # complex eigenvalues on stable branch? (damped oscillation only)
            J, _, _ = __import__("dormancy.bifurcation", fromlist=["x"])._branch_jacobian_arrays(sc.P[sc.n_unstable == 0][:50], p)
            if J.size:
                ev = np.linalg.eigvals(J)
                st["n_complex_eigs_on_stable_branch"] += bool(np.any(np.abs(ev.imag) > 1e-12))
    st = {"scan": label, "n_sampled": N, **st,
          "ranges": "; ".join(f"{k} in [{v[0]:g},{v[1]:g}] (log-uniform)" for k, v in ranges.items()),
          "fixed_at_baseline": ", ".join(f"{f}={getattr(p0, f)}" for f in ("K", "r", "delta_D", "s", "mu", "h") if f not in ranges)}
    scan_rows.append(st)
    log("   " + str({k: v for k, v in st.items() if k.startswith('n_')}))
pd.DataFrame(scan_rows).to_csv(TAB / "parameter_scan_summary.csv", index=False)
pd.DataFrame([{"kind": k, **e} for k, v in examples.items() for e in v]).to_csv(TAB / "parameter_scan_examples.csv", index=False)

# 2-D classification map in (eta, rho) at baseline for the other parameters
etas = np.linspace(0.0, 2.0, 81)
rhos = np.linspace(0.0, 0.2, 81)
code = np.zeros((len(rhos), len(etas)))
kfold_map = np.full_like(code, np.nan)
rows = []
for i, rho in enumerate(rhos):
    for j, eta in enumerate(etas):
        p = p0.with_(rho=rho, eta=eta)
        sc = scan_branch(p, n=301)
        # classes by which coexistence exists (independent of the direction of the transcritical point):
        # 0 none | 1 two stable NONZERO states only (k<k_c, tumour-free unstable) | 2 tumour-free + nonzero only
        # | 3 both kinds
        multi = coexisting_stable_nonzero(p, sc)
        c = (1 if multi else 0) + (2 if sc.bistable else 0)
        if multi and sc.bistable:
            c = 3
        elif sc.bistable:
            c = 2
        code[i, j] = c
        if sc.folds:
            kfold_map[i, j] = max(f[0] for f in sc.folds)
        rows.append({"eta": eta, "rho": rho, "class": c, "direction": sc.direction, "n_folds": len(sc.folds),
                     "tumour_free_plus_nonzero_bistable": sc.bistable, "two_stable_nonzero": multi, "k_prime0": sc.k_prime0,
                     "bistable_k_hi": sc.bistable_k_range[1] if sc.bistable_k_range else np.nan})
pd.DataFrame(rows).to_csv(TAB / "bistability_map_eta_rho.csv", index=False)
log(f"   (eta,rho) map: degenerate eta* at baseline = {degenerate_eta(p0):.4f}, rho* = {degenerate_rho(p0):.5f}")

# =========================================================================== 5. backward example + simulations
log("5. backward-bifurcation example (eta=0.6): continuation + simulation of bistability")
pB = p0.with_(eta=0.6)
kcB = critical_k(pB)
brB = trace_all_branches(pB, (0.0, 1.2), n_seed=31)
brB_df, evB = branch_frame(brB, pB, "eta=0.6"), event_frame(brB, "eta=0.6")
brB_df.to_csv(TAB / "bifurcation_branches_eta0p6.csv", index=False)
evB.to_csv(TAB / "bifurcation_events_eta0p6.csv", index=False)
fold_k = float(branch_k(fold_P(pB)[0], pB))
print(evB[evB.kind != "neutral saddle"][["branch", "kind", "k", "P"]].to_string())
rows = []
for kk in [0.5 * (kcB + fold_k), 0.9 * fold_k, 0.5 * kcB]:
    p = pB.with_(k=kk)
    summ = attractor_summary(p)
    for name, y0 in [("small residual disease", [1e-6, 1e-4, p.E_free]), ("large burden", [0.3, 0.3, 0.1])]:
        sol = solve_ivp(rhs, (0, 4e5), y0, args=(p,), method="Radau", rtol=1e-10, atol=1e-14)
        yf = sol.y[:, -1]
        rows.append({"k": kk, "k_over_kc": kk / kcB, "k_over_kfold": kk / fold_k, "n_stable_equilibria_incl_TF": summ["n_stable"],
                     "initial_condition": name, "y0": str(y0), "final_P": yf[0], "final_D": yf[1], "final_E": yf[2],
                     "final_residual": float(np.max(np.abs(rhs(0, yf, p)))),
                     "attractor": "tumour-free" if yf[0] + yf[1] < 1e-6 else "nonzero equilibrium"})
pd.DataFrame(rows).to_csv(TAB / "bistability_simulations.csv", index=False)
print(pd.DataFrame(rows)[["k_over_kc", "n_stable_equilibria_incl_TF", "initial_condition", "final_P", "attractor"]].to_string())

# =========================================================================== 6. search for sustained oscillations
log("6. direct-simulation search for sustained oscillations / limit cycles")
rows = []
N_OSC = 300
osc = {"a": loguni(N_OSC, 1e-4, 1.0), "q": loguni(N_OSC, 1e-3, 1.0), "rho": loguni(N_OSC, 1e-3, 2.0),
       "eta": loguni(N_OSC, 1e-3, 5.0), "r": loguni(N_OSC, 0.02, 0.5), "delta_D": loguni(N_OSC, 1e-4, 0.1),
       "mu": loguni(N_OSC, 0.01, 0.5), "h": loguni(N_OSC, 1e-3, 0.5)}
n_sustained = n_run = n_flagged = 0
for i in range(N_OSC):
    p = Params(**{k: float(v[i]) for k, v in osc.items()})
    if p.r_eff <= 0:
        continue
    kk = float(rng.uniform(0.0, 2.0) * critical_k(p))
    p = p.with_(k=kk)
    y0 = [float(10 ** rng.uniform(-5, -0.5)), float(10 ** rng.uniform(-5, -0.5)), float(p.E_free * 10 ** rng.uniform(-1, 0.7))]
    T = 1.5e4
    try:
        # states are clipped to the nonnegative orthant inside the RHS only to keep the solver away from the
        # pole of P/(h+P) at P = -h; the true flow leaves the orthant invariant, so this changes nothing physical
        sol = solve_ivp(lambda t, y, p: rhs(t, np.maximum(y, 0.0), p), (0, T), y0, args=(p,), method="LSODA",
                        rtol=1e-8, atol=1e-12, t_eval=np.linspace(0.6 * T, T, 600))
    except Exception:
        continue
    if not sol.success:
        continue
    n_run += 1
    N = sol.y[0] + sol.y[1]
    half = len(N) // 2
    amp1, amp2 = np.ptp(N[:half]), np.ptp(N[half:])
    rel = np.ptp(N) / max(np.max(N), 1e-12)
    flagged = bool(rel > 1e-3 and amp2 > 0.5 * amp1 and np.max(N) > 1e-6)     # stage 1: crude, over-inclusive
    if not flagged:
        continue
    # stage 2: a slow monotone relaxation is not an oscillation. Re-integrate very long (T = 3e6 months) and
    # require a non-vanishing residual |F| or many extrema of N of non-negligible relative amplitude.
    n_flagged += 1
    sol2 = solve_ivp(lambda t, y, p: rhs(t, np.maximum(y, 0.0), p), (0, 3e6), y0, args=(p,), method="Radau",
                     rtol=1e-10, atol=1e-14, dense_output=True)
    yf = sol2.y[:, -1]
    res_end = float(np.max(np.abs(rhs(0, yf, p))))
    tt = np.linspace(2e6, 3e6, 4000)
    Nl = sol2.sol(tt)[0] + sol2.sol(tt)[1]
    dN = np.diff(Nl)
    n_ext = int(np.sum(dN[1:] * dN[:-1] < 0))
    rel2 = float(np.ptp(Nl) / max(np.max(Nl), 1e-14))
    persistent = bool(res_end > 1e-10 or (n_ext >= 4 and rel2 > 1e-6 and np.max(Nl) > 1e-6))
    n_sustained += persistent
    rows.append({**{k: float(v[i]) for k, v in osc.items()}, "k": kk, "y0": str(y0), "stage1_relative_range_N": rel,
                 "residual_after_3e6_months": res_end, "n_extrema_late": n_ext, "late_relative_range_N": rel2,
                 "persistent_oscillation": persistent})
osc_summary = {"n_simulated": n_run, "n_flagged_stage1": n_flagged, "n_persistent_oscillation_after_long_integration": n_sustained}
pd.DataFrame([osc_summary]).to_csv(TAB / "oscillation_search_summary.csv", index=False)

# =========================================================================== FIGURES
log("7. figures")
note = baseline_note(p0)


def seg_plot(ax, k, y, lab, phys, **kw):
    """Plot a branch with solid=stable, dashed=unstable, dotted grey=nonphysical (P<0).

    Contiguous points of equal style are drawn as ONE line so that dashes are visible; marginal points
    (the located bifurcation points themselves) inherit the style of the neighbouring segment.
    """
    k, y, lab, phys = np.asarray(k), np.asarray(y), np.asarray(lab), np.asarray(phys)
    color = kw.get("color", BLUE)

    def style(i, j):
        if not (phys[i] and phys[j]):
            return "nonphys"
        for m in (i, j):
            if lab[m] != "marginal":
                return "stable" if lab[m].startswith("stable") else "unstable"
        return "stable"

    n = len(k)
    start = 0
    while start < n - 1:
        st = style(start, start + 1)
        end = start + 1
        while end < n - 1 and style(end, end + 1) == st:
            end += 1
        sl = slice(start, end + 1)
        if st == "nonphys":
            ax.plot(k[sl], y[sl], color=GREY, ls=":", lw=1.3)
        else:
            ax.plot(k[sl], y[sl], color=color, ls="-" if st == "stable" else (0, (5, 3)), lw=1.9)
        start = end


tf = df0[(df0.P.abs() < 1e-14)].sort_values("k")
nz = df0[(df0.P.abs() >= 1e-14)].sort_values("k")

# ---- A: branches
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
for ax, ycol, ylab in [(axes[0], "N", "equilibrium burden $N^*=P^*+D^*$ [fraction of K]"),
                       (axes[1], "E", "equilibrium immune level $E^*$ [E-units]")]:
    ax.plot(tf.k, tf[ycol] if ycol == "E" else tf.N, color=BLUE, lw=2.2, label="tumour-free branch")
    kk = np.linspace(0, kc * 1.12, 20)
    Pg = np.linspace(-0.45 * p0.h, p0.r_eff * p0.K / (p0.r * p0.phi) * 0.999, 4000)
    Bg = p0.mu + p0.eta * p0.phi * Pg - p0.rho * Pg / (p0.h + Pg)
    kg = branch_k(Pg, p0)
    ok = (Bg > 0)
    yy = (Pg * p0.phi if ycol == "N" else p0.s / Bg)
    phys = ok & (Pg >= 0) & (kg >= 0)
    ax.plot(kg[phys], yy[phys], color=ORANGE, lw=2.2, label="nonzero branch (admissible, $P>0$)")
    ext = ok & (Pg < 0) & (kg < KMAX)
    ax.plot(kg[ext], yy[ext], color=GREY, ls=":", lw=1.4, label="algebraic continuation, $P<0$ (nonphysical)")
    ax.axvline(kc, color=INK2, lw=0.8, ls="-.")
    ax.text(kc + 0.006, 0.55 if ycol == "N" else 0.58, f"$k_c$={kc:.4f}", fontsize=8, color=INK2)
    ax.set_xlim(0, kc * 1.6)
    ax.set_xlabel(K_LABEL); ax.set_ylabel(ylab)
axes[0].set_ylim(-0.03, 0.75); axes[1].set_ylim(0.1, 0.65)
axes[0].legend(loc="upper right")
axes[0].set_title("A. Equilibrium branches vs immune killing")
axes[1].set_title("immune level on the same branches")
fig.tight_layout(rect=(0, 0.04, 1, 1)); footer(fig, note + ". " + DISCLAIMER)
fig.savefig(FIG / "fig_A_equilibrium_branches.png"); plt.close(fig)

# ---- B: stability-classified diagram, baseline + backward example
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
ax = axes[0]
seg_plot(ax, tf.k, tf.N, tf.stability, tf.physical, color=BLUE)
for bi, b in nz.groupby("branch"):
    b = b.sort_values("k") if False else b            # keep arclength order
    seg_plot(ax, b.k.values, b.N.values, b.stability.values, b.physical.values, color=ORANGE)
ax.plot([kc], [0], "o", ms=8, mfc="white", mec=INK, mew=1.5, zorder=5)
ax.annotate("transcritical\n$k_c$ (zero eigenvalue)", (kc, 0), xytext=(kc + 0.03, 0.18), arrowprops=dict(arrowstyle="-", color=INK2), fontsize=8)
ax.set_xlim(0, kc * 1.5); ax.set_ylim(-0.04, 0.7)
ax.set_xlabel(K_LABEL); ax.set_ylabel("equilibrium burden $N^*$ [fraction of K]")
ax.set_title("B1. Baseline: single transcritical point, no other bifurcation")
from matplotlib.lines import Line2D
leg = [Line2D([], [], color=BLUE, lw=2, label="tumour-free"), Line2D([], [], color=ORANGE, lw=2, label="nonzero"),
       Line2D([], [], color=INK, lw=2, ls="-", label="solid: stable"), Line2D([], [], color=INK, lw=2, ls="--", label="dashed: unstable"),
       Line2D([], [], color=GREY, lw=1.2, ls=":", label="dotted: nonphysical ($P<0$)")]
ax.legend(handles=leg, loc="upper right")
ax = axes[1]
tfB = brB_df[brB_df.P.abs() < 1e-14].sort_values("k")
seg_plot(ax, tfB.k, tfB.N, tfB.stability, tfB.physical, color=BLUE)
for bi, b in brB_df[brB_df.P.abs() >= 1e-14].groupby("branch"):
    seg_plot(ax, b.k.values, b.N.values, b.stability.values, b.physical.values, color=ORANGE)
evp = evB[evB.kind.isin(["fold", "branch point (P=0, zero eigenvalue)"])]
ev_phys = evp[evp.P >= -1e-12]
for _, e in ev_phys.iterrows():
    ax.plot([e.k], [(e.P + e.D)], "o" if e.kind == "fold" else "s", ms=7, mfc="white", mec=INK, mew=1.5, zorder=5)
ax.axvspan(kcB, fold_k, color=VIOLET, alpha=0.10)
ax.annotate("bistable window\n(tumour-free and nonzero\nboth stable)", (0.5 * (kcB + fold_k), 0.10), ha="center", fontsize=8, color=INK2)
ax.annotate("fold", (fold_k, 0.058 * pB.phi), xytext=(fold_k + 0.04, 0.55), fontsize=8, arrowprops=dict(arrowstyle="-", color=INK2))
ax.annotate("transcritical", (kcB, 0), xytext=(kcB - 0.05, 0.25), fontsize=8, arrowprops=dict(arrowstyle="-", color=INK2))
ax.set_xlim(0, 1.15); ax.set_ylim(-0.04, 0.8)
ax.set_xlabel(K_LABEL); ax.set_ylabel("equilibrium burden $N^*$ [fraction of K]")
ax.set_title(r"B2. Verified example with $\eta=0.6$: backward bifurcation + fold")
fig.tight_layout(rect=(0, 0.04, 1, 1)); footer(fig, "Baseline except B2: eta=0.6 (1/(month*burden)). Circles: fold; squares: transcritical. Continuation points; " + DISCLAIMER)
fig.savefig(FIG / "fig_B_stability_bifurcation_diagram.png"); plt.close(fig)

# ---- C: largest real eigenvalue (global + zoom at k_c)
fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.1))
for ax, zoom in zip(axes, (False, True)):
    kfine = np.linspace(0.0, 1.6 * kc, 3000)      # exact eigenvalues of the tumour-free Jacobian on a fine grid
    ax.plot(kfine, [max_real_eigenvalue(jacobian(np.array([0, 0, p0.E_free]), p0.with_(k=kk))) for kk in kfine],
            color=BLUE, label="tumour-free")
    nzp = nz[nz.physical]
    for bi, b in nzp.groupby("branch"):
        ax.plot(b.k, b.max_re_eig, color=ORANGE, label="nonzero (admissible)" if bi == nzp.branch.min() else None)
    nzn = nz[~nz.physical]
    for bi, b in nzn.groupby("branch"):
        ax.plot(b.k, b.max_re_eig, color=GREY, ls=":", label="nonzero, $P<0$ (nonphysical)" if bi == nzn.branch.min() else None)
    ax.axhline(0, color=INK, lw=0.8); ax.axvline(kc, color=INK2, lw=0.8, ls="-.")
    ax.set_xlabel(K_LABEL)
    ax.set_ylabel(r"largest $\mathrm{Re}\,\lambda$ of the Jacobian [1/month]")
    if zoom:
        ax.set_xlim(0.9 * kc, 1.1 * kc); ax.set_ylim(-0.0025, 0.0025)
        ax.set_title("C2. Zoom at $k_c$: $\\lambda_{TF}\\approx\\alpha(k-k_c)$, $\\lambda_{nz}\\approx-\\alpha(k-k_c)$")
    else:
        ax.set_xlim(0, kc * 1.6); ax.set_ylim(-0.12, 0.12)
        ax.set_title("C1. Stability exchange at $k_c$")
        ax.legend(loc="upper right")
fig.tight_layout(rect=(0, 0.05, 1, 1)); footer(fig, note + ". " + DISCLAIMER)
fig.savefig(FIG / "fig_C_max_real_eigenvalue.png"); plt.close(fig)

# ---- D: dormant fraction
fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9))
ax = axes[0]
ax.plot(nzp.k, nzp.dormant_fraction, color=ORANGE, lw=2.2, label="nonzero equilibria (computed)")
ax.axhline(p0.q / (p0.q + p0.a + p0.delta_D), color=INK2, ls="--", lw=1, label=r"identity $q/(q+a+\delta_D)$=" + f"{p0.q/(p0.q+p0.a+p0.delta_D):.4f}")
ax.set_ylim(0.8, 0.95); ax.set_xlim(0, kc * 1.1)
ax.set_xlabel(K_LABEL); ax.set_ylabel("equilibrium dormant fraction $D^*/N^*$")
ax.set_title("D1. Equilibrium dormant fraction is a constant of the model, not an outcome")
ax.legend(loc="lower left")
ax = axes[1]
kgrid = np.linspace(0.02, 0.9, 30) * kc
fr = []
for kk in kgrid:
    s = simulate(p0.with_(k=kk), [1e-6, 1e-4, p0.E_free], 0.1, 2e4)
    fr.append(s.D[-1] / s.N[-1] if s.event else np.nan)
ax.plot(kgrid, fr, "o-", color=VIOLET, ms=4, label=r"$D/N$ at first crossing of $N_{det}=0.1$")
ax.axhline(p0.q / (p0.q + p0.a + p0.delta_D), color=INK2, ls="--", lw=1)
ax.axhline(1e-4 / (1e-4 + 1e-6), color=GREY, ls=":", lw=1, label="initial $D_0/N_0$=0.990")
ax.set_ylim(0.0, 1.02)
ax.set_xlabel(K_LABEL); ax.set_ylabel("dormant fraction at relapse time")
ax.set_title("D2. Dormant fraction along trajectories at detection (transient)")
ax.legend(loc="center left")
fig.tight_layout(rect=(0, 0.05, 1, 1)); footer(fig, note + f"; P0=1e-6, D0=1e-4, E0=s/mu. Missing points: no crossing within 2e4 months. " + DISCLAIMER)
fig.savefig(FIG / "fig_D_dormant_fraction.png"); plt.close(fig)

# ---- E: phase portraits near the critical boundary
fig, axes = plt.subplots(1, 3, figsize=(13, 4.3), sharey=True)
c_eq = p0.q / (p0.a + p0.delta_D)
for ax, f in zip(axes, [0.9, 1.0, 1.1]):
    p = p0.with_(k=f * kc)
    for N0 in [0.002, 0.01, 0.025, 0.045, 0.065]:
        for dfrac, col in [(c_eq / (1 + c_eq), BLUE), (0.0, ORANGE)]:
            for E0 in [0.34, 0.40, 0.47]:
                y0 = [N0 * (1 - dfrac), N0 * dfrac, E0]
                sol = solve_ivp(rhs, (0, 3000), y0, args=(p,), method="LSODA", rtol=1e-9, atol=1e-13, dense_output=True)
                yy = sol.sol(np.linspace(0, 3000, 1500))
                ax.plot(yy[0] + yy[1], yy[2], color=col, lw=0.7, alpha=0.6)
                ax.plot(y0[0] + y0[1], y0[2], ".", color=col, ms=3)
    for eq in enumerate_equilibria(p):
        lab = stability_label(np.linalg.eigvals(jacobian(eq.state, p)))
        ax.plot([eq.P + eq.D], [eq.E], "o", ms=9, mfc=ORANGE if lab.startswith("stable") else "white", mec=INK, mew=1.5, zorder=6, clip_on=False)
    tf_stable = f > 1
    ax.plot([0], [p0.E_free], "s", ms=9, mfc=BLUE if tf_stable else "white", mec=INK, mew=1.5, zorder=6, clip_on=False)
    tlab = {0.9: "tumour-free: saddle (unstable)", 1.0: "tumour-free: zero eigenvalue (marginal)", 1.1: "tumour-free: stable"}[f]
    ax.set_title(f"$k={f:.1f}\\,k_c$ = {f*kc:.4f}\n{tlab}")
    ax.set_xlim(-0.003, 0.075); ax.set_ylim(0.32, 0.5)
    ax.set_xlabel("burden $N=P+D$ [fraction of K]")
axes[0].set_ylabel("immune level $E$ [E-units]")
axes[0].text(0.03, 0.03, "orange: initial D/N=0;  blue: initial D/N=0.87\nsquare: tumour-free; circle: nonzero\nfilled = stable, open = unstable; dots = initial states",
             transform=axes[0].transAxes, fontsize=7, color=INK2, va="bottom")
fig.suptitle("E. Phase portraits (projection of the 3-D flow onto $(N,E)$) across the transcritical boundary", y=1.0, fontsize=10)
fig.tight_layout(rect=(0, 0.05, 1, 0.97)); footer(fig, note + "; trajectories to t=3000 months. Projection: crossings are not intersections in 3-D. " + DISCLAIMER)
fig.savefig(FIG / "fig_E_phase_portraits.png"); plt.close(fig)

# ---- F: additional (verified) bifurcation structure
fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.4))
ax = axes[0]
etas_c = np.linspace(0.0, 1.5, 76)
kfold_arr, bistab_hi, nn_lo, nn_hi = [], [], [], []
for eta in etas_c:
    p = p0.with_(eta=eta)
    fp = fold_P(p)
    kfold_arr.append(max([branch_k(P, p) for P in fp]) if len(fp) else np.nan)
    sc = scan_branch(p, n=801)
    bistab_hi.append(sc.bistable_k_range[1] if sc.bistable_k_range else np.nan)
    w = [w_ for w_ in stable_overlap_windows(sc) if w_[1] - w_[0] > 1e-6]
    nn_lo.append(w[0][0] if w else np.nan); nn_hi.append(w[0][1] if w else np.nan)
kfold_arr, bistab_hi = np.array(kfold_arr), np.array(bistab_hi)
nn_lo, nn_hi = np.array(nn_lo), np.array(nn_hi)
ax.plot(etas_c, kfold_arr, color=ORANGE, label="largest saddle-node $k_{SN}$")
ax.axhline(kc, color=BLUE, label="$k_c$ (transcritical; independent of $\\eta$)")
eta_star = degenerate_eta(p0)
ax.axvline(eta_star, color=INK2, ls="-.", lw=0.9)
ax.fill_between(etas_c, kc, np.where(np.isnan(bistab_hi), kc, bistab_hi), color=VIOLET, alpha=0.18,
                label="tumour-free + nonzero both stable")
ax.fill_between(etas_c, np.where(np.isnan(nn_lo), kc, nn_lo), np.where(np.isnan(nn_hi), kc, nn_hi), color=AQUA, alpha=0.45,
                label="two stable nonzero states (first window)")
ax.text(eta_star + 0.03, 0.05, f"$\\eta^*$={eta_star:.3f}: $k'(0)=0$\nbackward for $\\eta>\\eta^*$", fontsize=8, color=INK2)
ax.set_ylim(0, 2.2)
ax.set_xlabel(r"immune suppression $\eta$ [1/(month$\cdot$burden)]"); ax.set_ylabel(K_LABEL)
ax.set_title("F1. Saddle-node and coexistence windows vs $\\eta$"); ax.legend(loc="upper left", fontsize=7)
ax = axes[1]
cmap = matplotlib.colors.ListedColormap(["#e4e3df", "#9fd9c3", "#c9c4f0", "#6250d6"])
ax.imshow(code, origin="lower", extent=[etas[0], etas[-1], rhos[0], rhos[-1]], aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5)
rho_line = np.array([degenerate_rho(p0.with_(eta=e)) for e in np.linspace(0, 2, 200)])
ax.plot(np.linspace(0, 2, 200), rho_line, color=INK, lw=1.4)
ax.plot([p0.eta], [p0.rho], "*", color=ORANGE, ms=14, mec=INK)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color="#e4e3df", label="no coexistence of stable states"),
                   Patch(color="#9fd9c3", label="two stable NONZERO states (k<k_c; tumour-free unstable)"),
                   Patch(color="#c9c4f0", label="tumour-free + one nonzero stable"),
                   Patch(color="#6250d6", label="both kinds of coexistence"),
                   Line2D([], [], color=INK, lw=1.4, label="$k'(0)=0$: $\\rho^*(\\eta)$; backward (subcritical) below/right"),
                   Line2D([], [], marker="*", color="w", mfc=ORANGE, mec=INK, ms=11, label="baseline")],
          loc="lower right", fontsize=6, frameon=True, facecolor="white", framealpha=0.95)
ax.set_ylim(0, 0.2); ax.set_xlim(0, 2); ax.grid(False)
ax.set_xlabel(r"$\eta$ [1/(month$\cdot$burden)]"); ax.set_ylabel(r"$\rho$ [1/month]")
ax.set_title("F2. Regimes in $(\\eta,\\rho)$, other parameters at baseline")
ax = axes[2]
Pg = np.linspace(1e-6, 0.12, 400)
for eta in [0.0, 0.2, 0.4, 0.6, 1.0]:
    p = p0.with_(eta=eta)
    Bg = p.mu + p.eta * p.phi * Pg - p.rho * Pg / (p.h + Pg)
    kg = branch_k(Pg, p); ok = (Bg > 0) & (kg > 0)
    ax.plot(kg[ok], (Pg * p.phi)[ok], label=f"$\\eta$={eta}")
ax.axvline(kc, color=INK2, ls="-.", lw=0.8)
ax.set_xlim(0, 1.3); ax.set_ylim(0, 0.9)
ax.set_xlabel(K_LABEL); ax.set_ylabel("nonzero-branch burden $N^*$")
ax.set_title("F3. Nonzero branch $k(P)$ for several $\\eta$ (curves only)")
ax.legend(loc="upper right")
fig.tight_layout(rect=(0, 0.05, 1, 1))
footer(fig, "Other parameters at baseline. No Hopf point was found in the scanned ranges (parameter_scan_summary.csv); none is drawn. Stability along the curves: see Fig. B2. " + DISCLAIMER)
fig.savefig(FIG / "fig_F_additional_bifurcations.png"); plt.close(fig)

log("done")
