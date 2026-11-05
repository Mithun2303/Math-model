"""Phase 4, Tasks 1 and 3: structural identifiability (exact), event-only identifiability, Fisher spectra, profile likelihoods.

Run: .venv/bin/python experiments/identifiability_analysis.py
Outputs: results/tables/p4_structural*.csv, p4_event_only*.csv, p4_fisher*.csv, p4_profile*.csv; figures fig_p4_02, 03, 04, 11.
Structural results are exact statements about the model (generic LOCAL identifiability, see src/dormancy/identifiability.py);
Fisher spectra and profile likelihoods are NUMERICAL, practical statements about SYNTHETIC data under assumed noise.
"""
import math

from _p4common import *  # noqa: F401,F403
import matplotlib.pyplot as plt

from dormancy.identifiability import (PDEModel, OBS_SCENARIOS, event_time_gradient, linear_regime_symbolic, null_directions,
                                      numerical_spectrum, structural_report)
from dormancy.observation import ObservationAssumptions, generate_dataset, simulate_hidden
from dormancy.synthetic_validation import (FREE_SETS, FitProblem, MODEL_FREE, make_bounds, multistart, profile_likelihood,
                                           residuals_jacobian, wald_and_spectrum)
from dormancy.uncertainty import parallel_map

LF = CFG["fitting"]["bounds_log_factor"]
PARAMS_FIT = MODEL_FREE + ["N0", "fD"]
LAB = {"r": "r", "q": "q", "a": "a", "delta_D": "$\\delta_D$", "beta": "$\\beta$", "rho": "$\\rho$", "eta": "$\\eta$", "N0": "$N_0$", "fD": "$f_D$"}


def profile_worker(args):
    data, free, fixed, lo, hi, sigma, xbest, cost, name, hw, npts, seed = args
    from dormancy.synthetic_validation import FitResult
    prob = FitProblem(data, free, fixed, lo, hi, sigma)
    best = FitResult(xbest, cost, True, 0, np.zeros(len(free), bool), xbest)
    pr = profile_likelihood(prob, best, name, hw, npts, seed=seed)
    return pr


def main():
    seeds = CFG["seeds"]
    # ----------------------------------------------------------------------------- A. structural identifiability
    unknown_sets = {
        "A: 7 parameters (r,q,a,dD,beta,rho,eta); initial state known": ["r", "q", "a", "delta_D", "beta", "rho", "eta"],
        "B: A + P0,D0": ["r", "q", "a", "delta_D", "beta", "rho", "eta", "P0", "D0"],
        "C: B + h,mu": ["r", "q", "a", "delta_D", "beta", "rho", "eta", "P0", "D0", "h", "mu"],
        "D: C + e0": ["r", "q", "a", "delta_D", "beta", "rho", "eta", "P0", "D0", "h", "mu", "e0"],
        "E: D + signal scale c (K = 1 fixed)": ["r", "q", "a", "delta_D", "beta", "rho", "eta", "P0", "D0", "h", "mu", "e0", "c"],
        "F: E + K (burden scale unknown)": ["r", "q", "a", "delta_D", "beta", "rho", "eta", "P0", "D0", "h", "mu", "e0", "c", "K"],
    }
    rows = []
    for nm, U in unknown_sets.items():
        for sc in OBS_SCENARIOS:
            UU = list(U) + (["ce"] if sc in ("Ne", "PDe") and "c" in U else [])
            rep = structural_report(UU, sc, max_order=30, n_points=3, seed=CFG["identifiability"]["seed"])
            rows.append({"unknown_set": nm, "observed": sc, "n_unknown": rep["n_unknown"], "rank": rep["rank"], "deficiency": rep["deficiency"],
                         "ranks_at_3_random_points": str(rep["ranks_per_point"]), "rank_flat_over_last_6_orders": rep["flat_last_orders"],
                         "generic_local_structural_identifiability": rep["structurally_locally_identifiable"],
                         "rank_vs_order_1_5_10_15_20_25_30": ",".join(str(rep["rank_vs_order"][k]) for k in (1, 5, 10, 15, 20, 25, 30))})
    st = pd.DataFrame(rows)
    save_table(st, "p4_structural_identifiability")
    print(st[["unknown_set", "observed", "n_unknown", "rank", "deficiency"]].to_string())
    # exact null direction of the scale gauge (F, N-only)
    U = unknown_sets["F: E + K (burden scale unknown)"]
    pt, vecs, rk = null_directions(U, "N", order=18)
    gauge = pd.DataFrame([{"unknown": u, "null_vector_component": float(v), "point_value": float(pt[u]) if u in pt else float("nan"),
                           "component_over_value": float(v / pt[u]) if (u in pt and v != 0) else float("nan")} for u, v in zip(U, vecs[0]) if v != 0]) if vecs else pd.DataFrame()
    save_table(gauge, "p4_structural_scale_gauge_null_vector")
    print(gauge)

    # sub-models: which individual parameters are identifiable (null-space components) -- linear regime and frozen immunity
    sub = []
    for mname, model in {"full": PDEModel(), "logistic only (immune frozen)": PDEModel(logistic=True, immune=False),
                         "linear (no logistic, immune frozen)": PDEModel(logistic=False, immune=False)}.items():
        for sc in ("N", "PD"):
            Uu = ["r", "q", "a", "delta_D", "beta", "P0", "D0"] + (["rho", "eta"] if mname == "full" else [])
            pt, vecs, rk = null_directions(Uu, sc, model=model, order=18)
            nonid = [u for u in Uu if any(v[Uu.index(u)] != 0 for v in vecs)]
            sub.append({"model": mname, "observed": sc, "unknowns": ",".join(Uu), "rank": rk, "n_unknown": len(Uu), "null_dim": len(vecs),
                        "parameters_with_nonzero_null_component (not individually identifiable)": ",".join(nonid) if vecs else "none"})
    sub = pd.DataFrame(sub)
    save_table(sub, "p4_structural_submodels")
    print(sub.to_string())
    lin = linear_regime_symbolic()
    save_table(pd.DataFrame([{"quantity": "symbolic generating-series rank, linear regime, N-only, unknowns (m,a,q,s,P0,D0)", "value": lin["rank"]},
                             {"quantity": "null dimension", "value": lin["null_dim"]},
                             {"quantity": "interpretation", "value": "N(t)=c+ e^{lam+ t}+c- e^{lam- t}: 4 numbers for 6 unknowns; r and beta enter only as m=r-q-beta*e0"}]), "p4_structural_linear_regime")

    # ----------------------------------------------------------------------------- B. event-only (scenario 5)
    th = truth("G1_baseline"); th_e = {**th}
    T0, grad = event_time_gradient(th_e, 0.1)
    names = list(grad)
    g = np.array([grad[k] for k in names])
    from scipy.integrate import solve_ivp
    ev_rows = [{"quantity": "T (first crossing of 0.1), G1 [months]", "value": T0}]
    for k in names:
        ev_rows.append({"quantity": f"dT/dln({k}) [months]", "value": grad[k]})
    ev_rows.append({"quantity": "rank of the 1 x p event-time Jacobian", "value": 1})
    ev_rows.append({"quantity": "dimension of the iso-T manifold (p - 1) for p = 10 inputs", "value": len(names) - 1})
    save_table(pd.DataFrame(ev_rows), "p4_event_only_gradient")
    # iso-T demonstration: change beta by factors while re-solving r so that T stays equal; compare N(t) curves
    from dormancy.observation import detection_time
    iso = []
    for fb in (0.8, 0.9, 1.1, 1.25):
        th2 = {**th, "beta": th["beta"] * fb}
        lo_, hi_ = th["r"] * 0.3, th["r"] * 3.0
        f = lambda r_: detection_time({**th2, "r": r_}, 0.1, 5000.0) - T0
        from scipy.optimize import brentq
        try:
            r_star = brentq(f, lo_, hi_, xtol=1e-12)
        except ValueError:
            r_star = float("nan")
        if np.isfinite(r_star):
            th3 = {**th2, "r": r_star}
            tt = np.array([0.0, 60.0, 120.0, 180.0, T0 * 0.99])
            N1 = simulate_hidden(th, tt).N; N2 = simulate_hidden(th3, tt).N
            iso.append({"beta_factor": fb, "r_factor_needed_for_equal_T": r_star / th["r"], "T_equal": detection_time(th3, 0.1, 5000.0),
                        **{f"N(t={t:g}) ratio new/original": float(N2[i] / N1[i]) for i, t in enumerate(tt)}})
    iso = pd.DataFrame(iso)
    save_table(iso, "p4_event_only_iso_T_family")
    print(iso.round(4).to_string())
    # linear-regime confounding: scaling N0 shifts T by -ln(s)/lambda_plus (time-origin / threshold equivalence)
    from dormancy.relapse import unstable_mode
    p_, y0_ = None, None
    from dormancy.observation import theta_to_model
    p_, y0_ = theta_to_model(th)
    lam = unstable_mode(p_).lam
    conf = []
    for s_ in (0.1, 0.3, 3.0, 10.0):
        T1 = detection_time({**th, "N0": th["N0"] * s_}, 0.1, 5000.0)
        conf.append({"N0_scale_s": s_, "T_numeric": T1, "dT_numeric": T1 - T0, "dT_linear_prediction (-ln s / lambda_plus)": -math.log(s_) / lam,
                     "lambda_plus": lam})
        # threshold scaled instead (N0 fixed): T(N_d) with N_d = 0.1*s: dT = +ln(s)/lambda in the linear regime
        T2 = detection_time(th, 0.1 * s_ if 0.1 * s_ > th["N0"] else 0.11, 5000.0)
    conf = pd.DataFrame(conf)
    save_table(conf, "p4_event_only_N0_time_origin_confounding")
    print(conf.round(3).to_string())

    # ----------------------------------------------------------------------------- C. Fisher spectra
    fr = []
    dirs = []
    for obs in ("N", "PD", "Ne", "PDe"):
        for sig in (0.05, 0.10, 0.25):
            for trunc in (True, False):
                a = ObservationAssumptions(observed=obs, noise_sigma_ln=sig, lod=1e-7, truncate_at_detection=trunc)
                d = generate_dataset(th, a, seed=seeds["synthetic_data"])
                free = FREE_SETS["ic_unknown"]
                lo, hi = make_bounds({**th, "c": 1.0}, free, LF)
                prob = FitProblem(d, free, {"e0": th["e0"], "c": 1.0}, lo, hi)
                _, J = residuals_jacobian(prob, prob.to_x(th), True)
                sp = numerical_spectrum(J, 1.0)
                fr.append({"observed": obs, "sigma_ln": sig, "truncate_at_detection": trunc, "n_obs": len(d.times), "obs_end_month": float(d.times[-1]),
                           "n_dirs_s_gt_3": sp["n_above_3"], "n_dirs_s_gt_1": sp["n_above_1"],
                           **{f"s{i + 1}": float(v) for i, v in enumerate(sp["singular_values"])}})
                if sig == 0.10 and trunc:
                    for k in range(len(free)):
                        dirs.append({"observed": obs, "direction": k + 1, "singular_value": float(sp["singular_values"][k]),
                                     "sd_along_direction_ln": float(1 / sp["singular_values"][k]) if sp["singular_values"][k] > 0 else float("inf"),
                                     **{f"loading_{n}": float(c) for n, c in zip(free, sp["Vt"][k])}})
    fisher = pd.DataFrame(fr); save_table(fisher, "p4_fisher_spectra")
    save_table(pd.DataFrame(dirs), "p4_fisher_directions")
    print(fisher[["observed", "sigma_ln", "truncate_at_detection", "n_obs", "n_dirs_s_gt_3", "n_dirs_s_gt_1"]].to_string())

    # information on the r-beta separation vs initial burden (N-only, sigma 0.1)
    scal = []
    for s_ in (0.1, 1.0, 10.0, 100.0):
        th_s = {**th, "N0": th["N0"] * s_}
        a = ObservationAssumptions(observed="N", noise_sigma_ln=0.10, lod=1e-9, truncate_at_detection=True, detection_threshold=0.1)
        d = generate_dataset(th_s, a, seed=seeds["synthetic_data"])
        free = FREE_SETS["ic_known"]
        lo, hi = make_bounds({**th_s, "c": 1.0}, free, LF)
        prob = FitProblem(d, free, {"N0": th_s["N0"], "fD": th_s["fD"], "e0": 1.0, "c": 1.0}, lo, hi)
        _, J = residuals_jacobian(prob, prob.to_x(th_s), True)
        cov = np.linalg.pinv(J.T @ J, rcond=1e-14)
        i_r, i_b = free.index("r"), free.index("beta")
        # sd of the SEPARATION direction (ln r - ln beta) and of r - beta combination proxy
        v = np.zeros(len(free)); v[i_r], v[i_b] = 1, -1
        scal.append({"N0_scale": s_, "N0": th_s["N0"], "obs_end_month": float(d.times[-1]), "n_obs": len(d.times),
                     "sd_ln_r": float(np.sqrt(cov[i_r, i_r])), "sd_ln_beta": float(np.sqrt(cov[i_b, i_b])),
                     "sd_of_ln_r_minus_ln_beta": float(np.sqrt(v @ cov @ v)), "corr_r_beta": float(cov[i_r, i_b] / np.sqrt(cov[i_r, i_r] * cov[i_b, i_b]))})
    save_table(pd.DataFrame(scal), "p4_fisher_rbeta_vs_initial_burden")
    print(pd.DataFrame(scal).round(3).to_string())

    # ----------------------------------------------------------------------------- D. profile likelihoods (one data set per observed set)
    prof_rows, prof_store = [], {}
    for obs in ("N", "PD", "Ne", "PDe"):
        a = ObservationAssumptions(observed=obs, noise_sigma_ln=0.10, lod=1e-7)
        d = generate_dataset(th, a, seed=seeds["synthetic_data"] + 3)
        free = FREE_SETS["ic_unknown"]
        lo, hi = make_bounds({**th, "c": 1.0}, free, LF)
        fixed = {"e0": th["e0"], "c": 1.0}
        prob = FitProblem(d, free, fixed, lo, hi)
        ms = multistart(prob, 8, seeds["multistart"])
        best = ms[0]
        log(f"profiles: observed {obs}: best cost {best.cost:.2f}; n starts within 3.84: {sum(r.cost <= best.cost + 3.84 for r in ms)}/{len(ms)}")
        items = [(d, free, fixed, lo, hi, prob.sigma, best.x, best.cost, nm, math.log(10.0), CFG["fitting"]["profile_points"], seeds["profile"]) for nm in free]
        prs = parallel_map(profile_worker, items, chunksize=1)
        prof_store[obs] = (prob, best, prs)
        for pr in prs:
            lo_c = math.exp(pr["ci_lo"]) if pr["ci_lo"] is not None and pr["name"] != "fD" else (pr["ci_lo"] if pr["ci_lo"] is not None else float("nan"))
            prof_rows.append({"observed": obs, "parameter": pr["name"], "truth": th[pr["name"]],
                              "estimate": float(np.exp(best.x[free.index(pr["name"])])) if pr["name"] != "fD" else float(1 / (1 + np.exp(-best.x[free.index(pr["name"])]))),
                              "ci95_lo_log": pr["ci_lo"], "ci95_hi_log": pr["ci_hi"], "finite_both_sides": pr["finite_both"],
                              "open_low": pr["ci_lo"] is None, "open_high": pr["ci_hi"] is None,
                              "max_chi2_on_grid": float(np.nanmax(pr["chi2"])), "profile_window_half_width_ln": math.log(10.0)})
    pr_df = pd.DataFrame(prof_rows)
    save_table(pr_df, "p4_profile_likelihood_summary")
    print(pr_df[["observed", "parameter", "finite_both_sides", "open_low", "open_high", "max_chi2_on_grid"]].to_string())
    prof_curves = []
    for obs, (prob, best, prs) in prof_store.items():
        for pr in prs:
            for gx, c2 in zip(pr["grid"], pr["chi2"]):
                prof_curves.append({"observed": obs, "parameter": pr["name"], "grid_log": gx, "delta_chi2": c2})
    save_table(pd.DataFrame(prof_curves), "p4_profile_likelihood_curves")

    # ----------------------------------------------------------------------------- figures
    # Fig 3: profile likelihoods (observed N and observed P,D,e)
    fig, axes = plt.subplots(2, 9, figsize=(19, 5.6), sharey=True)
    for r_, obs in enumerate(("N", "PDe")):
        prob, best, prs = prof_store[obs]
        for c, pr in enumerate(prs):
            ax = axes[r_, c]
            x = pr["grid"] - best.x[prob.free.index(pr["name"])]
            ax.plot(x, pr["chi2"], "o-", ms=2.5, color=BLUE if obs == "N" else VIOLET)
            ax.axhline(3.841, color=MAGENTA, ls="--", lw=0.8)
            ax.axvline(0, color=GREY, lw=0.6)
            ax.set_ylim(-0.5, 12)
            ax.set_title(LAB[pr["name"]] + ("" if pr["finite_both"] else "  (open)"), fontsize=8)
            if r_ == 1:
                ax.set_xlabel("ln(theta/theta_hat)", fontsize=7)
            if c == 0:
                ax.set_ylabel(f"observed {obs}\n$\\Delta\\chi^2$", fontsize=8)
    fig.suptitle("Profile likelihoods (synthetic, G1 baseline; sigma_ln = 0.10; window +-ln 10; dashed: 95% threshold 3.84)", y=0.995, fontsize=10)
    fig.tight_layout(rect=(0, 0.06, 1, 0.95)); p4_footer(fig, "'(open)' = interval not closed within the profile window on at least one side: practically non-identifiable (Raue et al. 2009).")
    fig.savefig(FIG / "fig_p4_03_profile_likelihoods.png"); plt.close(fig)

    # Fig 2: parameter correlation heatmaps (Fisher information at the truth)
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    for ax, obs in zip(axes, ("N", "Ne", "PDe")):
        a = ObservationAssumptions(observed=obs, noise_sigma_ln=0.10, lod=1e-7)
        d = generate_dataset(th, a, seed=seeds["synthetic_data"])
        free = FREE_SETS["ic_unknown"]
        lo, hi = make_bounds({**th, "c": 1.0}, free, LF)
        prob = FitProblem(d, free, {"e0": th["e0"], "c": 1.0}, lo, hi)
        _, J = residuals_jacobian(prob, prob.to_x(th), True)
        cov = np.linalg.pinv(J.T @ J, rcond=1e-14)
        sd = np.sqrt(np.diag(cov)); corr = cov / np.outer(sd, sd)
        im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1)
        ax.set_xticks(range(len(free))); ax.set_xticklabels([LAB[n] for n in free], fontsize=8)
        ax.set_yticks(range(len(free))); ax.set_yticklabels([LAB[n] for n in free], fontsize=8)
        for i in range(len(free)):
            for j in range(len(free)):
                if i != j and abs(corr[i, j]) > 0.9:
                    ax.text(j, i, f"{corr[i, j]:.2f}", ha="center", va="center", fontsize=6)
        ax.set_title(f"observed {obs}: estimator correlation (Fisher, log-parameters)\nsd(ln): " + " ".join(f"{LAB[n].strip('$')}={s:.2g}" for n, s in zip(free, sd)), fontsize=7)
        ax.grid(False)
        fig.colorbar(im, ax=ax, fraction=0.046)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p4_footer(fig, "Linearised (Gauss-Newton) covariance at the truth, G1, sigma_ln = 0.10; |correlation| > 0.9 annotated; huge sd = non-identifiable direction.")
    fig.savefig(FIG / "fig_p4_02_parameter_correlation_heatmap.png"); plt.close(fig)

    # Fig 4: identifiability under observed-state combinations (singular spectra and estimable directions)
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    ax = axes[0]
    for obs, col in zip(("N", "PD", "Ne", "PDe"), (BLUE, AQUA, ORANGE, VIOLET)):
        s = fisher[(fisher.observed == obs) & (fisher.sigma_ln == 0.10) & fisher.truncate_at_detection].iloc[0]
        sv = [s[f"s{i + 1}"] for i in range(9)]
        ax.semilogy(range(1, 10), np.maximum(sv, 1e-6), "o-", color=col, label=f"observed {obs}")
    ax.axhline(3, color=MAGENTA, ls="--", lw=0.8); ax.text(8.9, 3.3, "s = 3", fontsize=7, ha="right", color=MAGENTA)
    ax.set_xlabel("singular value index"); ax.set_ylabel("singular value of the noise-scaled sensitivity matrix")
    ax.set_title("Fisher spectrum at the truth (9 log-parameters)", fontsize=9); ax.legend(fontsize=7)
    ax = axes[1]
    xs = np.arange(4); w = 0.26
    for k, (sig, col) in enumerate(((0.05, AQUA), (0.10, BLUE), (0.25, ORANGE))):
        vals = [fisher[(fisher.observed == o) & (fisher.sigma_ln == sig) & fisher.truncate_at_detection].n_dirs_s_gt_3.iloc[0] for o in ("N", "PD", "Ne", "PDe")]
        ax.bar(xs + (k - 1) * w, vals, w, color=col, label=f"sigma_ln {sig}")
    ax.axhline(9, color=INK2, ls=":"); ax.text(3.4, 9.1, "9 unknowns", fontsize=7, ha="right")
    ax.set_xticks(xs); ax.set_xticklabels(["N", "P, D", "N, e", "P, D, e"])
    ax.set_ylabel("directions with singular value > 3 (sd < 0.33 in ln)"); ax.set_title("Well-determined directions (practical)", fontsize=9); ax.legend(fontsize=7)
    ax = axes[2]
    for trunc, col, ls in ((True, BLUE, "-"), (False, ORANGE, "--")):
        vals = [fisher[(fisher.observed == o) & (fisher.sigma_ln == 0.10) & (fisher.truncate_at_detection == trunc)].n_dirs_s_gt_3.iloc[0] for o in ("N", "PD", "Ne", "PDe")]
        ax.plot(range(4), vals, "o" + ls, color=col, label="truncated at detection" if trunc else "NOT truncated (observe to 240 mo; unrealistic)")
    ax.set_xticks(range(4)); ax.set_xticklabels(["N", "P, D", "N, e", "P, D, e"]); ax.set_ylabel("directions with s > 3")
    ax.set_title("Effect of truncating observation at detection", fontsize=9); ax.legend(fontsize=7)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p4_footer(fig, "G1 baseline, 6-month visits, LOD 1e-7; structural rank is full (9/9 and 11/11) for every observed set: the shortfall is practical.")
    fig.savefig(FIG / "fig_p4_04_identifiability_by_observed_states.png"); plt.close(fig)

    # Fig 11: structural vs practical identifiability summary
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.8), gridspec_kw={"width_ratios": [1, 1.5]})
    ax = axes[0]
    sub_ = st[st.unknown_set.str.startswith(("A", "B", "D", "E", "F"))]
    piv = sub_.pivot(index="unknown_set", columns="observed", values="deficiency")[["N", "PD", "Ne", "PDe"]]
    ax.imshow((piv.values == 0).astype(float), cmap=matplotlib.colors.ListedColormap(["#f3c9b6", "#9fd9c3"]), vmin=0, vmax=1, aspect="auto")
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            ax.text(j, i, "identifiable" if piv.values[i, j] == 0 else f"deficiency {piv.values[i, j]}", ha="center", va="center", fontsize=7)
    ax.set_xticks(range(4)); ax.set_xticklabels(["N", "P, D", "N, e", "P, D, e"]); ax.set_yticks(range(piv.shape[0]))
    ax.set_yticklabels([s_.split(":")[0] + ": " + s_.split(":")[1][:34] for s_ in piv.index], fontsize=7); ax.grid(False)
    ax.set_title("STRUCTURAL (exact; generic local)", fontsize=9)
    ax = axes[1]
    status = np.zeros((9, 4))
    for j, obs in enumerate(("N", "PD", "Ne", "PDe")):
        for i, nm in enumerate(PARAMS_FIT):
            r_ = pr_df[(pr_df.observed == obs) & (pr_df.parameter == nm)].iloc[0]
            status[i, j] = 2 if r_.finite_both_sides else (1 if (not r_.open_low or not r_.open_high) else 0)
    cmap = matplotlib.colors.ListedColormap(["#f3c9b6", "#f6e3a1", "#9fd9c3"])
    ax.imshow(status, cmap=cmap, vmin=-0.5, vmax=2.5, aspect="auto")
    txt = {0: "not identifiable", 1: "one-sided", 2: "identifiable"}
    for i in range(9):
        for j in range(4):
            ax.text(j, i, txt[int(status[i, j])], ha="center", va="center", fontsize=7)
    ax.set_xticks(range(4)); ax.set_xticklabels(["N", "P, D", "N, e", "P, D, e"]); ax.set_yticks(range(9)); ax.set_yticklabels([LAB[n] for n in PARAMS_FIT]); ax.grid(False)
    ax.set_title("PRACTICAL (profile likelihood, G1, sigma_ln 0.10, one data set per column)", fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p4_footer(fig, "Structural identifiability does not imply practical identifiability: every column on the left is full rank, yet most parameters on the right are not identifiable from noisy truncated data.")
    fig.savefig(FIG / "fig_p4_11_structural_vs_practical_identifiability.png"); plt.close(fig)
    log("done")


import matplotlib  # noqa: E402

if __name__ == "__main__":
    main()
