"""Phase 4, Task 3: parameter redundancy and model reduction (nested reduced models fitted to synthetic data from the full model).

Run: .venv/bin/python experiments/model_reduction.py
The ORIGINAL full model stays the reference implementation; reduced models are special cases obtained by fixing parameters
(rho = eta = 0: immune frozen at e = 1; additionally delta_D = 0), and an alternative coordinate g = r - beta (net low-burden
growth after killing) in which the well-determined direction is explicit. Reduction is ACCEPTED only if it (i) fits full-model data
as well as the full model (AIC), (ii) has well-determined parameters (small Wald sd), (iii) predicts the crossing time equally well.
SYNTHETIC data, hypothetical observation model.
"""
from _p4common import *  # noqa: F401,F403
import matplotlib.pyplot as plt

from dormancy.observation import ObservationAssumptions, generate_dataset
from dormancy.synthetic_validation import (FitProblem, MODEL_FREE, make_bounds, multistart, prediction_errors, wald_and_spectrum)
from dormancy.uncertainty import parallel_map

MODELS = {
    "M9 full (r,q,a,dD,beta,rho,eta,N0,fD)": dict(free=["r", "q", "a", "delta_D", "beta", "rho", "eta", "N0", "fD"], fix={}),
    "M7 immune frozen (rho=eta=0)": dict(free=["r", "q", "a", "delta_D", "beta", "N0", "fD"], fix={"rho": 0.0, "eta": 0.0}),
    "M6 + delta_D = 0": dict(free=["r", "q", "a", "beta", "N0", "fD"], fix={"rho": 0.0, "eta": 0.0, "delta_D": 0.0}),
    "M5 + fD fixed at q/(q+a) (equilibrium ratio)": dict(free=["r", "q", "a", "beta", "N0"], fix={"rho": 0.0, "eta": 0.0, "delta_D": 0.0}),
}


def worker(args):
    gt_name, truth, obs, seed, mname = args
    spec = MODELS[mname]
    a = ObservationAssumptions(observed=obs, noise_sigma_ln=0.10, lod=1e-7)
    d = generate_dataset(truth, a, seed)
    free = spec["free"]
    fixed = {"e0": truth["e0"], "c": 1.0, **spec["fix"]}
    for k in MODEL_FREE + ["N0", "fD"]:
        if k not in free and k not in fixed:
            fixed[k] = truth[k]
    if mname.startswith("M5"):
        fixed["fD"] = truth["q"] / (truth["q"] + truth["a"])        # approx; refit q,a would change it, fixed at the TRUE ratio (generous)
    lo, hi = make_bounds({**truth, "c": 1.0}, free, 30.0)
    prob = FitProblem(d, free, fixed, lo, hi)
    ms = multistart(prob, 5, seed + 77)
    best = ms[0]
    w = wald_and_spectrum(prob, best)
    pe = prediction_errors(prob, best, {**truth, "c": 1.0})
    est = prob.to_theta(best.x)
    out = {"truth": gt_name, "observed": obs, "model": mname, "seed": seed, "k": len(free), "n_obs": len(d.times), "chi2": best.cost,
           "AIC": best.cost + 2 * len(free), "T_rel_error": pe["T_rel_error"], "rmse_lnN_in_window": pe["rmse_lnN_in_window"],
           "abs_lnN_error_at_240": pe["abs_lnN_error_at_240"], "n_bound_hits": int(best.bound_hit.sum()),
           "n_dirs_s_gt_3": w["spectrum"]["n_above_3"],
           "max_wald_sd_ln": float(np.nanmax(w["sd"])), "median_wald_sd_ln": float(np.nanmedian(w["sd"]))}
    g_true = truth["r"] - truth["beta"]
    g_est = est["r"] - est["beta"]
    out["err_ln_g(r-beta)"] = float(np.log(g_est / g_true)) if g_est > 0 else float("nan")
    # Wald sd of ln g via delta method on (ln r, ln beta)
    if "r" in free and "beta" in free:
        ir, ib = free.index("r"), free.index("beta")
        grad = np.zeros(len(free)); grad[ir] = est["r"] / (est["r"] - est["beta"]); grad[ib] = -est["beta"] / (est["r"] - est["beta"])
        out["wald_sd_ln_g"] = float(np.sqrt(max(grad @ w["cov"] @ grad, 0))) if g_est > 0 else float("nan")
    for nm in free:
        out[f"err_{nm}"] = float(np.log(est[nm] / truth[nm])) if nm != "fD" else float(np.log(est["fD"] / truth["fD"]))
        out[f"wald_sd_{nm}"] = float(w["sd"][free.index(nm)])
    return out


def main():
    seeds = CFG["seeds"]
    items = []
    for gt in ("G1_baseline", "G2_fast", "G5_lowdormant", "G6_stimulated"):
        th = truth(gt)
        for obs in ("N", "PDe"):
            for rep in range(6):
                for mname in MODELS:
                    items.append((gt, th, obs, seeds["synthetic_data"] + 500 + rep, mname))
    log(f"{len(items)} fits")
    out = parallel_map(worker, items, workers=4, chunksize=1)
    df = pd.DataFrame(out)
    save_table(df, "p4_model_reduction_replicates")
    summ = []
    for (gt, obs, m), g in df.groupby(["truth", "observed", "model"], sort=False):
        full = df[(df.truth == gt) & (df.observed == obs) & df.model.str.startswith("M9")].set_index("seed")
        gi = g.set_index("seed")
        dAIC = (gi.AIC - full.AIC.reindex(gi.index))
        summ.append({"truth": gt, "observed": obs, "model": m, "k": int(g.k.iloc[0]), "median_AIC_minus_AIC_full": float(dAIC.median()),
                     "frac_replicates_AIC_not_worse_than_full+2": float((dAIC <= 2).mean()), "median_chi2": float(g.chi2.median()),
                     "median_T_rel_error": float(g.T_rel_error.median()), "median_abs_lnN_err_240": float(g.abs_lnN_error_at_240.median()),
                     "median_err_ln_g": float(g["err_ln_g(r-beta)"].abs().median()), "median_wald_sd_ln_g": float(g["wald_sd_ln_g"].median()) if "wald_sd_ln_g" in g else float("nan"),
                     "median_max_wald_sd": float(g.max_wald_sd_ln.median()), "median_n_bound_hits": float(g.n_bound_hits.median()),
                     "median_n_dirs_s_gt_3": float(g.n_dirs_s_gt_3.median())})
    sm = pd.DataFrame(summ)
    save_table(sm, "p4_model_reduction_summary")
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(sm[["truth", "observed", "model", "k", "median_AIC_minus_AIC_full", "frac_replicates_AIC_not_worse_than_full+2", "median_T_rel_error", "median_err_ln_g", "median_wald_sd_ln_g", "median_max_wald_sd"]].round(3).to_string())

    # figure: reduction accept/reject (AIC difference and g precision)
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
    names = list(MODELS)
    short = [m.split(" ")[0] for m in names]
    cases = [(gt, obs) for gt in ("G1_baseline", "G2_fast", "G5_lowdormant", "G6_stimulated") for obs in ("N", "PDe")]
    for ax, (key, ttl, logy) in zip(axes, [("median_AIC_minus_AIC_full", "AIC relative to the full model (lower = reduced model preferred)", False),
                                           ("median_T_rel_error", "crossing-time relative error", True), ("median_wald_sd_ln_g", "Wald sd of ln(r - beta) (precision of the well-determined combination)", True)]):
        for i, (gt, obs) in enumerate(cases):
            vals = [sm[(sm.truth == gt) & (sm.observed == obs) & (sm.model == m)][key].iloc[0] for m in names]
            ax.plot(range(len(names)), vals, "o-", ms=4, label=f"{gt.split('_')[0]} obs {obs}", color=plt.cm.tab10(i), ls="-" if obs == "N" else "--")
        ax.set_xticks(range(len(names))); ax.set_xticklabels(short)
        if logy:
            ax.set_yscale("log")
        else:
            ax.axhline(0, color=INK, lw=0.8)
        ax.set_title(ttl, fontsize=8)
    axes[0].legend(fontsize=6, ncol=2)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p4_footer(fig, "Nested reduction (M9 full -> M7 immune frozen -> M6 delta_D=0 -> M5 fD at equilibrium ratio); data generated from the FULL model.")
    fig.savefig(FIG / "fig_p4_13_model_reduction.png"); plt.close(fig)
    log("done")


if __name__ == "__main__":
    main()
