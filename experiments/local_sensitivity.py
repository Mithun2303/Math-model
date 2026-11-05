"""Phase 3, Task 1: local sensitivity (central log-derivatives; event indicators handled by flip margins).

Run: .venv/bin/python experiments/local_sensitivity.py
Outputs: results/tables/p3_local_sensitivity*.csv, results/figures/fig_p3_01_local_sensitivity_*.png
All quantities are simulated; reference points and ranges are hypothetical.
"""
import math

from _p3common import *  # noqa: F401,F403
import matplotlib.pyplot as plt

from dormancy.equilibria import enumerate_equilibria
from dormancy.model import jacobian
from dormancy.normalized import NormalizedParams
from dormancy.sensitivity import flip_margin, log_derivative, perturb_row
from dormancy.uncertainty import (default_options, long_run, regime_of, row_to_model, run_trajectory)

INPUTS = ["r", "q", "a", "delta_D", "beta", "rho", "eta", "P0", "D0", "N_detect"]
LABEL = {"r": "r", "q": "q", "a": "a", "delta_D": "$\\delta_D$", "beta": "$\\beta$", "rho": "$\\rho$", "eta": "$\\eta$",
         "P0": "$P_0$", "D0": "$D_0$", "N_detect": "$N_{detect}$"}
H = 0.05                                    # log step; also H/2 for the agreement check
base = dict(CFG["baseline"])
opt = default_options(CFG)
# three reference points with different outcomes (beta varied only; every other input at baseline)
REFS = {"R_early (beta=0.020)": 0.020, "R_baseline (beta=0.040)": 0.040, "R_late (beta=0.0567)": 0.0567}
LOG_OUT = ["N_60", "N_120", "N_240", "T_event", "N_attractor", "dormancy_duration"]
RAW_OUT = ["lambda_TF", "lambda_attractor"]


def outputs(row: dict) -> dict:
    """All local-sensitivity outputs for one parameter dict (deterministic)."""
    p, y0, thr = row_to_model(row, CFG["fixed"])
    t = run_trajectory(p, y0, [thr], CFG["horizon_months"], main_threshold=thr,
                       dorm_eps=opt["dorm_eps"], dorm_min_frac=opt["dorm_min_frac"])
    k = f"{thr:g}".replace(".", "p")
    reg = regime_of(p)
    lr = long_run(p, y0, thr, with_event=False, regime=reg)
    out = {"N_60": t["N_60"], "N_120": t["N_120"], "N_240": t["N_240"], "T_event": t.get(f"T_{k}", float("nan")),
           "event": t.get(f"event_{k}", float("nan")), "dormancy_duration": t["dormancy_duration"],
           "dormancy_exited": t["dormancy_ever_exited"], "lambda_TF": reg["lambda_tf"],
           "N_attractor": lr["N_inf"], "attractor": lr["attractor"], "regime": reg["regime"]}
    # eigenvalue at the attractor reached (TF or the matching stable nonzero equilibrium)
    lam = float("nan")
    if lr["attractor"] == "TF":
        lam = reg["lambda_tf"]
    elif np.isfinite(lr["N_inf"]):
        for e in enumerate_equilibria(p):
            if abs(e.P + e.D - lr["N_inf"]) < 1e-6 * max(1, lr["N_inf"]):
                lam = float(np.max(np.linalg.eigvals(jacobian(e.state, p)).real))
    out["lambda_attractor"] = lam
    # fold proximity (N* jumps at folds): smallest |dk/dP|-type indicator replaced by a robust flag
    out["n_stable_nonzero"] = reg["n_stable_nonzero"]
    for hz in (60, 120, 240):
        out[f"event_{hz}"] = float(1 if (out["event"] == 1 and out["T_event"] <= hz) else 0) if np.isfinite(out["event"]) else float("nan")
    return out


def main():
    rows = []
    margins = []
    for rname, beta in REFS.items():
        ref = dict(base); ref["beta"] = beta
        r0 = outputs(ref)
        log(f"{rname}: event={r0['event']} T={r0['T_event']:.1f} N240={r0['N_240']:.4g} regime={r0['regime']} "
            f"attractor={r0['attractor']} N*={r0['N_attractor']:.4g}")
        for name in INPUTS:
            stencil = {s: outputs(perturb_row(ref, name, s)) for s in (H, -H, H / 2, -H / 2)}
            for o in LOG_OUT + RAW_OUT:
                vals = {s: stencil[s][o] for s in stencil}
                # event-time derivative only where the event exists at ALL stencil points (else 'undefined')
                if o == "T_event" and not all(stencil[s]["event"] == 1 for s in stencil):
                    d, st = float("nan"), "undefined"
                elif o == "dormancy_duration" and not all(stencil[s]["dormancy_exited"] for s in stencil):
                    d, st = float("nan"), "undefined"          # right-censored at the horizon: duration >= 240, not a number
                elif o in ("N_attractor", "lambda_attractor") and len({stencil[s]["attractor"] for s in stencil}) > 1:
                    d, st = float("nan"), "undefined"          # attractor switches inside the stencil: derivative meaningless
                elif o == "N_attractor" and r0["N_attractor"] == 0:
                    d, st = float("nan"), "undefined"
                else:
                    d, st = log_derivative(vals, H, log_output=(o in LOG_OUT), tol=0.2)
                if name == "N_detect" and o in ("N_60", "N_120", "N_240", "N_attractor", "lambda_TF", "lambda_attractor"):
                    st = "zero" if st == "ok" else st          # N_detect does not enter the trajectory: exact zero by construction
                    d = 0.0 if st == "zero" else d
                rows.append({"reference": rname, "beta_ref": beta, "input": name, "output": o, "derivative": d, "status": st,
                             "kind": "elasticity d ln(out)/d ln(theta)" if o in LOG_OUT else "d out / d ln(theta) [1/month]"})
            # flip margins for event indicators (never differentiated)
            for hz in (60, 120, 240):
                ind = lambda row, hz=hz: outputs(row)[f"event_{hz}"]
                fm = flip_margin(ind, ref, name, int(r0[f"event_{hz}"]), max_log=math.log(10.0), n_scan=60, n_bisect=40)
                margins.append({"reference": rname, "input": name, "horizon": hz, "base_event": int(r0[f"event_{hz}"]),
                                "flip_factor_up": fm["up_factor"], "flip_factor_down": fm["down_factor"]})
        log(f"   {rname}: stencils done")
    df = pd.DataFrame(rows)
    save_table(df, "p3_local_sensitivity")
    mg = pd.DataFrame(margins)
    save_table(mg, "p3_local_flip_margins")
    log("status counts: " + str(df.status.value_counts().to_dict()))

    # operational dormancy definition: sensitivity of the duration to its own thresholds (baseline)
    ref = dict(base)
    p, y0, thr = row_to_model(ref, CFG["fixed"])
    dd = []
    for eps in (0.005, 0.01, 0.02):
        for fr in (0.3, 0.5, 0.7):
            t = run_trajectory(p, y0, [thr], 240.0, main_threshold=thr, dorm_eps=eps, dorm_min_frac=fr)
            dd.append({"max_abs_dlnN_dt": eps, "min_dormant_fraction": fr, "dormancy_duration_months": t["dormancy_duration"],
                       "T_event": t.get(f"T_{'0p1'}", float("nan"))})
    save_table(pd.DataFrame(dd), "p3_dormancy_definition_sensitivity")

    # ---------------- figures
    # tornado of elasticities of N(240) and T_event, per reference
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.2))
    for c, rname in enumerate(REFS):
        for r_, (o, ttl) in enumerate([("N_240", "$d\\ln N(240)/d\\ln\\theta$"), ("T_event", "$d\\ln T_{first\\ crossing}/d\\ln\\theta$")]):
            ax = axes[r_, c]
            sub = df[(df.reference == rname) & (df.output == o)].set_index("input").reindex(INPUTS)
            order = sub.derivative.abs().fillna(-1).sort_values().index
            y = np.arange(len(order))
            for yi, nm in zip(y, order):
                d, st = sub.loc[nm, "derivative"], sub.loc[nm, "status"]
                if st in ("ok", "zero"):
                    ax.barh(yi, d, color=ORANGE if d > 0 else BLUE, height=0.65)
                else:
                    ax.text(0, yi, "undefined" if st == "undefined" else "unstable (near boundary)", va="center", ha="left",
                            fontsize=7, color=INK2, style="italic")
            ax.set_yticks(y); ax.set_yticklabels([LABEL[n] for n in order])
            ax.axvline(0, color=INK, lw=0.8)
            ax.set_title(f"{rname}\n{ttl}", fontsize=9)
            ax.set_xlabel("elasticity (dimensionless)")
    fig.suptitle("Local sensitivity (central log-derivatives, step 5% and 2.5% must agree within 20%)", y=0.995, fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97)); p3_footer(fig, "Italic: derivative not defined (event absent at a stencil point) or unstable (near event boundary).")
    fig.savefig(FIG / "fig_p3_01_local_sensitivity_tornado.png"); plt.close(fig)

    # heatmap of all outputs for the baseline reference + eigenvalue sensitivities
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    outs = ["N_60", "N_120", "N_240", "T_event", "N_attractor", "dormancy_duration", "lambda_TF", "lambda_attractor"]
    for ax, rname in zip(axes, REFS):
        M = np.full((len(INPUTS), len(outs)), np.nan)
        S = [[""] * len(outs) for _ in INPUTS]
        for i, nm in enumerate(INPUTS):
            for j, o in enumerate(outs):
                rr = df[(df.reference == rname) & (df.input == nm) & (df.output == o)].iloc[0]
                if rr.status in ("ok", "zero"):
                    M[i, j] = rr.derivative; S[i][j] = f"{rr.derivative:.2g}"
                else:
                    S[i][j] = "n/d" if rr.status == "undefined" else "~"
        # colour scale per column (outputs have different units)
        Mn = np.where(np.isnan(M), np.nan, M / np.nanmax(np.abs(M), axis=0, keepdims=True).clip(min=1e-12))
        ax.imshow(Mn, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        for i in range(len(INPUTS)):
            for j in range(len(outs)):
                ax.text(j, i, S[i][j], ha="center", va="center", fontsize=6.5, color=INK)
        ax.set_xticks(range(len(outs))); ax.set_xticklabels(["N(60)", "N(120)", "N(240)", "T", "N*", "dorm.\ndur.", "$\\lambda_{TF}$", "$\\lambda_{att}$"], fontsize=7)
        ax.set_yticks(range(len(INPUTS))); ax.set_yticklabels([LABEL[n] for n in INPUTS])
        ax.set_title(rname, fontsize=9); ax.grid(False)
    fig.suptitle("Local sensitivity map: elasticities (first 6 columns), $d\\lambda/d\\ln\\theta$ [1/month] (last 2); n/d = undefined, ~ = unstable", y=0.999, fontsize=9.5)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96)); p3_footer(fig, "Colour scaled within each column.")
    fig.savefig(FIG / "fig_p3_01b_local_sensitivity_map.png"); plt.close(fig)

    # flip-margin figure: smallest relative change (up or down) that flips 'event by 240 months'
    fig, ax = plt.subplots(figsize=(9.5, 4.4))
    sub = mg[(mg.reference == "R_baseline (beta=0.040)") & (mg.horizon == 240)]
    pct, direc = [], []
    for _, r_ in sub.iterrows():
        cand = []
        if np.isfinite(r_.flip_factor_up):
            cand.append(((r_.flip_factor_up - 1) * 100, "increase"))
        if np.isfinite(r_.flip_factor_down):
            cand.append(((1 - r_.flip_factor_down) * 100, "decrease"))
        if cand:
            v, d = min(cand)
            pct.append(v); direc.append(d)
        else:
            pct.append(np.nan); direc.append("none within x10")
    ax.bar(range(len(sub)), [v if np.isfinite(v) else 0 for v in pct],
           color=[(ORANGE if d == "increase" else BLUE) if np.isfinite(v) else GREY for v, d in zip(pct, direc)])
    for i, (v, d) in enumerate(zip(pct, direc)):
        ax.text(i, (v if np.isfinite(v) else 0) + 2, f"{'+' if d == 'increase' else '-'}{v:.0f}%" if np.isfinite(v) else "no flip\n(x10)", ha="center", fontsize=7.5)
    ax.set_xticks(range(len(sub))); ax.set_xticklabels([LABEL[n] for n in sub.input])
    ax.set_ylabel("smallest relative change of the input\nthat flips 'event by 240 months' [%]")
    ax.set_title("Distance of the baseline reference to the 240-month event boundary (orange: increase, blue: decrease)", fontsize=9)
    fig.tight_layout(rect=(0, 0.05, 1, 1)); p3_footer(fig, "Baseline reference beta=0.04: first crossing ~220 months (close to the horizon).")
    fig.savefig(FIG / "fig_p3_01c_local_flip_margins.png"); plt.close(fig)
    log("done")


if __name__ == "__main__":
    main()
