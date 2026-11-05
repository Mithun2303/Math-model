"""Phase 2 experiment: relapse-time approximations vs numerical first-crossing time (figure G).

Run:  .venv/bin/python experiments/relapse_time_approximation.py
Outputs: results/tables/relapse_time_*.csv, results/figures/fig_G_relapse_time_approximation.png
"""
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _style import *  # noqa: E402,F403
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from dormancy import Params, critical_k, load_config, simulate  # noqa: E402
from dormancy.relapse import (T_manifold, T_matched, T_matched_parts, T_original, T_projected, T_two_mode,  # noqa: E402
                              k_detection_limit, unstable_mode)

warnings.filterwarnings("ignore", category=RuntimeWarning)
cfg = load_config(ROOT / "configs" / "baseline.yaml")
p0 = cfg["params"]
kc = critical_k(p0)
N0 = 1.01e-4                       # baseline total initial burden P0 + D0
HOR = 1e5                          # months; long enough to resolve near-critical events
t0 = time.time()

METHODS = {"A0_original": T_original, "A1_projected": T_projected, "A1b_two_mode_linear": T_two_mode,
           "A2_manifold_2nd_order": T_manifold, "A3_matched": T_matched}


def y0_for(p, fD, N0=N0):
    return [N0 * (1 - fD), N0 * fD, p.E_free]


def numeric(p, y0, Nd):
    s = simulate(p, y0, Nd, HOR, n_points=400)
    return s.time if s.event else np.nan, s


rows = []
fracs = np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.72, 0.74, 0.75, 0.8, 0.85, 0.9, 0.95])
for Nd in (1e-3, 1e-2, 0.1):
    for fD in (0.0, 0.5, 0.99):
        for f in fracs:
            p = p0.with_(k=f * kc)
            y0 = y0_for(p, fD)
            T_num, sim = numeric(p, y0, Nd)
            m = unstable_mode(p)
            row = {"N_detect": Nd, "dormant_fraction_0": fD, "k": p.k, "k_over_kc": f, "T_numeric": T_num,
                   "lambda_plus": m.lam, "gamma": m.gamma, "P_det_over_h": Nd / p.phi / p.h,
                   "k_det_over_kc": k_detection_limit(p, Nd) / kc,
                   "E_at_event_over_Efree": (sim.E[-1] / p.E_free) if sim.event else np.nan}
            for name, fn in METHODS.items():
                T = fn(p, y0, Nd)
                row[name] = T if np.isfinite(T) else np.nan
                row[name + "_relerr"] = (T - T_num) / T_num if (np.isfinite(T) and np.isfinite(T_num) and T_num > 0) else np.nan
            tot, t_closed, tail = T_matched_parts(p, y0, Nd)
            row["A3_closed_form_part"] = t_closed
            row["A3_numerical_tail"] = tail
            row["A3_tail_fraction_of_T"] = tail / tot if np.isfinite(tot) and tot > 0 else np.nan
            rows.append(row)
    print(f"[{time.time() - t0:5.1f}s] N_detect={Nd} done", flush=True)
df = pd.DataFrame(rows)
df.to_csv(TAB / "relapse_time_approximation.csv", index=False)

# --- attribution of the baseline discrepancy (k = 0.5 k_c, N_det = 0.1, baseline IC)
attr = []
for label, kw in [("baseline", {}), ("rho=0, eta=0 (no immune feedback)", dict(rho=0.0, eta=0.0)),
                  ("rho=0 only", dict(rho=0.0)), ("eta=0 only", dict(eta=0.0))]:
    p = p0.with_(k=0.5 * kc, **kw)
    y0 = [1e-6, 1e-4, p.E_free]
    T_num, sim = numeric(p, y0, 0.1)
    attr.append({"variant": label, "T_numeric": T_num, "A0": T_original(p, y0, 0.1), "A1": T_projected(p, y0, 0.1),
                 "A2": T_manifold(p, y0, 0.1), "A3": T_matched(p, y0, 0.1),
                 "E_at_event_over_Efree": sim.E[-1] / p.E_free,
                 "lambda_plus": unstable_mode(p).lam})
attr = pd.DataFrame(attr)
attr["gap_A0_to_numeric"] = attr.T_numeric - attr.A0
attr["recovered_by_projection_A1_minus_A0"] = attr.A1 - attr.A0
attr["remaining_after_projection"] = attr.T_numeric - attr.A1
attr.to_csv(TAB / "relapse_time_attribution.csv", index=False)
print(attr.round(3).to_string())

# --- accuracy by method and threshold, restricted to the 'safe' region k <= 0.9 k_det (finite numeric time)
summ = []
for Nd, g in df.groupby("N_detect"):
    sel = g[(g.k_over_kc <= 0.9 * g.k_det_over_kc) & g.T_numeric.notna()]
    for name in METHODS:
        e = sel[name + "_relerr"]
        summ.append({"N_detect": Nd, "method": name, "n": int(e.notna().sum()), "n_predicts_no_event": int(e.isna().sum()),
                     "median_abs_relerr": float(np.nanmedian(np.abs(e))) if e.notna().any() else np.nan,
                     "max_abs_relerr": float(np.nanmax(np.abs(e))) if e.notna().any() else np.nan})
summ = pd.DataFrame(summ)
summ.to_csv(TAB / "relapse_time_accuracy_summary.csv", index=False)
print(summ.round(4).to_string())

# --- fine scan around k_det: is "no event beyond k_det" true? (monotone vs overshooting approach to N*)
rows = []
for Nd in (1e-3, 1e-2, 0.1):
    kd = k_detection_limit(p0, Nd)
    for fD in (0.0, 0.5, 0.99):
        for f in [0.98, 0.99, 0.995, 0.999, 1.0, 1.001, 1.005, 1.01, 1.02, 1.03, 1.04, 1.06]:
            p = p0.with_(k=f * kd)
            T_num, sim = numeric(p, y0_for(p, fD), Nd)
            rows.append({"N_detect": Nd, "k_det_over_kc": kd / kc, "k_over_kdet": f, "dormant_fraction_0": fD,
                         "T_numeric": T_num, "event": int(np.isfinite(T_num)), "max_N_along_trajectory": float(sim.N.max()),
                         "horizon_months": HOR})
kd_df = pd.DataFrame(rows)
kd_df.to_csv(TAB / "relapse_time_kdet_boundary.csv", index=False)
print("events at k > 1.001 k_det:", int(kd_df[kd_df.k_over_kdet > 1.001].event.sum()), "of", int((kd_df.k_over_kdet > 1.001).sum()))

# --- figure G
colors = {"A0_original": GREY, "A1_projected": BLUE, "A2_manifold_2nd_order": ORANGE, "A3_matched": AQUA,
          "A1b_two_mode_linear": VIOLET}
styles = {"A0_original": "--", "A1_projected": "-.", "A2_manifold_2nd_order": ":", "A3_matched": "-", "A1b_two_mode_linear": (0, (1, 1))}
labels = {"A0_original": "A0 original $\\ln(N_{det}/N_0)/\\lambda_+$", "A1_projected": "A1 eigen-projected",
          "A2_manifold_2nd_order": "A2 2nd-order manifold", "A3_matched": "A3 SEMI-analytic: closed form to N=1e-3 + numerical tail",
          "A1b_two_mode_linear": "A1b two-mode linear"}
fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
ax = axes[0]
g = df[(df.N_detect == 0.1) & (df.dormant_fraction_0 == 0.99)].sort_values("k_over_kc")
ax.plot(g.k_over_kc, g.T_numeric, "o-", color=INK, ms=5, label="numerical first-crossing time")
for name in ["A0_original", "A1_projected", "A2_manifold_2nd_order", "A3_matched"]:
    ax.plot(g.k_over_kc, g[name], color=colors[name], ls=styles[name], label=labels[name])
kd = g.k_det_over_kc.iloc[0]
ax.axvline(kd, color=MAGENTA, lw=1.2); ax.axvspan(kd, 1.0, color=MAGENTA, alpha=0.07)
ax.text(kd + 0.01, 110, "$k>k_{det}$: $N^*<N_{det}$,\nno event in simulation;\nA0/A1 still return finite $T$", fontsize=7.5, color=INK2)
ax.set_yscale("log"); ax.set_xlim(0, 1.0)
ax.set_xlabel("$k/k_c$ (dimensionless)"); ax.set_ylabel("relapse time [months]")
ax.set_title("G1. $N_{det}=0.1$, initial dormant fraction 0.99")
ax.legend(loc="upper left", fontsize=7)
ax = axes[1]
for Nd, ls in [(1e-3, "-"), (1e-2, "--"), (0.1, ":")]:
    gg = df[(df.N_detect == Nd) & (df.dormant_fraction_0 == 0.99) & (df.k_over_kc <= 0.9 * df.k_det_over_kc)]
    for name in ["A0_original", "A1_projected", "A2_manifold_2nd_order", "A3_matched"]:
        ax.plot(gg.k_over_kc, 100 * gg[name + "_relerr"], color=colors[name], ls=ls, marker="o", ms=2.5, lw=1.3)
ax.axhline(0, color=INK, lw=0.7)
ax.set_ylim(-60, 60); ax.set_xlim(0, 1.0)
ax.set_xlabel("$k/k_c$ (dimensionless)"); ax.set_ylabel("relative error (estimate − numeric)/numeric [%]")
ax.set_title("G2. Error by method; line style = threshold")
leg = [Line2D([], [], color=colors[n], lw=2, label=labels[n]) for n in ["A0_original", "A1_projected", "A2_manifold_2nd_order", "A3_matched"]]
leg += [Line2D([], [], color=INK, ls=ls, label=f"$N_{{det}}$={Nd:g}") for Nd, ls in [(1e-3, "-"), (1e-2, "--"), (0.1, ":")]]
ax.legend(handles=leg, fontsize=6.5, loc="upper right", ncol=2)
ax = axes[2]
gg = df[(df.k_over_kc == 0.5) & (df.N_detect == 0.1)].sort_values("dormant_fraction_0")
ax.plot(gg.dormant_fraction_0, gg.T_numeric, "o-", color=INK, label="numerical")
for name in ["A0_original", "A1_projected", "A3_matched"]:
    ax.plot(gg.dormant_fraction_0, gg[name], color=colors[name], ls=styles[name], marker="s", ms=4, label=labels[name])
ax.set_xlabel("initial dormant fraction $D_0/(P_0+D_0)$ (dimensionless)"); ax.set_ylabel("relapse time [months]")
ax.set_title("G3. Initial composition at $k=0.5k_c$, fixed $N_0=1.01\\times10^{-4}$")
ax.legend(fontsize=7)
fig.tight_layout(rect=(0, 0.05, 1, 1))
footer(fig, "Baseline parameters (see configs/baseline.yaml), E0=s/mu, numerical = DOP853 first upward crossing (rtol 1e-9). Hypothetical mathematical comparison, not clinical.")
fig.savefig(FIG / "fig_G_relapse_time_approximation.png"); plt.close(fig)
print(f"[{time.time() - t0:5.1f}s] done")
