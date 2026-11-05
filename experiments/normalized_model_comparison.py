"""Phase 2 experiment (M1): original vs normalised model, and invariance under the E-gauge (figure H).

Run:  .venv/bin/python experiments/normalized_model_comparison.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _style import *  # noqa: E402,F403
import matplotlib.pyplot as plt  # noqa: E402

from dormancy import NormalizedParams, Params, critical_k, load_config, simulate  # noqa: E402
from dormancy.normalized import (gauge_transform, simulate_normalized, to_normalized_state,  # noqa: E402
                                 to_original_state)

cfg = load_config(ROOT / "configs" / "baseline.yaml")
p0 = cfg["params"]
import yaml  # noqa: E402
with open(ROOT / "configs" / "normalized.yaml") as fh:
    ny = yaml.safe_load(fh)
n0 = NormalizedParams(**ny["parameters"])
ic = ny["initial_conditions"]
kc = critical_k(p0)
H = 3000.0

# 1) equivalence of the shipped configs: baseline.yaml (k=0.1, s=0.02) <-> normalized.yaml (beta=0.04)
assert abs(NormalizedParams.from_params(p0).beta - n0.beta) < 1e-15, "configs disagree on beta"

rows = []
trajs = {}
for f in [0.25, 0.5, 0.75, 0.98]:
    for s in [0.002, 0.02, 0.2, 2.0]:                 # gauge family: same beta, different (k, s)
        p = gauge_transform(p0.with_(k=f * kc), s / p0.s)
        n = NormalizedParams.from_params(p)
        y0 = [1e-6, 1e-4, p.E_free]                    # E0 = tumour-free level in each gauge
        so = simulate(p, y0, 0.1, H, n_points=3000)
        t, z, ev, T = simulate_normalized(n, to_normalized_state(y0, p), 0.1, H, n_points=3000)
        m = min(len(so.t), len(t))
        # compare on each solution's own accepted grid via the event time and the end state
        dstate = np.max(np.abs(so.y[:2, -1] - z[:2, -1])) if (so.event and ev) else np.nan
        rows.append({"k_over_kc": f, "s": s, "k_original": p.k, "beta": n.beta, "T_original": so.time if so.event else np.nan,
                     "T_normalized": T if ev else np.nan, "abs_diff_T": abs(so.time - T) if (so.event and ev) else np.nan,
                     "event_state_abs_diff_PD": dstate,
                     "E_end_over_E_free_orig": so.E[-1] / p.E_free, "e_end_norm": z[2, -1]})
        trajs[(f, s)] = (so, p, (t, z), n)
df = pd.DataFrame(rows)
df.to_csv(TAB / "normalized_model_equivalence.csv", index=False)
print(df.to_string())

# figure H
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
f = 0.5
cols = {0.002: BLUE, 0.02: ORANGE, 0.2: AQUA, 2.0: VIOLET}
ax = axes[0]
for s, c in cols.items():
    so, p, (t, z), n = trajs[(f, s)]
    ax.semilogy(so.t, so.N, color=c, lw=2.2, alpha=0.8, label=f"original, s={s:g}, k={p.k:.3g}")
tn, zn = trajs[(f, 0.02)][2]
ax.semilogy(tn, zn[0] + zn[1], "k--", lw=1.2, label=r"normalised ($\beta$=" + f"{trajs[(f,0.02)][3].beta:.4f})")
ax.set_xlabel("time [months]"); ax.set_ylabel("burden $N=P+D$ [fraction of K]")
ax.set_title("H1. Burden: all gauges = normalised")
ax.legend(fontsize=6.5, loc="lower right")
ax = axes[1]
for s, c in cols.items():
    so, p, (t, z), n = trajs[(f, s)]
    ax.plot(so.t, so.E, color=c, lw=2.2, alpha=0.8, label=f"E, s={s:g}")
    ax.plot(t, z[2] * p.s / p.mu, "k--", lw=0.9)
ax.set_yscale("log")
ax.set_xlabel("time [months]"); ax.set_ylabel("immune level $E$ [E-units; scale is arbitrary]")
ax.set_title(r"H2. $E$ scales with $s$; $e=\mu E/s$ (dashed) matches")
ax.legend(fontsize=7)
ax = axes[2]
for s, c in cols.items():
    g = df[df.s == s]
    ax.semilogy(g.k_over_kc, g.abs_diff_T.clip(lower=1e-16), "o-", color=c, ms=5, label=f"s={s:g}")
ax.set_xlabel(r"$\beta/\beta_c$ = $k/k_c$ in every gauge (dimensionless)"); ax.set_ylabel("|relapse time original − normalised| [months]")
ax.set_title("H3. Relapse-time difference (solver tolerance level)")
ax.legend(fontsize=7)
fig.tight_layout(rect=(0, 0.05, 1, 1))
footer(fig, "Baseline parameters; (k,s)->(k/c,cs) leaves beta=ks/mu fixed. Hypothetical illustrative parameters; mathematical equivalence check, not clinical.")
fig.savefig(FIG / "fig_H_original_vs_normalized.png"); plt.close(fig)
print("max |dT| =", np.nanmax(df.abs_diff_T), " max |dPD| =", np.nanmax(df.event_state_abs_diff_PD))
