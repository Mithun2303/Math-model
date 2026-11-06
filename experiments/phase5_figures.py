"""Phase 5 Figure 10: clinical comparison summary (built only from saved Phase 5 tables; no new computation)."""
import sys, textwrap
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from _style import *  # noqa
import pandas as pd, numpy as np, matplotlib.pyplot as plt

h = pd.read_csv(TAB / "p5_horizon_comparison.csv", comment="#")
rf = pd.read_csv(TAB / "p5_consistency_robust_findings.csv", comment="#")
p1 = h[(h.prior == "P1_base") & (h.horizon_years.isin([5, 10, 15, 20]))].sort_values("horizon_years")
fig = plt.figure(figsize=(11, 5.6))
ax = fig.add_axes([0.07, 0.17, 0.40, 0.70])
x = np.arange(len(p1)); w = 0.32
ax.bar(x - w/2, p1.F_km, w, color="#000000", label="METABRIC RFS event (KM)", yerr=[p1.F_km - p1.km_ci_lo, p1.km_ci_hi - p1.F_km], capsize=3)
ax.bar(x + w/2, p1.F_sim_reference, w, color="#0072B2", label="Model, reference scenario (P1, θ=0.1, τ=d=0)")
for xi, lo, hi in zip(x + w/2, p1.F_sim_grid_min_prior, p1.F_sim_grid_max_prior):
    ax.plot([xi, xi], [lo, hi], color="#D55E00", lw=2.5)
ax.plot([], [], color="#D55E00", lw=2.5, label="Model range over P1 threshold/τ/d/death-first grid")
ax.set_xticks(x); ax.set_xticklabels([f"{int(y)} y" + ("*" if y == 15 else "") for y in p1.horizon_years])
ax.set_ylabel("Cumulative incidence (1 − event-free)"); ax.set_ylim(0, 0.95)
ax.set_title("(a) Observed vs simulated horizon estimates"); ax.legend(fontsize=7.5, loc="upper left")
ax2 = fig.add_axes([0.52, 0.17, 0.46, 0.70]); ax2.axis("off")
lines = ["(b) Pre-specified statements (1,764 scenarios)"]
for _, r in rf.iterrows():
    tag = "ROBUST" if str(r.verdict).startswith("ROBUST") else f"{100*r.frac_true:.1f}% of scenarios"
    lines += textwrap.wrap(f"{r.statement} [{tag}]: {r.text}", 78, subsequent_indent="    ")
lines += ["", "Reading: the model is below METABRIC at 5 and 10 y in all or",
          "nearly all scenarios and overshoots its late increment; this is",
          "expected (no treatment/early-recurrence mechanism, threshold",
          "event vs composite RFS, unknown origin). Consistency at one",
          "horizon is not evidence of mechanism. No parameter was fitted."]
ax2.text(0, 1, "\n".join(lines), va="top", fontsize=7.4, family="monospace")
fig.text(0.01, 0.01, "* 15 y has no Phase 4 reliability criterion. MODEL (hypothetical) vs CLINICAL (descriptive); not validation. Design hash 9081048b8863.", fontsize=7, color="#1c4f9b")
fig.text(0.995, 0.99, "MODEL (hypothetical) vs CLINICAL (descriptive)", fontsize=8, fontweight="bold", color="#1c4f9b", ha="right", va="top")
fig.savefig(FIG / "fig_p5_10_clinical_comparison_summary.png"); print("ok")
