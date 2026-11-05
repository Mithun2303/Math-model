"""Experiment 1/5: baseline trajectories for weak, near-boundary and strong immune killing."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from dormancy import load_config, simulate, critical_k, nonzero_equilibria, eigenvalues, classify

ROOT = Path(__file__).resolve().parents[1]
cfg = load_config(ROOT / "configs" / "baseline.yaml")
p0, ic, rel, sv = cfg["params"], cfg["initial_conditions"], cfg["relapse"], cfg["solver"]
y0 = [ic["P0"], ic["D0"], ic["E0"] if ic["E0"] is not None else p0.E_free]
kc = critical_k(p0)
scenarios = {"weak (0.25 kc)": 0.25 * kc, "moderate (0.5 kc)": 0.5 * kc,
             "near boundary (0.98 kc)": 0.98 * kc, "strong (2 kc)": 2 * kc}

rows = []
fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
for (label, k), ax in zip(scenarios.items(), axes.flat):
    p = p0.with_(k=k)
    sim = simulate(p, y0, rel["N_detect"], rel["horizon"], sv["method"], sv["rtol"], sv["atol"])
    ax.semilogy(sim.t, sim.P, label="P (proliferating)")
    ax.semilogy(sim.t, sim.D, label="D (dormant)")
    ax.semilogy(sim.t, sim.N, "k", label="N = P + D")
    ax.axhline(rel["N_detect"], color="r", ls="--", lw=0.8, label="N_detect")
    ax.set_title(f"{label}: k = {k:.3f}" + (f", relapse at {sim.time:.0f} mo" if sim.event else ", no event"))
    ax.set_ylim(1e-8, 2)
    eqs = nonzero_equilibria(p)
    rows.append({"scenario": label, "k": k, "event": sim.event, "time_months": sim.time,
                 "nonzero_equilibria": [(round(e.P, 5), round(e.D, 5), round(e.E, 4),
                                         classify(eigenvalues(e.state, p))) for e in eqs]})
for ax in axes[1]:
    ax.set_xlabel("time (months)")
axes[0, 0].legend(fontsize=7)
fig.suptitle("Simulated (hypothetical) dynamics, illustrative parameters")
fig.tight_layout()
(ROOT / "results" / "figures").mkdir(parents=True, exist_ok=True)
fig.savefig(ROOT / "results" / "figures" / "baseline_trajectories.png", dpi=150)
df = pd.DataFrame(rows)
df.to_csv(ROOT / "results" / "tables" / "baseline_scenarios.csv", index=False)
print(f"critical k = {kc:.4f}")
print(df.to_string())
