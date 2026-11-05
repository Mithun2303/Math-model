"""Phase 3, Task 3: basins of attraction for the two kinds of bistability found in Phase 2, and their link to relapse events.

Run: .venv/bin/python experiments/basin_analysis.py
Cases (hypothetical parameter sets):
  A   tumour-free + nonzero  : baseline with eta = 0.6, k = 0.572 (beta = k s/mu = 0.2288)       [Phase 2, Fig B2]
  B1  two nonzero (baseline-adjacent): eta = 0.164, k = 0.27 (beta = 0.108)                        [Phase 2, sect. 4.2]
  B2  two nonzero (well separated)   : a=0.306, q=1.58e-3, rho=0.240, eta=1.637, k=0.15 (beta=0.06) [Phase 2 scan example]
Outputs: results/tables/p3_basin_*.csv, results/figures/fig_p3_04/05*.png
"""
from _p3common import *  # noqa: F401,F403
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D

from dormancy.basins import (basin_map, bisect_boundary, boundary_mask, classify_ic, jitter_flip_rate,
                            random_initial_conditions, slice_initial_conditions, tight_options)
from dormancy.equilibria import enumerate_equilibria
from dormancy.model import jacobian
from dormancy.normalized import NormalizedParams
from dormancy.uncertainty import default_options, parallel_map, regime_of

CASES = {
    "A": dict(title="A: tumour-free + nonzero", kind="TF+NZ", eta=0.6, beta=0.4 * 0.572083602, thresholds=[0.3, 0.5]),
    "B1": dict(title="B1: two nonzero (eta=0.164)", kind="NZ+NZ", eta=0.164, beta=0.4 * 0.27, thresholds=[0.3, 0.5]),
    "B2": dict(title="B2: two nonzero (strong stimulation)", kind="NZ+NZ", a=0.3062386609, q=0.0015841984, rho=0.2404914173,
               eta=1.6373830589, beta=0.4 * 0.15, thresholds=[0.01, 0.3]),
}
COL = {"TF": BLUE, "NZ": ORANGE, "NZ-low": AQUA, "NZ-high": VIOLET, "unresolved": "#444444"}
NG = 60
ND = 0.1


def model_for(case: dict):
    base = dict(r=0.12, q=0.04, a=0.005, delta_D=0.001, beta=0.04, rho=0.04, eta=0.08)
    base.update({k: v for k, v in case.items() if k in base})
    return NormalizedParams(**base).to_params()


def stable_and_saddles(p):
    st, sd = [], []
    for e in enumerate_equilibria(p):
        ev = np.linalg.eigvals(jacobian(e.state, p))
        (st if np.all(ev.real < 0) else sd).append(e)
    return st, sd


def draw_map(ax, df, kind, order, p, case, a_vals, b_vals):
    codes = df.attractor.map({l: i for i, l in enumerate(order)}).fillna(len(order)).values.reshape(NG, NG)
    cols = [COL[o] for o in order] + [COL["unresolved"]]
    ax.pcolormesh(a_vals, b_vals, codes.T, cmap=ListedColormap(cols), vmin=-0.5, vmax=len(cols) - 0.5, shading="nearest", alpha=0.55, rasterized=True)
    bm = boundary_mask(codes.astype(int))
    ia, ib = np.where(bm)
    ax.plot(a_vals[ia], b_vals[ib], ".", color=INK, ms=1.3, alpha=0.7, rasterized=True)
    # events by 240 months for the extra thresholds (initial states strictly below the threshold)
    for th, mk, cc in zip(case["thresholds"], ["o", "^"], [INK, MAGENTA]):
        key = f"event240_{th:g}".replace(".", "p")
        if key in df:
            m = (df[key] == 1).values.reshape(NG, NG)
            ia, ib = np.where(m)
            ax.plot(a_vals[ia], b_vals[ib], mk, mfc="none", mec=cc, ms=3.2, mew=0.7, rasterized=True)
    st, sd = stable_and_saddles(p)
    for e in st:
        x = (e.P, e.D) if kind == "PD" else (e.P + e.D, e.E / p.E_free)
        ax.plot(*x, "*", color="white", mec=INK, ms=11, mew=1.0, zorder=6)
    for e in sd:
        x = (e.P, e.D) if kind == "PD" else (e.P + e.D, e.E / p.E_free)
        ax.plot(*x, "X", color="white", mec=INK, ms=8, mew=1.0, zorder=6)
    if kind == "PD":
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("initial $P_0$ [fraction of K]"); ax.set_ylabel("initial $D_0$ [fraction of K]")
        pp = np.geomspace(a_vals.min(), a_vals.max(), 50)
        ax.plot(pp, np.maximum(ND - pp, 1e-12), color=INK2, ls="--", lw=0.9)         # N0 = N_detect line
    else:
        ax.set_xscale("log")
        ax.set_xlabel("initial burden $N_0=P_0+D_0$ [fraction of K]"); ax.set_ylabel("initial immune level $e_0=E_0/E_{free}$")
        ax.axvline(ND, color=INK2, ls="--", lw=0.9)
    ax.set_xlim(a_vals.min(), a_vals.max()); ax.set_ylim(b_vals.min(), b_vals.max()); ax.grid(False)


def make_figures(results, cross):
    # ------------------------------------------------------------------ figures
    def fig_for(cnames, fname, title):
        fig, axes = plt.subplots(len(cnames), 2, figsize=(11.5, 4.6 * len(cnames)))
        axes = np.atleast_2d(axes)
        for r_, cname in enumerate(cnames):
            p, maps, order, _ = results[cname]
            for c, kind in enumerate(("PD", "NE")):
                a_vals, b_vals, ys, df = maps[kind]
                draw_map(axes[r_, c], df, kind, order, p, CASES[cname], a_vals, b_vals)
                axes[r_, c].set_title(f"{CASES[cname]['title']} - slice {'$(P_0,D_0)$, $e_0=1$' if kind=='PD' else '$(N_0,e_0)$, $f_D=q/(q+a+\\delta_D)$'}", fontsize=8.5)
        for r_, cname in enumerate(cnames):
            order_ = results[cname][2]
            leg = [Line2D([], [], marker="s", ls="", color=COL[o], alpha=0.6, label=f"basin of {o}") for o in order_]
            leg += [Line2D([], [], marker="*", ls="", color="white", mec=INK, ms=10, label="stable equilibrium"),
                    Line2D([], [], marker="X", ls="", color="white", mec=INK, ms=7, label="saddle"),
                    Line2D([], [], marker=".", ls="", color=INK, label="boundary cells"),
                    Line2D([], [], ls="--", color=INK2, label="$N_0=0.1$ (=N_detect)")]
            for th, mk, cc in zip(CASES[cname]["thresholds"], ["o", "^"], [INK, MAGENTA]):
                leg.append(Line2D([], [], marker=mk, ls="", mfc="none", mec=cc, label=f"event by 240 mo, $N_{{detect}}$={th:g}"))
            axes[r_, 1].legend(handles=leg, fontsize=6, loc="lower right", frameon=True, facecolor="white", framealpha=0.92)
        fig.suptitle(title, y=0.999, fontsize=10)
        fig.tight_layout(rect=(0, 0.035, 1, 0.97)); p3_footer(fig, "Attractor = long-horizon end state (LSODA, >=1e5 months) matched to stable equilibria; 60x60 grid per slice.")
        fig.savefig(FIG / fname); plt.close(fig)

    fig_for(["A"], "fig_p3_04_basins_tumour_free_nonzero.png", "Basins of attraction: tumour-free + nonzero bistability (hypothetical parameters, $\\eta=0.6$)")
    fig_for(["B1", "B2"], "fig_p3_05_basins_two_nonzero.png", "Basins of attraction: two stable nonzero states (hypothetical parameters)")

    # event vs attractor summary figure
    ct = cross
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, cname in zip(axes, CASES):
        sub = ct[ct.case == cname]
        ths = sorted(sub.threshold.unique())
        order = results[cname][2]
        w = 0.8 / len(order)
        for k, l in enumerate(order):
            vals = [sub[(sub.threshold == t) & (sub.attractor == l)] for t in ths]
            ntot = [int(v.n_starts_below_threshold.sum()) for v in vals]
            frac_ev = [float(v.event_by_240.sum() / n_) if n_ >= 10 else 0.0 for v, n_ in zip(vals, ntot)]     # n < 10: bar omitted
            ax.bar(np.arange(len(ths)) + k * w, frac_ev, width=w * 0.9, color=COL[l], label=f"ends in {l}")
            for i, (f_, n_) in enumerate(zip(frac_ev, ntot)):
                ax.text(i + k * w, f_ + 0.02, f"n={n_}" + (" (omitted)" if n_ < 10 else ""), ha="center", fontsize=6.5)
        ax.set_xticks(np.arange(len(ths)) + w * (len(order) - 1) / 2); ax.set_xticklabels([f"{t:g}" for t in ths])
        ax.set_xlabel("detection threshold $N_{detect}$ (model-defined)"); ax.set_ylabel("fraction with event by 240 months")
        ax.set_ylim(0, 1.35); ax.set_title(CASES[cname]["title"], fontsize=9); ax.legend(fontsize=7, loc="upper center", ncol=2)
    fig.suptitle("Event outcome vs attractor reached (random 3-D initial states strictly below the threshold; groups with n<10 omitted)", y=0.999, fontsize=10)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96)); p3_footer(fig)
    fig.savefig(FIG / "fig_p3_05b_basin_attractor_vs_event.png"); plt.close(fig)


def load_results():
    """Rebuild the plotting inputs from the saved grid tables (no simulation): used by --replot."""
    results = {}
    for cname, case in CASES.items():
        p = model_for(case)
        order = ["TF", "NZ"] if case["kind"] == "TF+NZ" else ["NZ-low", "NZ-high"]
        maps = {}
        for kind in ("PD", "NE"):
            a_vals, b_vals, ys = slice_initial_conditions(kind, p, NG, NG, lo=1e-6, hi=0.9)
            df = pd.read_csv(TAB / f"p3_basin_grid_{cname}_{kind}.csv", comment="#")
            maps[kind] = (a_vals, b_vals, ys, df)
        results[cname] = (p, maps, order, None)
    return results, pd.read_csv(TAB / "p3_basin_attractor_vs_event.csv", comment="#")


def main():
    import sys
    if "--replot" in sys.argv:
        results, cross = load_results()
        make_figures(results, cross)
        log("replotted from saved tables")
        return
    opts = default_options(CFG, long_floor=1e5)
    seeds = CFG["seeds"]
    summary, cross, lines, checks, ens_rows = [], [], [], [], []
    results = {}
    for cname, case in CASES.items():
        p = model_for(case)
        reg = regime_of(p)
        st, sd = stable_and_saddles(p)
        log(f"{cname}: regime={reg['regime']} stable N={['%.4g' % (e.P+e.D) for e in st]} saddle N={['%.4g' % (e.P+e.D) for e in sd]} TFstable={reg['tf_stable']}")
        o = dict(opts, extra_thresholds=case["thresholds"])
        order = ["TF", "NZ"] if case["kind"] == "TF+NZ" else ["NZ-low", "NZ-high"]
        maps = {}
        for kind in ("PD", "NE"):
            a_vals, b_vals, ys = slice_initial_conditions(kind, p, NG, NG, lo=1e-6, hi=0.9)
            df = basin_map(p, ys, ND, o)
            maps[kind] = (a_vals, b_vals, ys, df)
            df.insert(0, "case", cname); df.insert(1, "slice", kind)
            save_table(df, f"p3_basin_grid_{cname}_{kind}")
            lab = df.attractor.value_counts().to_dict()
            log(f"   {cname} {kind}: attractors {lab}, unresolved={(df.attractor=='unresolved').sum()}, failures={(~df.solver_ok).sum()}")
            below = df[~df.start_above]
            row = {"case": cname, "slice": kind, "n": len(df), "n_start_below_0.1": len(below)}
            for l in order:
                row[f"frac_{l}"] = float(np.mean(df.attractor == l))
                row[f"frac_{l}_given_start_below_0.1"] = float(np.mean(below.attractor == l)) if len(below) else float("nan")
            for th in case["thresholds"]:
                key = f"event240_{th:g}".replace(".", "p")
                row[f"event240_theta{th:g}_frac_of_valid_starts"] = float(np.nanmean(df[key])) if key in df else float("nan")
            row["event240_theta0.1_frac_of_valid_starts"] = float(np.nanmean(df.event240))
            summary.append(row)
        results[cname] = (p, maps, order, o)

        # ---- random 3-D initial conditions (P0, D0, E0)
        ys = random_initial_conditions(p, 3000, seeds["basins"], lo=1e-6, hi=0.9)
        dfr = basin_map(p, ys, ND, o)
        dfr.insert(0, "case", cname)
        save_table(dfr, f"p3_basin_random_ic_{cname}")
        for l in order + ["unresolved"]:
            m = dfr.attractor == l
            if m.sum():
                ens_rows.append({"case": cname, "attractor": l, "fraction_all_ics": float(m.mean()),
                                 "fraction_of_ics_starting_below_0.1": float(((dfr.attractor == l) & ~dfr.start_above).sum() / max((~dfr.start_above).sum(), 1)),
                                 "event240_theta0.1": float(np.nanmean(dfr.loc[m, "event240"])),
                                 **{f"event240_theta{th:g}": float(np.nanmean(dfr.loc[m, f"event240_{th:g}".replace('.', 'p')])) for th in case["thresholds"]}})
        # attractor vs event crosstab (below-threshold starts, theta = case's last threshold)
        for th in [0.1] + case["thresholds"]:
            key = "event240" if th == 0.1 else f"event240_{th:g}".replace(".", "p")
            below = dfr[dfr.N0 < th]
            for l in order:
                sub = below[below.attractor == l]
                cross.append({"case": cname, "threshold": th, "attractor": l, "n_starts_below_threshold": len(sub),
                              "event_by_240": int(np.nansum(sub[key])), "no_event_by_240": int(len(sub) - np.nansum(sub[key])),
                              "transient_exceedance_by_240": int(np.nansum(sub[f"transient240" if th == 0.1 else f"transient240_{th:g}".replace(".", "p")]))})

        # ---- basin boundary by bisection along lines
        eq = st[0] if case["kind"] == "TF+NZ" else st[0]
        high = st[-1] if case["kind"] != "TF+NZ" else st[0]
        small = np.array([1e-4 * 0.13, 1e-4 * 0.87, p.E_free])
        big = np.array([0.5 * 0.13, 0.5 * 0.87, p.E_free])
        pairs = {"N0: 1e-4 -> 0.5 (e0=1)": (small, big),
                 "e0: 1 -> 0.05 at N0=0.3": (np.array([0.3 * 0.13, 0.3 * 0.87, p.E_free]), np.array([0.3 * 0.13, 0.3 * 0.87, 0.05 * p.E_free])),
                 "N0: 1e-3 -> 0.9, e0=0.3": (np.array([1e-3 * 0.13, 1e-3 * 0.87, 0.3 * p.E_free]), np.array([0.9 * 0.13, 0.9 * 0.87, 0.3 * p.E_free]))}
        for nm, (ya, yb) in pairs.items():
            r = bisect_boundary(p, ya, yb, ND, o)
            row = {"case": cname, "line": nm, "found": r["found"], "label_a": r["label_a"], "label_b": r["label_b"]}
            if r["found"]:
                ys_ = r["y_star"]
                row.update(t_star=r["t_star"], bracket=r["bracket"], N0_at_boundary=ys_[0] + ys_[1], e0_at_boundary=ys_[2] / p.E_free,
                           saddle_N=[e.P + e.D for e in sd][0] if sd else float("nan"))
            lines.append(row)

        # ---- numerical verification of labels: boundary cells + random interior cells, tight options, jitter
        for kind in ("PD", "NE"):
            a_vals, b_vals, ys, df = maps[kind]
            codes = df.attractor.map({l: i for i, l in enumerate(order)}).fillna(len(order)).values.reshape(NG, NG).astype(int)
            bm = boundary_mask(codes).ravel()
            rng = np.random.default_rng(seeds["basins"])
            idx_b = np.where(bm)[0]
            idx_i = rng.choice(np.where(~bm)[0], size=min(200, int((~bm).sum())), replace=False)
            idx_b = idx_b[:300]
            sel = np.concatenate([idx_b, idx_i])
            tight = basin_map(p, [ys[i] for i in sel], ND, tight_options(o))
            agree = (tight.attractor.values == df.attractor.values[sel])
            jr_b = jitter_flip_rate(p, [ys[i] for i in idx_b], list(df.attractor.values[idx_b]), ND, o, 1e-6, seeds["basins"])
            jr_i = jitter_flip_rate(p, [ys[i] for i in idx_i], list(df.attractor.values[idx_i]), ND, o, 1e-6, seeds["basins"])
            checks.append({"case": cname, "slice": kind, "n_boundary_cells_checked": len(idx_b), "n_interior_cells_checked": len(idx_i),
                           "label_agreement_boundary_tight_numerics": float(agree[:len(idx_b)].mean()) if len(idx_b) else float("nan"),
                           "label_agreement_interior_tight_numerics": float(agree[len(idx_b):].mean()),
                           "jitter_1e-6_flip_rate_boundary": jr_b, "jitter_1e-6_flip_rate_interior": jr_i,
                           "n_unresolved_tight": int((tight.attractor == "unresolved").sum())})
            log(f"   {cname} {kind}: tight-numerics agreement boundary={checks[-1]['label_agreement_boundary_tight_numerics']:.3f} "
                f"interior={checks[-1]['label_agreement_interior_tight_numerics']:.3f}; 1e-6 jitter flips: boundary={jr_b:.3f}, interior={jr_i:.3f}")

    save_table(pd.DataFrame(summary), "p3_basin_fractions")
    save_table(pd.DataFrame(ens_rows), "p3_basin_random_ic_summary")
    save_table(pd.DataFrame(cross), "p3_basin_attractor_vs_event")
    save_table(pd.DataFrame(lines), "p3_basin_boundary_bisection")
    save_table(pd.DataFrame(checks), "p3_basin_numerical_checks")

    make_figures(results, pd.DataFrame(cross))
    log("done")


if __name__ == "__main__":
    main()
