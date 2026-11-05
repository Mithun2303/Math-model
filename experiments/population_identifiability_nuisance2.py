"""Phase 4, Task 8 (pre-registered criterion B, part 2): same population test as population_identifiability_nuisance.py but with the
unmeasurable initial-burden scale s0 (N0 -> s0*N0) ALSO free, in addition to tau. Smoothed Monte Carlo pushforward through a 2-D table
T(ln beta, ln s0). SYNTHETIC. Run: .venv/bin/python experiments/population_identifiability_nuisance2.py
"""
from _p4common import *  # noqa: F401,F403
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize
from scipy.stats import norm

from dormancy.observation import detection_time, interval_censor, visit_schedule
from dormancy.uncertainty import parallel_map

SCALES = np.array([0.1, 0.3, 1.0, 3.0, 10.0])


def worker(args):
    th, b, s0 = args
    return detection_time({**th, "beta": float(b), "N0": th["N0"] * s0}, 0.1, 3000.0)


def main():
    seeds = CFG["seeds"]
    th = truth("G1_baseline")
    lb = np.linspace(np.log(0.012), np.log(0.115), 36)
    items = [(th, np.exp(b), s) for b in lb for s in SCALES]
    Tt = np.array(parallel_map(worker, items, workers=3, chunksize=4)).reshape(len(lb), len(SCALES))
    Tt = np.where(np.isfinite(Tt), Tt, 5000.0)
    interp = RegularGridInterpolator((lb, np.log(SCALES)), Tt, bounds_error=False, fill_value=None)
    log("2-D table done")
    rng = np.random.default_rng(seeds["populations"] + 1)
    Z = np.random.default_rng(12345).standard_normal(4000)
    h = 2.0
    m0, s0_, tau0, n = np.log(0.04), 0.3, 12.0, 2380
    visits = visit_schedule(12.0, 240.0)

    def model_T(m, s, ls0):
        lbeta = np.clip(m + s * Z, lb[0], lb[-1])
        pts = np.column_stack([lbeta, np.full_like(lbeta, ls0)])
        T = interp(pts)
        T[(m + s * Z) > lb[-1]] = 5000.0            # beta beyond the table: no event
        return T

    def nll(x, L, R):
        m, ls, tau, ls0 = x
        s = np.exp(ls)
        T = model_T(m, s, ls0)
        def F(t):
            t = np.atleast_1d(np.asarray(t, float))
            out = np.empty(len(t))
            fin = np.isfinite(t)
            out[~fin] = np.mean(T < 4999)
            if fin.any():
                out[fin] = np.mean(norm.cdf((t[fin, None] - tau - T[None, :]) / h), axis=1)
            return out
        ev = np.isfinite(R)
        FL = F(L); FR = F(np.where(ev, R, 0.0))
        p = np.where(ev, FR - FL, 1.0 - FL)
        return float(-np.sum(np.log(np.maximum(p, 1e-300))))

    # synthetic data: true N0 scale 1, tau 12; individual T from the s0 = 1 slice
    reps = []
    for rep in range(8):
        beta_l = m0 + s0_ * rng.standard_normal(n)
        Tm = np.where(beta_l <= lb[-1], RegularGridInterpolator((lb, np.log(SCALES)), Tt, bounds_error=False, fill_value=None)(np.column_stack([np.clip(beta_l, lb[0], lb[-1]), np.zeros(n)])), 5000.0)
        Tc = Tm + tau0
        L = np.empty(n); R = np.empty(n)
        for i, t in enumerate(Tc):
            l, r_, e = interval_censor(t if t < 4999 else np.inf, visits, 240.0)
            L[i], R[i] = l, r_
        # bin the data into unique (L, R) cells with counts (speed)
        cells, counts = np.unique(np.column_stack([L, R]), axis=0, return_counts=True)
        Lc, Rc = cells[:, 0], cells[:, 1]

        def nllw(x):
            m, ls, tau, ls0 = x
            if not (-8 < ls < 1) or not (-30 < tau < 120) or not (np.log(0.05) < ls0 < np.log(20)):
                return 1e12
            T = model_T(m, np.exp(ls), ls0)
            def F(t):
                t = np.atleast_1d(np.asarray(t, float)); out = np.empty(len(t)); fin = np.isfinite(t)
                out[~fin] = np.mean(T < 4999)
                out[fin] = np.mean(norm.cdf((t[fin, None] - tau - T[None, :]) / h), axis=1)
                return out
            ev = np.isfinite(Rc)
            FL = F(Lc); FR = F(np.where(ev, Rc, 0.0))
            p = np.where(ev, FR - FL, 1.0 - FL)
            return float(-np.sum(counts * np.log(np.maximum(p, 1e-300))))

        best = None
        for x0 in ([m0, np.log(s0_), tau0, 0.0], [m0 + 0.2, np.log(0.4), 0.0, 0.7], [m0 - 0.2, np.log(0.2), 25.0, -0.7]):
            r = minimize(nllw, x0, method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-5, "maxiter": 600})
            if best is None or r.fun < best.fun:
                best = r
        reps.append({"rep": rep, "m_hat": best.x[0], "s_hat": float(np.exp(best.x[1])), "tau_hat": best.x[2], "N0_scale_hat": float(np.exp(best.x[3])), "nll": best.fun})
        log(f"rep {rep}: m_hat {best.x[0]:.3f} s_hat {np.exp(best.x[1]):.3f} tau_hat {best.x[2]:.1f} N0scale_hat {np.exp(best.x[3]):.2f}")
        if rep == 0:
            prof = {}
            xb = best.x
            for name, idx, grid in (("m", 0, np.linspace(xb[0] - 0.4, xb[0] + 0.4, 17)), ("ln s", 1, np.linspace(xb[1] - 0.5, xb[1] + 0.5, 17)),
                                    ("tau", 2, np.linspace(-25, 70, 17)), ("ln N0 scale", 3, np.linspace(-1.5, 1.5, 17))):
                vals = []
                for gv in grid:
                    free = [i for i in range(4) if i != idx]
                    def f(z, gv=gv, idx=idx, free=free):
                        x = xb.copy(); x[idx] = gv; x[free] = z
                        return nllw(x)
                    r = minimize(f, xb[free], method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-4, "maxiter": 300})
                    vals.append(2 * (r.fun - best.fun))
                prof[name] = (grid, np.array(vals), xb[idx])
    df = pd.DataFrame(reps)
    save_table(df, "p4_population_nuisance2_replicates")
    rows = [{"quantity": "truth", "value": f"m={m0:.3f}, s={s0_}, tau={tau0}, N0 scale=1"},
            {"quantity": "sd across replicates: m_hat / s_hat / tau_hat / N0scale_hat", "value": f"{df.m_hat.std():.3f} / {df.s_hat.std():.3f} / {df.tau_hat.std():.1f} / {df.N0_scale_hat.std():.2f}"},
            {"quantity": "mean: m_hat / s_hat / tau_hat / N0scale_hat", "value": f"{df.m_hat.mean():.3f} / {df.s_hat.mean():.3f} / {df.tau_hat.mean():.1f} / {df.N0_scale_hat.mean():.2f}"},
            {"quantity": "corr(tau_hat, ln N0scale_hat) across replicates", "value": f"{np.corrcoef(df.tau_hat, np.log(df.N0_scale_hat))[0, 1]:.3f}"},
            {"quantity": "corr(m_hat, tau_hat) across replicates", "value": f"{np.corrcoef(df.m_hat, df.tau_hat)[0, 1]:.3f}"}]
    pc = []
    for name, (grid, vals, xb) in prof.items():
        inside = grid[vals <= 3.841]
        closed = bool(vals[0] > 3.841 and vals[-1] > 3.841)
        rows.append({"quantity": f"profile rep 0 (tau and N0 scale free): 95% interval for {name}", "value": f"[{inside.min():.3f}, {inside.max():.3f}]  closed within window: {closed}  (window [{grid[0]:.2f}, {grid[-1]:.2f}])"})
        for g_, v_ in zip(grid, vals):
            pc.append({"parameter": name, "grid": g_, "delta_chi2": v_, "estimate": xb})
    save_table(pd.DataFrame(rows), "p4_population_nuisance2_summary")
    save_table(pd.DataFrame(pc), "p4_population_nuisance2_profiles")
    print(pd.DataFrame(rows).to_string(index=False))
    log("done")


if __name__ == "__main__":
    main()
