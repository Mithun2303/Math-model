"""Phase 4, Task 8 (pre-registered criterion B): can a population parameter combination be identified from INTERVAL-CENSORED event
times when the time-origin offset tau is left FREE as a nuisance? (Threshold / N0 nuisances are exactly confounded with tau in the
linear regime, see p4_event_only_N0_time_origin_confounding.csv, so tau stands for all three.)

Synthetic population: ln(beta) ~ N(m, s^2), tau_true = 12 months, annual visits in CLINICAL time to 240 months, n = 2380.
Run: .venv/bin/python experiments/population_identifiability_nuisance.py       SYNTHETIC.
"""
from _p4common import *  # noqa: F401,F403
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import norm

from dormancy.observation import detection_time, interval_censor, visit_schedule
from dormancy.uncertainty import parallel_map


def T_of_beta(args):
    th, b = args
    return detection_time({**th, "beta": float(b)}, 0.1, 3000.0)


def make_F(betas_ev, T_ev, beta_det, m, s, tau):
    """Vectorised model CDF of the clinical event time T_clin = tau + T(beta), ln beta ~ N(m, s^2); no event for beta > beta_det."""
    lt = np.log(betas_ev)
    top = norm.cdf((np.log(beta_det) - m) / s)

    def F(t):
        t = np.atleast_1d(np.asarray(t, float)) - tau
        out = np.empty(len(t))
        fin = np.isfinite(t)
        x = np.where(fin, t, T_ev[-1])
        out[:] = norm.cdf((np.interp(x, T_ev, lt) - m) / s)
        out[~fin] = top
        out[x <= T_ev[0]] = 0.0
        out[x >= T_ev[-1]] = top
        return out

    return F


def nll(params, tab, L, R):
    m, ls, tau = params
    s = np.exp(ls)
    F = make_F(*tab, m, s, tau)
    FL, FR = F(L), F(np.where(np.isfinite(R), R, np.inf))
    ev = np.isfinite(R)
    p = np.where(ev, FR - FL, 1.0 - FL)
    return float(-np.sum(np.log(np.maximum(p, 1e-300))))


def fit(tab, L, R, starts, fix=None):
    """Minimise over (m, ln s, tau); ``fix`` = {index: value} holds components fixed (profiling)."""
    best = None
    for x0 in starts:
        x0 = np.array(x0, float)
        free = [i for i in range(3) if not (fix and i in fix)]

        def f(z):
            x = x0.copy()
            if fix:
                for i, v in fix.items():
                    x[i] = v
            x[free] = z
            if not (-8 < x[1] < 1) or not (-30 < x[2] < 120):
                return 1e12
            return nll(x, tab, L, R)

        r = minimize(f, x0[free], method="Nelder-Mead", options={"xatol": 1e-5, "fatol": 1e-6, "maxiter": 800})
        if best is None or r.fun < best[0]:
            xx = x0.copy()
            if fix:
                for i, v in fix.items():
                    xx[i] = v
            xx[free] = r.x
            best = (r.fun, xx)
    return best


def main():
    seeds = CFG["seeds"]
    th = truth("G1_baseline")
    betas = np.geomspace(0.012, 0.11, 60)
    Ts = parallel_map(T_of_beta, [(th, b) for b in betas], workers=3, chunksize=2)
    Ts = np.array(Ts)
    ok = np.isfinite(Ts)
    beta_det = float(betas[ok][-1])
    tab = (betas[ok], Ts[ok], beta_det)
    log(f"T(beta) table: {ok.sum()} finite points; beta_det ~ {beta_det:.4f}; monotone increasing: {bool(np.all(np.diff(Ts[ok]) > 0))}")
    m0, s0, tau0, n = np.log(0.04), 0.3, 12.0, 2380
    visits = visit_schedule(12.0, 240.0)
    rng = np.random.default_rng(seeds["populations"])
    reps = []
    prof_store = None
    for rep in range(30):
        beta = np.exp(m0 + s0 * rng.standard_normal(n))
        Tm = np.array([np.interp(b, tab[0], tab[1], right=np.inf) if b <= tab[0][-1] else np.inf for b in beta])
        Tm = np.where(beta < tab[0][0], np.inf, Tm)
        Tc = Tm + tau0
        L = np.empty(n); R = np.empty(n)
        for i, t in enumerate(Tc):
            l, r_, e = interval_censor(t, visits, 240.0)
            L[i], R[i] = l, r_
        starts = [[m0, np.log(s0), tau0], [m0 + 0.3, np.log(0.5), 0.0], [m0 - 0.3, np.log(0.15), 30.0]]
        f_free = fit(tab, L, R, starts)
        f_tau0 = fit(tab, L, R, [[m0, np.log(s0), 0.0]], fix={2: 0.0})          # analyst assumes tau = 0 (wrong)
        reps.append({"rep": rep, "m_hat_tau_free": f_free[1][0], "s_hat_tau_free": float(np.exp(f_free[1][1])), "tau_hat": f_free[1][2], "nll_tau_free": f_free[0],
                     "m_hat_tau_fixed0": f_tau0[1][0], "s_hat_tau_fixed0": float(np.exp(f_tau0[1][1])), "nll_tau_fixed0": f_tau0[0],
                     "m_hat_minus_truth_tau_free": f_free[1][0] - m0, "s_hat_minus_truth_tau_free": float(np.exp(f_free[1][1])) - s0,
                     "m_hat_minus_truth_tau_fixed0": f_tau0[1][0] - m0, "s_hat_minus_truth_tau_fixed0": float(np.exp(f_tau0[1][1])) - s0})
        if rep == 0:
            prof = {}
            xb = f_free[1]
            for name, idx, grid in (("m", 0, np.linspace(xb[0] - 0.5, xb[0] + 0.5, 21)), ("ln s", 1, np.linspace(xb[1] - 0.6, xb[1] + 0.6, 21)), ("tau", 2, np.linspace(-20, 60, 21))):
                vals = []
                x_start = xb.copy()
                for gv in grid:
                    st = [x_start]
                    r_ = fit(tab, L, R, st, fix={idx: gv})
                    vals.append(2 * (r_[0] - f_free[0]))
                prof[name] = (grid, np.array(vals), xb[idx])
            prof_store = prof
        log(f"rep {rep}: tau_free m_hat {f_free[1][0]:.3f} (truth {m0:.3f}) s_hat {np.exp(f_free[1][1]):.3f} tau_hat {f_free[1][2]:.1f}; tau=0 m_hat {f_tau0[1][0]:.3f}")
    df = pd.DataFrame(reps)
    save_table(df, "p4_population_nuisance_replicates")
    summ = {"truth_m": m0, "truth_s": s0, "truth_tau": tau0, "n": n,
            "tau_free: sd(m_hat)": float(df.m_hat_tau_free.std()), "tau_free: sd(s_hat)": float(df.s_hat_tau_free.std()), "tau_free: sd(tau_hat)": float(df.tau_hat.std()),
            "tau_free: mean bias m": float(df.m_hat_minus_truth_tau_free.mean()), "tau_free: mean bias s": float(df.s_hat_minus_truth_tau_free.mean()),
            "tau_free: mean tau_hat": float(df.tau_hat.mean()),
            "tau=0 assumed: mean bias m": float(df.m_hat_minus_truth_tau_fixed0.mean()), "tau=0 assumed: mean bias s": float(df.s_hat_minus_truth_tau_fixed0.mean()),
            "tau=0 assumed: mean delta nll vs tau-free (2x)": float((2 * (df.nll_tau_fixed0 - df.nll_tau_free)).mean())}
    crit = {}
    for name, (grid, vals, xb) in prof_store.items():
        inside = grid[vals <= 3.841]
        crit[name] = (float(inside.min()) if len(inside) else None, float(inside.max()) if len(inside) else None, bool(vals[0] > 3.841 and vals[-1] > 3.841))
    rows = [{"quantity": k, "value": v} for k, v in summ.items()]
    for name, (lo, hi, closed) in crit.items():
        rows.append({"quantity": f"profile (rep 0, tau free): 95% interval for {name}", "value": f"[{lo}, {hi}]  closed within the window: {closed}"})
    save_table(pd.DataFrame(rows), "p4_population_nuisance_summary")
    pc = []
    for name, (grid, vals, xb) in prof_store.items():
        for g_, v_ in zip(grid, vals):
            pc.append({"parameter": name, "grid": g_, "delta_chi2": v_, "estimate": xb})
    save_table(pd.DataFrame(pc), "p4_population_nuisance_profiles")
    print(pd.DataFrame(rows).to_string())
    log("done")


if __name__ == "__main__":
    main()
