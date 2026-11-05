"""Phase 4, Task 7 (+ fig 8): synthetic survival validation and population-level hyper-parameter recovery.

Run: .venv/bin/python -u experiments/survival_validation.py
SYNTHETIC: the 'true' population is the pooled prior-P1 Monte Carlo population of Phase 3 (hypothetical priors); synthetic cohorts
are drawn by resampling individuals and censoring them with KNOWN mechanisms. Outputs: results/tables/p4_survival_*.csv and
results/figures/fig_p4_08_synthetic_event_free_survival_recovery.png.
"""
import math

from _p4common import *  # noqa: F401,F403
import matplotlib.pyplot as plt
from scipy import stats
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import brentq, minimize, minimize_scalar

from dormancy.observation import detection_time, interval_censor, visit_schedule
from dormancy.synthetic_validation import (kaplan_meier, km_at, km_cumulative_incidence_ci, lognormal_beta_nll,
                                           population_T_table, turnbull)
from dormancy.uncertainty import parallel_map

HZ = (60, 120, 240)
GRID_T = np.arange(12.0, 241.0, 12.0)
NWORK = 3
NREP = 200
LOSS_MEAN = 500.0                       # exponential loss-to-follow-up mean (months): P(censored before 240) ~ 38%
TB_BOOT = 100
TB_BOOT_REPS = 50
SOLVER = dict(method="DOP853", rtol=1e-9, atol=1e-12)


# ----------------------------------------------------------------------------------------------- part 1: cohorts
def load_pool():
    names = ["p3_mc_P1_samples.csv", "p3_mc_P1_base_run_B_samples.csv", "p3_mc_P1_batches_run_C_samples.csv"]
    frames, used = [], []
    for nm in names:
        p = TAB / nm
        if p.exists():
            frames.append(pd.read_csv(p, comment="#", usecols=["event_0p1", "T_0p1", "transient_0p1", "persistent_0p1"]))
            used.append(nm)
    d = pd.concat(frames, ignore_index=True)
    T = np.where(d.event_0p1.values == 1, d.T_0p1.values, np.inf)
    pers = (d.event_0p1.values == 1) & (d.persistent_0p1.values == 1)
    return T, pers, used


def km_F(time, event, ts):
    ut, S, se = kaplan_meier(time, event)
    return np.array([1 - km_at(ut, S, se, t)[0] for t in ts])


def tb_F(L, R, ts):
    _, _, F = turnbull(L, R, max_iter=300, tol=1e-7)
    return np.asarray(F(np.asarray(ts, float)), float)


def naive_F(T_obs_event, n, ts):
    """Naive: censored individuals treated as event-free; F(t) = (# recorded events with time <= t) / n."""
    s = np.sort(T_obs_event[np.isfinite(T_obs_event)])
    return np.searchsorted(s, ts, side="right") / n


def cohort_rep(args):
    """One synthetic cohort: all schemes and estimators. Returns a list of result rows."""
    T_pool, pers_pool, w, n, seed, pop, do_boot, F_true, F_true_pers, ts_eval = args
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(T_pool), n, p=w)
    T = T_pool[idx]
    Tp = np.where(pers_pool[idx], T, np.inf)               # persistent-only event definition
    rows = []

    def add(scheme, est_name, est_at, truth_curve, h_list, ci=None, boot_ci=None):
        sup = float(np.max(np.abs(est_at["grid"] - truth_curve["grid"])))
        for h in h_list:
            row = {"population": pop, "n": n, "rep": seed, "scheme": scheme, "estimator": est_name, "horizon": h,
                   "estimate": float(est_at[h]), "truth": float(truth_curve[h]), "sup_norm": sup}
            if ci is not None:
                row["covered"] = bool(ci[h][0] <= truth_curve[h] <= ci[h][1])
            if boot_ci is not None:
                row["covered"] = bool(boot_ci[h][0] <= truth_curve[h] <= boot_ci[h][1])
            rows.append(row)

    def curve(F_fn, ts_grid):
        out = {"grid": np.asarray(F_fn(ts_grid))}
        for h in HZ:
            out[h] = float(F_fn(np.array([float(h)]))[0])
        return out

    tc = {"grid": F_true(GRID_T)}
    for h in HZ:
        tc[h] = float(F_true(np.array([float(h)]))[0])

    # (a) administrative right censoring at 60 / 120 / 240
    for H in HZ:
        time, ev = np.minimum(T, H), (T <= H).astype(int)
        hs = [h for h in HZ if h <= H]
        gr = GRID_T[GRID_T <= H]
        est = {"grid": km_F(time, ev, gr)}
        tcH = {"grid": tc["grid"][:len(gr)], **{h: tc[h] for h in hs}}
        for h in hs:
            est[h] = float(km_F(time, ev, np.array([float(h)]))[0])
        cis = {h: km_cumulative_incidence_ci(time, ev, float(h))[1:] for h in hs}
        add(f"admin censoring at {H} mo", "Kaplan-Meier", est, tcH, hs, ci=cis)
    # (b) administrative 240 + random loss to follow-up
    C = rng.exponential(LOSS_MEAN, n)
    time, ev = np.minimum(np.minimum(T, C), 240.0), ((T <= C) & (T <= 240)).astype(int)
    est = {"grid": km_F(time, ev, GRID_T)}
    for h in HZ:
        est[h] = float(km_F(time, ev, np.array([float(h)]))[0])
    cis = {h: km_cumulative_incidence_ci(time, ev, float(h))[1:] for h in HZ}
    add("admin 240 + random loss (exp, mean 500 mo)", "Kaplan-Meier", est, tc, HZ, ci=cis)
    Tn = np.where(ev == 1, T, np.inf)
    nest = {"grid": naive_F(Tn, n, GRID_T), **{h: float(naive_F(Tn, n, np.array([float(h)]))[0]) for h in HZ}}
    add("admin 240 + random loss (exp, mean 500 mo)", "naive (censored = event-free)", nest, tc, HZ)
    # (d) wrong event definition: analyst records only persistent exceedances, but the target is first crossing
    timep, evp = np.minimum(Tp, 240.0), (Tp <= 240).astype(int)
    estp = {"grid": km_F(timep, evp, GRID_T), **{h: float(km_F(timep, evp, np.array([float(h)]))[0]) for h in HZ}}
    tcp = {"grid": F_true_pers(GRID_T), **{h: float(F_true_pers(np.array([float(h)]))[0]) for h in HZ}}
    add("persistent-only events, target = persistent (right definition)", "Kaplan-Meier", estp, tcp, HZ)
    add("persistent-only events, target = first crossing (WRONG definition)", "Kaplan-Meier", estp, tc, HZ)
    # (c) interval censoring by visits
    for step, label in ((12.0, "annual visits"), (6.0, "6-monthly visits")):
        visits = visit_schedule(step, 240.0)
        LR = [interval_censor(float(t), visits, 240.0) for t in T]
        L = np.array([x[0] for x in LR]); R = np.array([x[1] for x in LR]); evi = np.array([x[2] for x in LR])
        est = {"grid": tb_F(L, R, GRID_T), **{h: float(tb_F(L, R, np.array([float(h)]))[0]) for h in HZ}}
        boot = None
        if do_boot and step == 12.0:
            Fb = np.empty((TB_BOOT, len(HZ)))
            for b in range(TB_BOOT):
                j = rng.integers(0, n, n)
                Fb[b] = tb_F(L[j], R[j], np.array(HZ, float))
            boot = {h: (float(np.quantile(Fb[:, k], 0.025)), float(np.quantile(Fb[:, k], 0.975))) for k, h in enumerate(HZ)}
        add(f"interval censoring, {label}", "Turnbull NPMLE", est, tc, HZ, boot_ci=boot)
        mid = np.where(evi == 1, 0.5 * (L + R), np.inf)
        mest = {"grid": naive_F(mid, n, GRID_T), **{h: float(naive_F(mid, n, np.array([float(h)]))[0]) for h in HZ}}
        add(f"interval censoring, {label}", "midpoint imputation (naive)", mest, tc, HZ)
    return rows


# ----------------------------------------------------------------------------------------------- part 2: hyper-parameters
def T_of_beta_factory(theta):
    return lambda b: detection_time({**theta, "beta": float(b)}, 0.1, 5000.0, **SOLVER)


def table_1d(theta, betas):
    f = T_of_beta_factory(theta)
    return population_T_table(f, betas)


def make_T_interp(be, Te):
    lb = np.log(be)
    return lambda beta_arr: np.interp(np.log(beta_arr), lb, Te, right=np.inf, left=Te[0])


def grouped_nll(params, rows, cnt, be, Te, bdet):
    s = 0.0
    for (L, R), c in zip(rows, cnt):
        s += c * lognormal_beta_nll(params, be, Te, np.array([L]), np.array([R]), bdet)
    return s


def hyper_rep(args):
    th_hat_tab, gen_tab, m_true, s_true, n, seed = args
    be, Te, bdet = th_hat_tab
    Tgen = make_T_interp(*gen_tab)
    rng = np.random.default_rng(seed)
    beta = np.exp(m_true + s_true * rng.standard_normal(n))
    T = Tgen(beta)
    visits = visit_schedule(12.0, 240.0)
    LR = [interval_censor(float(t), visits, 240.0) for t in T]
    uniq, cnt = np.unique(np.array(LR)[:, :2], axis=0, return_counts=True)
    rows = [(float(a), float(b)) for a, b in uniq]
    f = lambda p: grouped_nll(p, rows, cnt, be, Te, bdet)
    best = None
    for st in ((math.log(0.05), math.log(0.5)), (math.log(0.03), math.log(0.2))):
        r = minimize(f, st, method="Nelder-Mead", options={"xatol": 1e-5, "fatol": 1e-7, "maxiter": 400})
        if best is None or r.fun < best.fun:
            best = r
    m_hat, ls_hat = best.x
    # Wald from numerical Hessian
    h = 1e-3
    H = np.zeros((2, 2))
    x0 = best.x
    for i in range(2):
        for j in range(2):
            ei, ej = np.eye(2)[i] * h, np.eye(2)[j] * h
            H[i, j] = (f(x0 + ei + ej) - f(x0 + ei - ej) - f(x0 - ei + ej) + f(x0 - ei - ej)) / (4 * h * h)
    cov = np.linalg.pinv(H)
    sd = np.sqrt(np.maximum(np.diag(cov), 0))
    # profile CIs
    def profile(which, grid):
        out = []
        for g in grid:
            if which == "m":
                rr = minimize_scalar(lambda ls: f((g, ls)), bounds=(math.log(0.02), math.log(2.0)), method="bounded", options={"xatol": 1e-5})
            else:
                rr = minimize_scalar(lambda m: f((m, g)), bounds=(math.log(0.01), math.log(0.2)), method="bounded", options={"xatol": 1e-6})
            out.append(rr.fun)
        return 2 * (np.array(out) - best.fun)

    res = {"seed": seed, "m_hat": m_hat, "ln_s_hat": ls_hat, "sd_m": sd[0], "sd_ln_s": sd[1], "nll": best.fun}
    for which, ctr, sdv, truth in (("m", m_hat, sd[0], m_true), ("s", ls_hat, sd[1], math.log(s_true))):
        grid = ctr + np.linspace(-4, 4, 17) * max(sdv, 1e-3)
        d = np.maximum(profile(which, grid), 0)
        lo, hi = None, None
        i0 = int(np.argmin(np.abs(grid - ctr)))
        for i in range(i0, len(grid)):
            if d[i] > 3.841:
                hi = grid[i - 1] + (3.841 - d[i - 1]) / (d[i] - d[i - 1]) * (grid[i] - grid[i - 1]) if i > i0 else grid[i]
                break
        for i in range(i0, -1, -1):
            if d[i] > 3.841:
                lo = grid[i + 1] + (3.841 - d[i + 1]) / (d[i] - d[i + 1]) * (grid[i] - grid[i + 1]) if i < i0 else grid[i]
                break
        res[f"{which}_profile_ci"] = (lo, hi)
        res[f"{which}_profile_covered"] = bool(lo is not None and hi is not None and lo <= truth <= hi)
        res[f"{which}_wald_covered"] = bool(abs(ctr - truth) <= 1.96 * sdv)
        res[f"{which}_profile_grid"] = grid.tolist()
        res[f"{which}_profile_d"] = d.tolist()
    return res


def table_2d(theta, betas, rs):
    T = np.full((len(betas), len(rs)), 1e4)
    for i, b in enumerate(betas):
        for j, r_ in enumerate(rs):
            t = detection_time({**theta, "beta": float(b), "r": float(r_)}, 0.1, 5000.0, **SOLVER)
            if np.isfinite(t):
                T[i, j] = t
    return T


def table_row_worker(args):
    theta, b, rs = args
    out = []
    for r_ in rs:
        t = detection_time({**theta, "beta": float(b), "r": float(r_)}, 0.1, 5000.0, **SOLVER)
        out.append(t if np.isfinite(t) else 1e4)
    return out


def pushforward_nll_factory(interp, Z, rows, cnt, h_smooth=2.0):
    ts_u = np.unique(np.concatenate([np.array([a for a, _ in rows]), np.array([b for _, b in rows if np.isfinite(b)])]))
    idx = {t: i for i, t in enumerate(ts_u)}

    def nll(p):
        mb, lsb, mr, lsr = p
        sb, sr = np.exp(lsb), np.exp(lsr)
        beta = np.exp(mb + sb * Z[:, 0]); r_ = np.exp(mr + sr * Z[:, 1])
        pts = np.column_stack([np.clip(np.log(np.maximum(beta, 1e-300)), interp.grid[0][0], interp.grid[0][-1]),
                               np.clip(np.log(np.maximum(r_, 1e-300)), interp.grid[1][0], interp.grid[1][-1])])
        T = interp(pts)
        F = stats.norm.cdf((ts_u[:, None] - T[None, :]) / h_smooth).mean(1)
        s = 0.0
        for (L, R), c in zip(rows, cnt):
            pr = (F[idx[R]] - F[idx[L]]) if np.isfinite(R) else (1.0 - F[idx[L]])
            s -= c * math.log(max(pr, 1e-12))
        return s
    return nll, ts_u


def main():
    seeds = CFG["seeds"]
    T_pool, pers_pool, used = load_pool()
    N = len(T_pool)
    log(f"pool: {N} individuals from {used}")
    ev = np.isfinite(T_pool)
    Fp = lambda w: (lambda t: np.array([np.sum(w * (T_pool <= x)) / w.sum() for x in np.atleast_1d(t)]))
    # truth populations (weights): base and tilted mixtures
    pops = {"base": np.ones(N), "early-heavy (T<80 x3)": np.where(ev & (T_pool < 80), 3.0, 1.0),
            "late-heavy (160<T<=240 x3)": np.where(ev & (T_pool > 160), 3.0, 1.0), "mostly-no-crossing (no event x3)": np.where(~ev, 3.0, 1.0)}
    comp = pd.DataFrame([{"quantity": "individuals pooled", "value": N}, {"quantity": "runs used", "value": ", ".join(used)},
                         {"quantity": "no crossing within 240 mo", "value": int((~ev).sum())},
                         {"quantity": "early crossers (T < 80)", "value": int((ev & (T_pool < 80)).sum())},
                         {"quantity": "late crossers (160 < T <= 240)", "value": int((ev & (T_pool > 160)).sum())},
                         {"quantity": "events recorded as transient (cross then fall back within 240 mo)", "value": int((ev & ~pers_pool).sum())},
                         {"quantity": "F(60) / F(120) / F(240) first crossing", "value": " / ".join(f"{Fp(np.ones(N))(h)[0]:.4f}" for h in HZ)},
                         {"quantity": "F(60) / F(120) / F(240) persistent-only", "value": " / ".join(f"{(np.sum(pers_pool & (T_pool <= h)) / N):.4f}" for h in HZ)}])
    save_table(comp, "p4_survival_truth_population")

    # ---------------------------------------------------------------- cohorts
    items = []
    for pi, (pop, w) in enumerate(pops.items()):
        wn = w / w.sum()
        Ft = (lambda wn_: (lambda t: np.array([np.sum(wn_ * (T_pool <= x)) for x in np.atleast_1d(t)])))(wn)
        Ftp = (lambda wn_: (lambda t: np.array([np.sum(wn_ * (pers_pool & (T_pool <= x))) for x in np.atleast_1d(t)])))(wn)
        for n in (1000, 2380):
            nrep = NREP if pop == "base" else 100
            for rep in range(nrep):
                items.append((T_pool, pers_pool, wn, n, seeds["populations"] + 100000 * pi + 10000 * (n == 2380) + rep, pop,
                              bool(pop == "base" and n == 1000 and rep < TB_BOOT_REPS), Ft, Ftp, GRID_T))
    # closures are not picklable by spawn: rebuild truth functions inside the worker from weights instead
    items = [(a[0], a[1], a[2], a[3], a[4], a[5], a[6], None, None, a[9]) for a in items]
    log(f"{len(items)} synthetic cohorts")
    import sys
    if "--skip-cohorts" in sys.argv:
        summ = pd.read_csv(TAB / "p4_survival_estimator_performance.csv", comment="#")
        log("cohort part loaded from saved tables")
        df = None
    else:
        out = parallel_map(cohort_rep_wrapper, items, workers=NWORK, chunksize=2)
        df = pd.DataFrame([r for rows in out for r in rows])
        save_table(df, "p4_survival_cohort_estimates")
    g = df.groupby if df is not None else None
    if df is not None:
      g = df.groupby(["population", "n", "scheme", "estimator", "horizon"])
      summ = g.apply(lambda x: pd.Series({"truth": x.truth.iloc[0], "mean_estimate": x.estimate.mean(), "bias": (x.estimate - x.truth).mean(),
                                           "mc_se_bias": (x.estimate - x.truth).std(ddof=1) / math.sqrt(len(x)), "rmse": math.sqrt(((x.estimate - x.truth) ** 2).mean()),
                                           "mean_sup_norm": x.sup_norm.mean(), "n_rep": len(x),
                                           "coverage95": x.covered.mean() if "covered" in x and x.covered.notna().any() else np.nan,
                                           "n_cov_reps": int(x.covered.notna().sum()) if "covered" in x else 0}), include_groups=False).reset_index()
      save_table(summ, "p4_survival_estimator_performance")
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 20); pd.set_option("display.max_rows", 200)
    print(summ[(summ.population == "base") & (summ.n == 2380)][["scheme", "estimator", "horizon", "truth", "bias", "mc_se_bias", "rmse", "coverage95", "mean_sup_norm"]].round(4).to_string())
    print(summ[(summ.population != "base") & (summ.n == 2380) & (summ.estimator.isin(["Kaplan-Meier", "Turnbull NPMLE"]))][["population", "scheme", "estimator", "horizon", "bias", "rmse"]].round(4).to_string())

    # ---------------------------------------------------------------- hyper-parameters (1-D)
    th = truth("G1_baseline")
    m_true, s_true = math.log(0.04), 0.3
    be_fit, Te_fit = table_1d(th, np.geomspace(0.008, 0.2, 70))
    be_gen, Te_gen = table_1d(th, np.geomspace(0.008, 0.2, 300))
    log(f"T(beta) tables: fit {len(be_fit)} points, gen {len(be_gen)}; T range {Te_fit[0]:.1f}..{Te_fit[-1]:.0f}")
    # monotonicity check
    mono = bool(np.all(np.diff(Te_gen) > 0))
    hitems = [((be_fit, Te_fit, be_fit[-1]), (be_gen, Te_gen), m_true, s_true, 2380, seeds["populations"] + 900000 + r) for r in range(50)]
    hres = parallel_map(hyper_rep, hitems, workers=NWORK, chunksize=1)
    hdf = pd.DataFrame([{k: v for k, v in r.items() if not k.endswith(("_grid", "_d"))} for r in hres])
    save_table(hdf.astype({"m_profile_ci": str, "s_profile_ci": str}), "p4_survival_hyperparam_replicates")
    hs = pd.DataFrame([{"parameter": "m = mean ln(beta)", "truth": m_true, "mean_estimate": hdf.m_hat.mean(), "bias": (hdf.m_hat - m_true).mean(), "rmse": math.sqrt(((hdf.m_hat - m_true) ** 2).mean()),
                        "wald_coverage": hdf.m_wald_covered.mean(), "profile_coverage": hdf.m_profile_covered.mean(), "median_sd_wald": hdf.sd_m.median()},
                       {"parameter": "ln s (sd of ln beta)", "truth": math.log(s_true), "mean_estimate": hdf.ln_s_hat.mean(), "bias": (hdf.ln_s_hat - math.log(s_true)).mean(),
                        "rmse": math.sqrt(((hdf.ln_s_hat - math.log(s_true)) ** 2).mean()), "wald_coverage": hdf.s_wald_covered.mean(), "profile_coverage": hdf.s_profile_covered.mean(),
                        "median_sd_wald": hdf.sd_ln_s.median()}])
    hs["T_of_beta_monotone_increasing"] = mono
    hs["n_replicates"] = len(hdf)
    save_table(hs, "p4_survival_hyperparam_summary")
    print(hs.round(4).to_string())
    prof_rows = []
    for which in ("m", "s"):
        r0 = hres[0]
        for gx, d in zip(r0[f"{which}_profile_grid"], r0[f"{which}_profile_d"]):
            prof_rows.append({"parameter": which, "grid": gx, "delta_chi2": d, "truth": m_true if which == "m" else math.log(s_true), "seed": r0["seed"]})

    # ---------------------------------------------------------------- 2-D: r also varies -> non-identifiability of m_beta
    betas2 = np.geomspace(0.008, 0.14, 24); rs2 = np.geomspace(0.06, 0.22, 16)
    betas2g = np.geomspace(0.008, 0.14, 48); rs2g = np.geomspace(0.06, 0.22, 32)
    tab = lambda bs, rs: np.array(parallel_map(table_row_worker, [(th, b, rs) for b in bs], workers=NWORK, chunksize=2))
    T2f, T2g = tab(betas2, rs2), tab(betas2g, rs2g)
    I_fit = RegularGridInterpolator((np.log(betas2), np.log(rs2)), T2f, bounds_error=False, fill_value=1e4)
    I_gen = RegularGridInterpolator((np.log(betas2g), np.log(rs2g)), T2g, bounds_error=False, fill_value=1e4)
    mr_true, sr_true = math.log(0.12), 0.15
    rng = np.random.default_rng(seeds["populations"] + 777)
    n = 2380
    Zg = rng.standard_normal((n, 2))
    Tdat = I_gen(np.column_stack([np.clip(m_true + s_true * Zg[:, 0], np.log(betas2g[0]), np.log(betas2g[-1])), np.clip(mr_true + sr_true * Zg[:, 1], np.log(rs2g[0]), np.log(rs2g[-1]))]))
    visits = visit_schedule(12.0, 240.0)
    LR = [interval_censor(float(t), visits, 240.0) for t in Tdat]
    uq, cnt = np.unique(np.array(LR)[:, :2], axis=0, return_counts=True)
    rows = [(float(a), float(b)) for a, b in uq]
    Zfit = np.random.default_rng(seeds["populations"] + 888).standard_normal((20000, 2))
    nll, ts_u = pushforward_nll_factory(I_fit, Zfit, rows, cnt)
    starts = [(math.log(0.05), math.log(0.4), math.log(0.1), math.log(0.3)), (math.log(0.03), math.log(0.2), math.log(0.15), math.log(0.1)),
              (math.log(0.04), math.log(0.3), math.log(0.12), math.log(0.15))]
    fits = [minimize(nll, s0, method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 600}) for s0 in starts]
    best = min(fits, key=lambda r: r.fun)
    log(f"2-D fit: best nll {best.fun:.3f}; params {np.round(best.x, 3)} (truth {np.round([m_true, math.log(s_true), mr_true, math.log(sr_true)], 3)}); all nll {[round(f.fun, 2) for f in fits]}")
    # profile of mb with (ln sb, mr, ln sr) free
    grid2 = m_true + np.linspace(-0.8, 0.8, 17)
    prof2, ridge_x = [], []
    x_prev = best.x.copy()
    order = np.argsort(np.abs(grid2 - best.x[0]))
    pv = {}
    for i in order:
        g_ = grid2[i]
        cand = []
        for s0 in (x_prev[1:], np.array([math.log(0.3), math.log(0.12), math.log(0.15)])):
            r = minimize(lambda q: nll(np.concatenate([[g_], q])), s0, method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 500})
            cand.append(r)
        b = min(cand, key=lambda r: r.fun)
        x_prev = np.concatenate([[g_], b.x]); pv[i] = (2 * (b.fun - best.fun), x_prev.copy())
    for i in range(len(grid2)):
        prof2.append(pv[i][0]); ridge_x.append(pv[i][1])
    prof2 = np.maximum(np.array(prof2), 0)
    # fitted CDF along the ridge vs truth: event-time distribution is well estimated while m_beta is not
    tt = np.arange(12.0, 241.0, 12.0)
    def F_from(p, Z=Zfit, interp=I_fit):
        mb, lsb, mr, lsr = p
        beta = np.exp(np.clip(mb + np.exp(lsb) * Z[:, 0], np.log(betas2[0]), np.log(betas2[-1]))); r_ = np.exp(np.clip(mr + np.exp(lsr) * Z[:, 1], np.log(rs2[0]), np.log(rs2[-1])))
        T = interp(np.column_stack([np.log(beta), np.log(r_)]))
        return stats.norm.cdf((tt[:, None] - T[None, :]) / 2.0).mean(1)
    F_best = F_from(best.x)
    F_true2 = F_from(np.array([m_true, math.log(s_true), mr_true, math.log(sr_true)]), Z=Zfit, interp=I_gen) if False else None
    Tt = I_gen(np.column_stack([np.clip(m_true + s_true * Zfit[:, 0], np.log(betas2g[0]), np.log(betas2g[-1])), np.clip(mr_true + sr_true * Zfit[:, 1], np.log(rs2g[0]), np.log(rs2g[-1]))]))
    F_true2 = stats.norm.cdf((tt[:, None] - Tt[None, :]) / 2.0).mean(1)
    sup_best = float(np.max(np.abs(F_best - F_true2)))
    ridge_F = np.array([F_from(x) for x in ridge_x])
    ridge_sup_vs_best = ridge_F.__sub__(F_best[None, :])
    ok = prof2 < 3.841
    ridge_spread = float(np.max(np.abs(ridge_F[ok] - F_best[None, :]))) if ok.any() else float("nan")
    m_ok = grid2[ok]
    ut, S, se = kaplan_meier(*(lambda L, R: (np.where(np.isfinite(R), (L + R) / 2, L), np.isfinite(R).astype(int)))(np.array(LR)[:, 0], np.array(LR)[:, 1]))
    prof_df = pd.DataFrame({"m_beta_grid": grid2, "delta_chi2_2D_model(r varies; 4 unknowns)": prof2, "within_95pct": ok})
    save_table(prof_df, "p4_survival_hyperparam_profile_2D")
    save_table(pd.DataFrame(prof_rows), "p4_survival_hyperparam_profile_1D")
    summ2 = pd.DataFrame([{"quantity": "best-fit nll", "value": best.fun}, {"quantity": "estimated (m_beta, ln s_beta, m_r, ln s_r)", "value": str(np.round(best.x, 3))},
                          {"quantity": "truth", "value": str(np.round([m_true, math.log(s_true), mr_true, math.log(sr_true)], 3))},
                          {"quantity": "range of m_beta with profile delta-chi2 < 3.84 (grid truth+-0.8)", "value": f"{m_ok.min():.3f} .. {m_ok.max():.3f}" if ok.any() else "none"},
                          {"quantity": "n grid points within 95% (of 17) -> ridge", "value": int(ok.sum())},
                          {"quantity": "sup-norm |F_hat - F_true| for the event-time CDF at the best fit", "value": sup_best},
                          {"quantity": "max |F along the 95% ridge - F_best| (event-time CDF barely moves along the ridge)", "value": ridge_spread},
                          {"quantity": "same 1-D (r known) profile: n points within 95% of 17 (window +-4 sd)", "value": int(np.sum(np.array(hres[0]["m_profile_d"]) < 3.841))}])
    save_table(summ2, "p4_survival_hyperparam_2D_summary")
    print(summ2.to_string())

    # ---------------------------------------------------------------- figure 8
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax = axes[0, 0]
    bw = pops["base"] / pops["base"].sum()
    tt_f = np.linspace(0, 240, 241)
    Ftr = np.array([np.mean(T_pool <= x) for x in tt_f])
    rng = np.random.default_rng(seeds["populations"] + 4242)
    idx = rng.choice(N, 2380); Tc = T_pool[idx]
    C = rng.exponential(LOSS_MEAN, 2380)
    time, evn = np.minimum(np.minimum(Tc, C), 240.0), ((Tc <= C) & (Tc <= 240)).astype(int)
    ut, S, se = kaplan_meier(time, evn)
    Fk = np.array([1 - km_at(ut, S, se, x)[0] for x in tt_f]); sek = np.array([km_at(ut, S, se, x)[1] for x in tt_f])
    ax.plot(tt_f, Ftr, color=INK, lw=2.2, label="true F(t) (pooled population)")
    ax.step(tt_f, Fk, color=BLUE, lw=1.6, where="post", label="Kaplan-Meier (admin 240 + random loss)")
    ax.fill_between(tt_f, np.clip(Fk - 1.96 * sek, 0, 1), np.clip(Fk + 1.96 * sek, 0, 1), color=BLUE, alpha=0.2, step="post")
    ev_naive = np.where(evn == 1, Tc, np.inf)
    ax.plot(tt_f, naive_F(ev_naive, 2380, tt_f), color=ORANGE, ls="--", label="naive (censored = event-free): biased low")
    visits = visit_schedule(12.0, 240.0)
    LRc = [interval_censor(float(t), visits, 240.0) for t in Tc]
    Lc = np.array([x[0] for x in LRc]); Rc = np.array([x[1] for x in LRc]); evc = np.array([x[2] for x in LRc])
    ax.plot(tt_f, tb_F(Lc, Rc, tt_f), color=VIOLET, lw=1.5, label="Turnbull (annual visits)")
    ax.plot(tt_f, naive_F(np.where(evc == 1, 0.5 * (Lc + Rc), np.inf), 2380, tt_f), color=AQUA, ls=":", label="midpoint imputation (naive)")
    ax.set_xlabel("time [months]"); ax.set_ylabel("cumulative event probability F(t)"); ax.legend(fontsize=7, loc="upper left")
    ax.set_title("(a) One synthetic cohort (n = 2380): estimators vs truth", fontsize=9)
    ax = axes[0, 1]
    sub = summ[(summ.population == "base") & (summ.n == 2380) & (summ.horizon.isin([120, 240]))]
    combos = [("admin censoring at 240 mo", "Kaplan-Meier", "KM, admin"), ("admin 240 + random loss (exp, mean 500 mo)", "Kaplan-Meier", "KM, random loss"),
              ("admin 240 + random loss (exp, mean 500 mo)", "naive (censored = event-free)", "naive, random loss"), ("interval censoring, annual visits", "Turnbull NPMLE", "Turnbull, annual"),
              ("interval censoring, 6-monthly visits", "Turnbull NPMLE", "Turnbull, 6-mo"), ("interval censoring, annual visits", "midpoint imputation (naive)", "midpoint, annual"),
              ("persistent-only events, target = first crossing (WRONG definition)", "Kaplan-Meier", "KM persistent-only\nvs first crossing")]
    xs = np.arange(len(combos)); w = 0.38
    for k, (h, col) in enumerate(((120, BLUE), (240, ORANGE))):
        vals, errs = [], []
        for sc, es, _ in combos:
            r = sub[(sub.scheme == sc) & (sub.estimator == es) & (sub.horizon == h)]
            vals.append(float(r.bias.iloc[0]) if len(r) else np.nan); errs.append(1.96 * float(r.mc_se_bias.iloc[0]) if len(r) else 0)
        ax.bar(xs + (k - 0.5) * w, vals, w, yerr=errs, color=col, label=f"horizon {h} mo", capsize=2)
    ax.axhline(0, color=INK, lw=0.8); ax.set_xticks(xs); ax.set_xticklabels([c[2] for c in combos], rotation=35, ha="right", fontsize=7)
    ax.set_ylabel("bias of F(h) estimate (mean over 200 cohorts; 95% MC band)"); ax.legend(fontsize=7)
    ax.set_title("(b) Bias by censoring scheme and estimator (n = 2380)", fontsize=9)
    ax = axes[1, 0]
    r0 = hres[0]
    for which, col, tr, lab in (("m", BLUE, m_true, "m = mean ln(beta)"), ("s", ORANGE, math.log(s_true), "ln s")):
        gx = np.array(r0[f"{which}_profile_grid"]) - tr
        ax.plot(gx, r0[f"{which}_profile_d"], "o-", ms=3, color=col, label=f"{lab} (one data set)")
    ax.axhline(3.841, color=MAGENTA, ls="--", lw=0.8); ax.axvline(0, color=GREY, lw=0.8)
    ax.set_ylim(-0.5, 14); ax.set_xlabel("estimate - truth (log scale)"); ax.set_ylabel("profile $\\Delta\\chi^2$")
    ax.set_title(f"(c) One varying parameter: ln beta ~ N(m, s^2) is identifiable from interval-censored events\n(50 replicates: profile coverage m {hs.profile_coverage.iloc[0]:.2f}, ln s {hs.profile_coverage.iloc[1]:.2f})", fontsize=8.5)
    ax.legend(fontsize=7)
    ax = axes[1, 1]
    gx1 = np.array(r0["m_profile_grid"]) - m_true
    ax.plot(gx1, r0["m_profile_d"], "o-", ms=3, color=BLUE, label="r known: profile of m_beta")
    ax.plot(grid2 - m_true, prof2, "s-", ms=3, color=VIOLET, label="r also varies (4 unknowns): profile of m_beta")
    ax.axhline(3.841, color=MAGENTA, ls="--", lw=0.8); ax.axvline(0, color=GREY, lw=0.8); ax.set_ylim(-0.5, 14)
    ax.set_xlabel("m_beta - truth (log scale)"); ax.set_ylabel("profile $\\Delta\\chi^2$")
    ax.set_title(f"(d) Two varying parameters: m_beta is a ridge (non-identifiable) while the event-time CDF is fitted well\n(sup-norm error of F: {sup_best:.3f}; max change of F along the 95% ridge: {ridge_spread:.3f})", fontsize=8.5)
    ax.legend(fontsize=7)
    fig.suptitle("Synthetic survival validation: known censoring is handled by KM/Turnbull; event times identify a distribution, not mechanisms", y=0.995, fontsize=10)
    fig.tight_layout(rect=(0, 0.05, 1, 0.97)); p4_footer(fig, f"Truth: pooled prior-P1 population (n = {N}); N_detect 0.1; horizon 240 mo; hyper-parameter part uses G1 baseline with ln beta ~ N(ln 0.04, 0.3^2).")
    fig.savefig(FIG / "fig_p4_08_synthetic_event_free_survival_recovery.png"); plt.close(fig)
    log("done")


# worker that rebuilds the truth CDFs (closures are not picklable under spawn)
def cohort_rep_wrapper(args):
    T_pool, pers_pool, w, n, seed, pop, do_boot, _a, _b, ts_eval = args
    srt = np.argsort(T_pool)
    cw = np.cumsum(w[srt]); cwp = np.cumsum((w * pers_pool)[srt]); Ts = T_pool[srt]

    def make(c):
        return lambda t: np.array([c[np.searchsorted(Ts, x, side="right") - 1] if np.searchsorted(Ts, x, side="right") > 0 else 0.0 for x in np.atleast_1d(t)])
    return cohort_rep((T_pool, pers_pool, w, n, seed, pop, do_boot, make(cw), make(cwp), ts_eval))


if __name__ == "__main__":
    main()
