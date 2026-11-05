"""Phase 5 (descriptive only): Kaplan-Meier curves, log-rank tests and landmark summaries for right-censored event-time data.

Used for the METABRIC RFS endpoint as exported by cBioPortal.  Nothing here fits or overlays a mechanistic model.  The RFS time
origin is undocumented (Phase 4), so every statement is "time since the (unknown) RFS origin".

Estimators
----------
* Kaplan-Meier product-limit with Greenwood variance, log-log (Kalbfleisch-Prentice) pointwise CI  (`km_curve`, `km_point`).
* Mantel-Cox log-rank test for k >= 2 groups, optionally administratively censored at a horizon (`logrank_test`), pairwise
  version with Holm step-down adjustment (`pairwise_logrank`, `holm`).
* 5-year landmark analysis: conditional cumulative incidence among patients still event-free and under observation at the
  landmark (`landmark_incidence`).
* Reliability filter inherited from the Phase 4 criteria (`reliable_horizon_mask`).

Not implemented (deliberately): Aalen-Johansen cumulative incidence.  The RFS event is a composite of relapse and disease-specific
death and the export does not separate the first event (Phase 4), so competing-risk curves cannot be formed honestly.
"""
from __future__ import annotations

import numpy as np
from scipy import stats

from .synthetic_validation import kaplan_meier, km_at


def administrative_censor(time, event, horizon):
    """Censor at ``horizon`` (events strictly after the horizon become censored at the horizon)."""
    time = np.asarray(time, float).copy(); event = np.asarray(event, int).copy()
    late = time > horizon
    time[late] = horizon; event[late] = 0
    return time, event


def _loglog_ci(s, sd, z):
    if s <= 0 or s >= 1 or sd == 0:
        return max(0.0, s - z * sd), min(1.0, s + z * sd)
    se_ll = sd / (s * abs(np.log(s)))
    return s ** np.exp(z * se_ll), s ** np.exp(-z * se_ll)  # (lo, hi) of S


def km_point(time, event, t, conf: float = 0.95):
    """Return dict with S(t), se, CI of S, cumulative incidence 1-S and its CI, number at risk (time >= t) and events <= t."""
    time = np.asarray(time, float); event = np.asarray(event, int)
    ut, S, se = kaplan_meier(time, event)
    s, sd = km_at(ut, S, se, t)
    z = stats.norm.ppf(0.5 + conf / 2)
    lo, hi = _loglog_ci(s, sd, z)
    return {"S": s, "se": sd, "S_lo": lo, "S_hi": hi, "cum_inc": 1 - s, "ci_lo": 1 - hi, "ci_hi": 1 - lo,
            "n_at_risk": int(np.sum(time >= t)), "n_events_by_horizon": int(np.sum((time <= t) & (event == 1)))}


def km_curve(time, event, conf: float = 0.95, t_max: float | None = None):
    """Step curve (t, S, S_lo, S_hi) with t starting at 0; truncated at t_max if given."""
    ut, S, se = kaplan_meier(time, event)
    z = stats.norm.ppf(0.5 + conf / 2)
    lo, hi = zip(*[_loglog_ci(s, d, z) for s, d in zip(S, se)]) if len(S) else ([], [])
    t = np.concatenate([[0.0], ut]); S = np.concatenate([[1.0], S])
    lo = np.concatenate([[1.0], lo]); hi = np.concatenate([[1.0], hi])
    if t_max is not None:
        keep = t <= t_max
        t, S, lo, hi = t[keep], S[keep], lo[keep], hi[keep]
    return t, S, lo, hi


def logrank_test(time, event, group, horizon: float | None = None):
    """Mantel-Cox log-rank test for k >= 2 groups.  Returns (chi2, df, p, observed, expected).

    Uses the hypergeometric variance with the full covariance matrix and drops the last group (the statistic is invariant to which).
    """
    time = np.asarray(time, float); event = np.asarray(event, int); group = np.asarray(group)
    if horizon is not None:
        time, event = administrative_censor(time, event, horizon)
    labels = list(dict.fromkeys(group.tolist()))
    k = len(labels)
    if k < 2:
        raise ValueError("need at least two groups")
    gi = np.array([labels.index(g) for g in group])
    ut = np.unique(time[event == 1])
    O = np.zeros(k); E = np.zeros(k); V = np.zeros((k, k))
    for u in ut:
        at = time >= u
        n_j = np.array([np.sum(at & (gi == j)) for j in range(k)], float)
        n = n_j.sum()
        dm = (time == u) & (event == 1)
        d = dm.sum()
        d_j = np.array([np.sum(dm & (gi == j)) for j in range(k)], float)
        O += d_j
        E += d * n_j / n
        if n > 1:
            c = d * (n - d) / (n - 1) / n ** 2
            V += c * (n * np.diag(n_j) - np.outer(n_j, n_j))
    U = (O - E)[:-1]; Vr = V[:-1, :-1]
    chi2 = float(U @ np.linalg.solve(Vr, U))
    df = k - 1
    return chi2, df, float(stats.chi2.sf(chi2, df)), O, E


def holm(pvals):
    """Holm step-down adjusted p-values (same order as input)."""
    p = np.asarray(pvals, float); m = len(p)
    order = np.argsort(p)
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * p[i])
        adj[i] = min(1.0, running)
    return adj


def pairwise_logrank(time, event, group, horizon: float | None = None):
    """All pairwise log-rank tests; returns list of dicts (a, b, chi2, p_raw) -- Holm applied by the caller (within a family)."""
    group = np.asarray(group)
    labels = list(dict.fromkeys(group.tolist()))
    rows = []
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            m = (group == labels[i]) | (group == labels[j])
            chi2, _, p, _, _ = logrank_test(np.asarray(time)[m], np.asarray(event)[m], group[m], horizon)
            rows.append({"group_a": labels[i], "group_b": labels[j], "chi2": chi2, "p_raw": p})
    if rows:
        adj = holm([r["p_raw"] for r in rows])
        for r, a in zip(rows, adj):
            r["p_holm"] = float(a)
    return rows


def landmark_incidence(time, event, landmark: float, ahead: float, conf: float = 0.95):
    """Among patients with time > landmark: KM cumulative incidence of the event within ``ahead`` months after the landmark."""
    time = np.asarray(time, float); event = np.asarray(event, int)
    m = time > landmark
    t2 = time[m] - landmark; e2 = event[m]
    out = km_point(t2, e2, ahead, conf)
    out["n_landmark"] = int(m.sum())
    return out


def reliable_horizon_mask(n_events, n_at_risk, horizon_years: int):
    """Phase 4 reliability criterion: 5 y and 10 y need >= 50 events (whole group) and >= 100 at risk; 20 y needs >= 50 events and >= 30 at risk."""
    need = {5: 100, 10: 100, 20: 30}
    if horizon_years not in need:
        raise ValueError("Phase 4 defined reliability only at 5, 10 and 20 years")
    return bool(n_events >= 50 and n_at_risk >= need[horizon_years])
