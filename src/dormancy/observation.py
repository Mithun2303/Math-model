"""Phase 4: explicit observation model.

Five levels are kept apart (nothing here is a validated clinical measurement model):

1. **Hidden biological states**: x(t) = (P, D, e) from the P-D-e ODE (`HiddenTrajectory`).
2. **Mathematically defined tumour burden**: N(t) = P(t) + D(t), a fraction of K.
3. **Hypothetical measurement signals**: y_N = c * N * exp(sigma * eps) (multiplicative log-normal noise, unknown/known scale c),
   y_P, y_D likewise, y_e = c_e * e * exp(sigma * eps); values below a limit of detection (LOD) are LEFT-CENSORED, not dropped.
4. **Threshold-defined recurrence event**: first upward crossing of N = N_detect (the project's event definition), observed only
   through the visit schedule, i.e. INTERVAL-censored (`interval_censor`).
5. **Observed clinical outcome**: a further mapping (time origin, detection delay, competing death) in `endpoint_mapping.py`.

Relapse is diagnosed (and treated) once N reaches N_detect, so the UNTREATED course is not observed after that time:
``truncate_at_detection`` stops the observation window there (long follow-up means observing up to detection, not N(240)).
Noise, LOD, schedules and signal scales are assumptions.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.integrate import solve_ivp

from .normalized import NormalizedParams
from .uncertainty import run_trajectory

CHANNELS = {"N": ["N"], "PD": ["P", "D"], "Ne": ["N", "e"], "PDe": ["P", "D", "e"]}


@dataclass(frozen=True)
class ObservationAssumptions:
    """Hypothetical observation assumptions (all are assumptions, see configs/phase4.yaml)."""
    observed: str = "N"                  # key of CHANNELS
    signal_scale_c: float = 1.0          # signal = c * N (also applied to P, D)
    immune_scale_ce: float = 1.0         # immune signal = c_e * e
    noise_sigma_ln: float = 0.10         # multiplicative log-normal noise, natural-log sd
    lod: float = 1e-6                    # limit of detection of the (scaled) signal
    schedule_months: float = 6.0         # regular visits
    follow_up_months: float = 240.0
    detection_threshold: float = 0.1     # N_detect: relapse diagnosed -> untreated observation stops
    truncate_at_detection: bool = True
    include_baseline: bool = True        # measurement at t = 0 (post-treatment baseline)

    def with_(self, **kw):
        return replace(self, **kw)


def theta_to_model(theta: dict, mu: float = 0.05, h: float = 0.02, K: float = 1.0):
    """(Params in the s = mu gauge, initial state) from a normalised-parameter dict."""
    n = NormalizedParams(K=K, r=theta["r"], q=theta["q"], a=theta["a"], delta_D=theta["delta_D"], beta=theta["beta"],
                         mu=mu, rho=theta["rho"], h=h, eta=theta["eta"])
    p = n.to_params()
    N0, fD, e0 = theta["N0"], theta["fD"], theta.get("e0", 1.0)
    return p, np.array([N0 * (1 - fD), N0 * fD, e0 * p.E_free])


@dataclass
class HiddenTrajectory:
    """Hidden states of the ground truth, generated with a TIGHT, DIFFERENT solver than the fitter (no inverse crime)."""
    times: np.ndarray
    X: np.ndarray                         # (T, 3): P, D, e

    @property
    def N(self):
        return self.X[:, 0] + self.X[:, 1]


def simulate_hidden(theta: dict, times, method: str = "Radau", rtol: float = 1e-12, atol: float = 1e-15, mu: float = 0.05,
                    h: float = 0.02, K: float = 1.0) -> HiddenTrajectory:
    from .model import rhs
    p, x0 = theta_to_model(theta, mu, h, K)
    times = np.asarray(times, float)
    sol = solve_ivp(rhs, (0.0, float(times.max())), x0, args=(p,), method=method, rtol=rtol, atol=atol, t_eval=times)
    return HiddenTrajectory(sol.t, sol.y.T)


def detection_time(theta: dict, N_detect: float, horizon: float, method: str = "Radau", rtol: float = 1e-12, atol: float = 1e-15,
                   mu: float = 0.05, h: float = 0.02, K: float = 1.0) -> float:
    """First upward crossing of N_detect (np.inf if none within the horizon); uses the verified Phase-3 trajectory runner."""
    p, x0 = theta_to_model(theta, mu, h, K)
    if x0[0] + x0[1] >= N_detect:
        raise ValueError("initial burden must be below the detection threshold")
    t = run_trajectory(p, x0, [N_detect], horizon, method=method, rtol=rtol, atol=atol, main_threshold=N_detect)
    k = f"{N_detect:g}".replace(".", "p")
    T = t.get(f"T_{k}", float("nan"))
    return float(T) if np.isfinite(T) else float("inf")


def observation_times(a: ObservationAssumptions, T_detect: float) -> np.ndarray:
    """Visit times; stops at the detection time when ``truncate_at_detection`` (the visit at/after detection is not an untreated observation)."""
    end = a.follow_up_months
    if a.truncate_at_detection and np.isfinite(T_detect):
        end = min(end, T_detect)
    ts = np.arange(a.schedule_months, end + 1e-9, a.schedule_months)
    if a.include_baseline:
        ts = np.concatenate([[0.0], ts])
    return ts


@dataclass
class Dataset:
    theta_true: dict
    assumptions: ObservationAssumptions
    times: np.ndarray
    y: dict                               # channel -> observed signal (NaN where censored)
    censored: dict                        # channel -> bool array: True = '< LOD'
    T_detect: float                       # true first-crossing time of N_detect (inf if none) -- a hidden truth
    seed: int = 0
    meta: dict = field(default_factory=dict)


def generate_dataset(theta: dict, a: ObservationAssumptions, seed: int, mu: float = 0.05, h: float = 0.02, K: float = 1.0,
                     solver=None) -> Dataset:
    """Synthetic data from known ground truth: hidden states -> signals with noise and LOD -> truncated at detection."""
    solver = solver or {"method": "Radau", "rtol": 1e-12, "atol": 1e-15}
    T_det = detection_time(theta, a.detection_threshold, max(a.follow_up_months, 1.0) * 4, mu=mu, h=h, K=K, **solver)
    times = observation_times(a, T_det)
    hid = simulate_hidden(theta, times, mu=mu, h=h, K=K, **solver)
    rng = np.random.default_rng(seed)
    y, cens = {}, {}
    for ch in CHANNELS[a.observed]:
        true = {"N": hid.N, "P": hid.X[:, 0], "D": hid.X[:, 1], "e": hid.X[:, 2]}[ch]
        scale = a.immune_scale_ce if ch == "e" else a.signal_scale_c
        sig = scale * true * np.exp(a.noise_sigma_ln * rng.standard_normal(len(times)))
        c = sig < a.lod
        y[ch] = np.where(c, np.nan, sig)
        cens[ch] = c
    return Dataset(dict(theta), a, times, y, cens, T_det, seed)


# ----------------------------------------------------------------------------------------------- event observation
def visit_schedule(step: float, horizon: float, jitter_sd: float = 0.0, rng=None) -> np.ndarray:
    v = np.arange(step, horizon + 1e-9, step)
    if jitter_sd > 0 and rng is not None:
        v = np.sort(np.clip(v + jitter_sd * rng.standard_normal(len(v)), 1e-6, horizon))
    return v


def interval_censor(T: float, visits: np.ndarray, horizon: float):
    """Event observed through visits: returns (L, R, event) with L < T <= R (R = inf, event = 0 if no crossing by the last visit).

    A relapse at time T is first noticed at the first visit v_j >= T (L = previous visit or 0). If T exceeds the follow-up or
    never occurs, the record is right-censored at the last visit before the horizon.
    """
    visits = np.asarray(visits, float)
    visits = visits[visits <= horizon]
    if np.isfinite(T) and T <= (visits[-1] if len(visits) else 0.0):
        j = int(np.searchsorted(visits, T, side="left"))
        L = visits[j - 1] if j > 0 else 0.0
        return float(L), float(visits[j]), 1
    return float(visits[-1] if len(visits) else 0.0), float("inf"), 0


def exact_censor(T: float, horizon: float):
    """Exact event time or right-censoring at the horizon: (time, event)."""
    return (float(T), 1) if (np.isfinite(T) and T <= horizon) else (float(horizon), 0)
