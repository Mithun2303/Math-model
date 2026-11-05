"""Phase 4: hypothetical mapping from model-defined recurrence events to a clinical RFS-type endpoint.

Everything here is a NUISANCE ASSUMPTION used to quantify how much the apparent event-time alignment can move, not an estimate
of any clinical delay distribution (there are no data to support one). The model clock starts at the post-treatment state
(t = 0); the clinical clock (METABRIC) has an undocumented origin (diagnosis for OS; RFS undocumented) so

    t_clin = tau + t_model + d          (recurrence recorded at threshold crossing + detection delay d)
    t_clin = tau + t_model + L          (event recorded as disease-specific death without a recorded recurrence, lag L)

tau = unknown diagnosis(-or-surgery)-to-model-origin offset; d = biological-to-clinically-recorded delay; L = relapse-to-death lag.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ClinicalMapping:
    tau: float = 0.0                 # months: clinical origin -> model t = 0
    detect_delay: float = 0.0        # months added at recurrence (recorded later than the model crossing)
    death_first_fraction: float = 0.0    # fraction of events recorded as disease-specific death with no earlier recorded recurrence
    death_lag_mean: float = 24.0     # mean (exponential) lag from model crossing to death when death is the first recorded event


def map_event_times(T_model: np.ndarray, m: ClinicalMapping, rng: np.random.Generator | None = None) -> np.ndarray:
    """Clinical event times for events (finite T_model); non-events (inf/NaN) pass through as inf."""
    T = np.asarray(T_model, float)
    out = np.where(np.isfinite(T), T + m.tau + m.detect_delay, np.inf)
    if m.death_first_fraction > 0:
        rng = rng or np.random.default_rng(0)
        is_death = (rng.random(len(T)) < m.death_first_fraction) & np.isfinite(T)
        lag = rng.exponential(m.death_lag_mean, len(T))
        out = np.where(is_death, T + m.tau + lag, out)
    return out


def cumulative_incidence(T: np.ndarray, times: np.ndarray) -> np.ndarray:
    """Empirical F(t) = P(T <= t) over ``times`` (T may contain inf = no event)."""
    T = np.sort(np.asarray(T, float))
    return np.searchsorted(T, np.asarray(times, float), side="right") / len(T)


def time_to_level(F: np.ndarray, times: np.ndarray, level: float) -> float:
    """First time at which F reaches ``level`` (np.nan if never)."""
    idx = np.where(F >= level)[0]
    return float(times[idx[0]]) if len(idx) else float("nan")


def sup_distance(F1: np.ndarray, F2: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(F1) - np.asarray(F2))))


def dkw_epsilon(n: int, alpha: float = 0.05) -> float:
    """Dvoretzky-Kiefer-Wolfowitz sup-norm half-width: sampling noise of an empirical CDF of size n."""
    return float(np.sqrt(np.log(2.0 / alpha) / (2.0 * n)))
