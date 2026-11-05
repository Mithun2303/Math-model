"""ODE integration with first-upward-crossing relapse detection."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .model import rhs
from .parameters import Params


@dataclass
class Simulation:
    t: np.ndarray
    y: np.ndarray  # shape (3, len(t)): P, D, E
    event: int     # 1 if N reached N_detect within horizon, else 0 (right-censored)
    time: float    # relapse time if event == 1, else the horizon

    @property
    def P(self): return self.y[0]
    @property
    def D(self): return self.y[1]
    @property
    def E(self): return self.y[2]
    @property
    def N(self): return self.y[0] + self.y[1]


def simulate(p: Params, y0, N_detect: float = 0.1, horizon: float = 240.0,
             method: str = "DOP853", rtol: float = 1e-9, atol: float = 1e-12,
             n_points: int = 2000) -> Simulation:
    """Integrate to the horizon, stopping at the first upward crossing of N = N_detect.

    If the initial burden is already >= N_detect the relapse time is 0.
    A trajectory that never crosses is censored at the horizon (no invented time).
    """
    y0 = np.asarray(y0, dtype=float)
    if y0[0] + y0[1] >= N_detect:
        return Simulation(np.array([0.0]), y0[:, None], 1, 0.0)

    def crossing(t, y, p):
        return y[0] + y[1] - N_detect
    crossing.terminal = True
    crossing.direction = 1

    sol = solve_ivp(rhs, (0.0, horizon), y0, args=(p,), method=method, rtol=rtol, atol=atol,
                    events=crossing, t_eval=np.linspace(0.0, horizon, n_points))
    if not sol.success:
        raise RuntimeError(sol.message)
    if sol.t_events[0].size:
        t_ev = float(sol.t_events[0][0])
        y_ev = sol.y_events[0][0]
        keep = sol.t < t_ev
        t = np.append(sol.t[keep], t_ev)
        y = np.column_stack([sol.y[:, keep], y_ev])
        return Simulation(t, y, 1, t_ev)
    return Simulation(sol.t, sol.y, 0, horizon)
