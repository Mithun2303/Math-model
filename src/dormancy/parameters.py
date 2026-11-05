"""Parameter definitions, units, validation and config loading."""
from __future__ import annotations

from dataclasses import dataclass, fields, replace
from pathlib import Path

import yaml


@dataclass(frozen=True)
class Params:
    """Model parameters. Time in months; burdens normalised to K.

    Defaults are illustrative mathematical choices, not fitted estimates.
    """
    K: float = 1.0
    r: float = 0.12
    q: float = 0.04
    a: float = 0.005
    delta_D: float = 0.001
    k: float = 0.1
    s: float = 0.02
    mu: float = 0.05
    rho: float = 0.04
    h: float = 0.02
    eta: float = 0.08

    def __post_init__(self):
        for f in fields(self):
            if getattr(self, f.name) < 0:
                raise ValueError(f"parameter {f.name} must be nonnegative")
        for name in ("K", "mu", "h"):
            if getattr(self, name) <= 0:
                raise ValueError(f"parameter {name} must be positive")

    def with_(self, **changes) -> "Params":
        return replace(self, **changes)

    @property
    def E_free(self) -> float:
        """Immune level at the tumour-free equilibrium, s/mu."""
        return self.s / self.mu

    @property
    def beta(self) -> float:
        """Gauge-invariant effective immune killing, beta = k s / mu (see normalized.py)."""
        return self.k * self.s / self.mu

    @property
    def r_eff(self) -> float:
        """Net low-burden growth rate r - q delta_D/(a+delta_D) (tumour-free boundary: beta = r_eff)."""
        d = self.a + self.delta_D
        return self.r - (self.q * self.delta_D / d if d > 0 else 0.0)

    @property
    def phi(self) -> float:
        """N/P on the dormancy manifold at equilibrium, 1 + q/(a+delta_D)."""
        return 1.0 + self.q / (self.a + self.delta_D)


def load_config(path: str | Path) -> dict:
    """Load a YAML config; returns dict with a Params under key 'params'."""
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    cfg["params"] = Params(**cfg["parameters"])
    return cfg
