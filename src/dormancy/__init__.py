"""Tumour dormancy / relapse ODE model: proliferating (P), dormant (D), immune (E)."""
from .parameters import Params, load_config
from .model import rhs, jacobian
from .solver import simulate, Simulation
from .equilibria import tumour_free, nonzero_equilibria, Equilibrium, enumerate_equilibria
from .stability import eigenvalues, classify, tumour_free_stable, critical_k
from .normalized import NormalizedParams

__all__ = [
    "Params", "load_config", "rhs", "jacobian", "simulate", "Simulation",
    "tumour_free", "nonzero_equilibria", "Equilibrium", "enumerate_equilibria",
    "eigenvalues", "classify", "tumour_free_stable", "critical_k", "NormalizedParams",
]
