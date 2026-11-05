"""Phase 4: structural identifiability (exact) and numerical identifiability (forward sensitivities).

**Proof standard.** For the unknown set U (parameters and/or initial states) and an observed output y(t) = h(x(t), theta),
the *identifiability matrix* has rows grad_U y_k (k = 0, 1, ...), where y_k are the Taylor coefficients of y(t) at t = 0.
They are computed EXACTLY with forward-mode automatic differentiation through truncated power series, in a prime field
GF(p) (p = 2^61 - 1) or in the rationals (Fraction), at a random point with strictly positive states and generic parameters.
  * full column rank over GF(p) at one point  =>  full rank over Q at that point  =>  the generating-series map is locally
    injective  =>  generic LOCAL structural identifiability of U (finite-to-one). It does NOT prove global identifiability.
  * rank deficiency is only EVIDENCE of non-identifiability (a random point could be non-generic). It is certified here by
    taking the exact rational null vector and checking the corresponding symmetry analytically (see report section 3).
  * a numerically full-rank sensitivity matrix is NOT a structural proof (see `numerical_spectrum`).

Time series of y_k are in the observation window t -> 0 only: this is an analytic-germ statement about the model, not a
statement about what a noisy, sparse, truncated data set can resolve (that is practical identifiability, synthetic_validation.py).
"""
from __future__ import annotations

import random
from dataclasses import dataclass
from fractions import Fraction

import numpy as np

P_MOD = (1 << 61) - 1


# ----------------------------------------------------------------------------------------------- fields
class GF:
    """Element of GF(p), p = 2^61 - 1 (enough operations for +, -, *, / used below)."""
    __slots__ = ("v",)

    def __init__(self, v):
        self.v = int(v) % P_MOD

    def __add__(self, o): return GF(self.v + _v(o))
    __radd__ = __add__
    def __sub__(self, o): return GF(self.v - _v(o))
    def __rsub__(self, o): return GF(_v(o) - self.v)
    def __mul__(self, o): return GF(self.v * _v(o))
    __rmul__ = __mul__
    def __neg__(self): return GF(-self.v)
    def inv(self):
        if self.v == 0:
            raise ZeroDivisionError("division by zero in GF(p)")
        return GF(pow(self.v, P_MOD - 2, P_MOD))
    def __truediv__(self, o): return self * GF(_v(o)).inv()
    def __rtruediv__(self, o): return GF(_v(o)) * self.inv()
    def __eq__(self, o): return self.v == _v(o)
    def __hash__(self): return hash(self.v)
    def __repr__(self): return f"GF({self.v})"


def _v(o):
    return o.v if isinstance(o, GF) else int(o) % P_MOD


def to_gf(fr: Fraction) -> GF:
    return GF(fr.numerator) / GF(fr.denominator)


# ----------------------------------------------------------------------------------------------- dual numbers
class Dual:
    """value + gradient (list or None for a gradient-free constant) over a field."""
    __slots__ = ("v", "g")

    def __init__(self, v, g=None):
        self.v, self.g = v, g

    def __add__(self, o):
        if not isinstance(o, Dual):
            return Dual(self.v + o, self.g)
        return Dual(self.v + o.v, _gadd(self.g, o.g, 1, 1))

    __radd__ = __add__

    def __sub__(self, o):
        if not isinstance(o, Dual):
            return Dual(self.v - o, self.g)
        return Dual(self.v - o.v, _gadd(self.g, o.g, 1, -1))

    def __rsub__(self, o):
        return Dual(o - self.v, None if self.g is None else [-x for x in self.g])

    def __neg__(self):
        return Dual(-self.v, None if self.g is None else [-x for x in self.g])

    def __mul__(self, o):
        if not isinstance(o, Dual):
            return Dual(self.v * o, None if self.g is None else [x * o for x in self.g])
        g = None
        if self.g is not None or o.g is not None:
            ga = [] if self.g is None else [o.v * x for x in self.g]
            gb = [] if o.g is None else [self.v * x for x in o.g]
            n = len(self.g if self.g is not None else o.g)
            g = [(ga[i] if ga else 0) + (gb[i] if gb else 0) for i in range(n)]
        return Dual(self.v * o.v, g)

    __rmul__ = __mul__

    def inv(self):
        iv = 1 / self.v if not isinstance(self.v, GF) else self.v.inv()
        g = None if self.g is None else [-(x * iv * iv) for x in self.g]
        return Dual(iv, g)

    def __truediv__(self, o):
        if not isinstance(o, Dual):
            return self * (1 / o if not isinstance(o, GF) else o.inv())
        return self * o.inv()


def _gadd(a, b, ca, cb):
    if a is None and b is None:
        return None
    if a is None:
        return [cb * x for x in b]
    if b is None:
        return [ca * x for x in a]
    return [ca * x + cb * y for x, y in zip(a, b)]


# ----------------------------------------------------------------------------------------------- series model
PARAMS = ["r", "q", "a", "delta_D", "beta", "rho", "eta", "h", "mu", "K"]
ICS = ["P0", "D0", "e0"]
SCALES = ["c", "ce"]            # unknown signal scales: y_N = c * N, y_e = ce * e
OBS_SCENARIOS = {"N": ["N"], "PD": ["P", "D"], "Ne": ["N", "e"], "PDe": ["P", "D", "e"]}


@dataclass
class PDEModel:
    """Normalised P-D-e model as a power-series ODE (Taylor-coefficient recurrence).

    dP = r P (1 - N/K) - q P + a D - beta e P ;  dD = q P - (a + delta_D) D ;  de = mu (1 - e) + rho g e - eta N e,  g = P/(h+P).
    Flags allow the linear / frozen-immune sub-models used in the reduction analysis.
    """
    logistic: bool = True
    immune: bool = True

    def coefficients(self, order, par, x0, one, zero):
        """Return lists P, D, e of series coefficients (Dual) up to ``order``."""
        r, q, a, dd, beta, rho, eta, h, mu, K = (par[k] for k in PARAMS)
        invK = K.inv()
        P, D, e = [x0[0]], [x0[1]], [x0[2]]
        PN, eP, eN, G, gE, HP = [], [], [], [], [], []

        def conv(A, B, k):
            s = A[0] * B[k]
            for i in range(1, k + 1):
                s = s + A[i] * B[k - i]
            return s

        for k in range(order):
            N_k = P[k] + D[k]
            Nl = [P[i] + D[i] for i in range(k + 1)]
            PN.append(conv(P, Nl, k))
            eP.append(conv(e, P, k))
            eN.append(conv(e, Nl, k))
            # g = P / (h + P)
            HP.append((h + P[k]) if k == 0 else P[k])
            if k == 0:
                G.append(P[0] * HP[0].inv())
            else:
                s = P[k]
                for i in range(1, k + 1):
                    s = s - HP[i] * G[k - i]
                G.append(s * HP[0].inv())
            gE.append(conv(G, e, k))
            FP = r * P[k] - q * P[k] + a * D[k] - beta * eP[k]
            if self.logistic:
                FP = FP - r * invK * PN[k]
            FD = q * P[k] - (a + dd) * D[k]
            if self.immune:
                Fe = (mu * one if k == 0 else mu * zero) - mu * e[k] + rho * gE[k] - eta * eN[k]
            else:
                Fe = zero * mu
            sc = GF(k + 1).inv() if isinstance(one.v, GF) else Fraction(1, k + 1)
            P.append(FP * sc)
            D.append(FD * sc)
            e.append(Fe * sc)
        return P, D, e


def _outputs(P, D, e, scenario, c, ce):
    out = []
    for name in OBS_SCENARIOS[scenario]:
        if name == "N":
            out.append([(P[k] + D[k]) * c for k in range(len(P))])
        elif name == "P":
            out.append([x * c for x in P])
        elif name == "D":
            out.append([x * c for x in D])
        elif name == "e":
            out.append([x * ce for x in e])
    return out


def random_point(seed: int, positive: bool = True) -> dict:
    """Random rational point (generic parameters, strictly positive states)."""
    rng = random.Random(seed)
    fr = lambda lo, hi: Fraction(rng.randint(lo, hi), rng.randint(97, 997))
    pt = {"r": fr(20, 90), "q": fr(5, 30), "a": fr(3, 25), "delta_D": fr(1, 15), "beta": fr(5, 40), "rho": fr(5, 40),
          "eta": fr(5, 60), "h": fr(4, 40), "mu": fr(10, 60), "K": Fraction(rng.randint(150, 300), rng.randint(97, 201)),
          "P0": fr(3, 30), "D0": fr(3, 30), "e0": fr(30, 120), "c": fr(40, 160), "ce": fr(40, 160)}
    return pt


def identifiability_matrix(model: PDEModel, unknown: list[str], point: dict, order: int, scenario: str = "N",
                           field: str = "gf"):
    """Rows: gradient of every output Taylor coefficient w.r.t. the unknowns. Columns follow ``unknown``."""
    conv = to_gf if field == "gf" else (lambda x: Fraction(x))
    n = len(unknown)

    def make(name):
        val = conv(point[name])
        g = [conv(1) if unknown[i] == name else conv(0) for i in range(n)] if name in unknown else None
        return Dual(val, g)

    zero_v, one_v = conv(0), conv(1)
    par = {k: make(k) for k in PARAMS}
    x0 = [make("P0"), make("D0"), make("e0")]
    c, ce = make("c"), make("ce")
    one, zero = Dual(one_v, None), Dual(zero_v, None)
    P, D, e = model.coefficients(order, par, x0, one, zero)
    rows = []
    for series in _outputs(P, D, e, scenario, c, ce):
        for k in range(order + 1):
            d = series[k]
            rows.append(list(d.g) if d.g is not None else [conv(0)] * n)
    return rows


def rank_gf(rows) -> int:
    """Rank over GF(p) (Gaussian elimination)."""
    M = [[x.v if isinstance(x, GF) else int(x) % P_MOD for x in r] for r in rows]
    rank, ncols = 0, len(M[0]) if M else 0
    for col in range(ncols):
        piv = None
        for i in range(rank, len(M)):
            if M[i][col] % P_MOD:
                piv = i
                break
        if piv is None:
            continue
        M[rank], M[piv] = M[piv], M[rank]
        inv = pow(M[rank][col], P_MOD - 2, P_MOD)
        M[rank] = [(x * inv) % P_MOD for x in M[rank]]
        for i in range(len(M)):
            if i != rank and M[i][col]:
                f = M[i][col]
                M[i] = [(x - f * y) % P_MOD for x, y in zip(M[i], M[rank])]
        rank += 1
    return rank


def nullspace_q(rows):
    """Exact rational null space (list of Fraction vectors) and rank over Q."""
    from sympy import Matrix, Rational
    M = Matrix([[Rational(x.numerator, x.denominator) for x in r] for r in rows])
    ns = M.nullspace()
    return [[Fraction(int(x.p), int(x.q)) for x in v] for v in ns], M.rank()


def structural_report(unknown: list[str], scenario: str, model: PDEModel | None = None, max_order: int = 30,
                      n_points: int = 3, seed: int = 20261201, flat_orders: int = 6) -> dict:
    """Generic local structural identifiability (rank mod p) with rank-vs-order saturation, over several random points."""
    model = model or PDEModel()
    ranks, curves = [], []
    for j in range(n_points):
        pt = random_point(seed + 17 * j)
        rows = identifiability_matrix(model, unknown, pt, max_order, scenario, "gf")
        n_out = len(OBS_SCENARIOS[scenario])
        curve = {}
        for k in range(1, max_order + 1):
            sub = []
            for o in range(n_out):
                sub += rows[o * (max_order + 1): o * (max_order + 1) + k + 1]
            curve[k] = rank_gf(sub)
        curves.append(curve)
        ranks.append(curve[max_order])
    flat = all(len(set(list(c.values())[-flat_orders:])) == 1 for c in curves)
    return {"unknown": list(unknown), "scenario": scenario, "n_unknown": len(unknown), "rank": max(ranks), "ranks_per_point": ranks,
            "deficiency": len(unknown) - max(ranks), "rank_vs_order": curves[0], "flat_last_orders": flat,
            "structurally_locally_identifiable": max(ranks) == len(unknown) and flat}


def null_directions(unknown: list[str], scenario: str, model: PDEModel | None = None, order: int = 16, seed: int = 20261201):
    """Exact rational null vectors (non-identifiable directions in the tangent space of the unknowns) at one random point."""
    model = model or PDEModel()
    pt = random_point(seed)
    rows = identifiability_matrix(model, unknown, pt, order, scenario, "q")
    vecs, rk = nullspace_q(rows)
    return pt, vecs, rk


# ----------------------------------------------------------------------------------------------- numerical rank (supplement)
def _aug_rhs_factory(th: np.ndarray, idx: list, n_model: int, n_extra: int):
    """Right-hand side of the augmented (state + forward sensitivity) system; analytic Jacobians (hand-derived)."""
    r, q, a, dd, beta, rho, eta, h, mu, K = th
    n = n_model + n_extra
    idx = np.asarray(idx, int)

    def rhs(t, y):
        P, D, e = y[0], y[1], y[2]
        N = P + D
        hp = h + P
        g = P / hp
        f1 = r * P * (1 - N / K) - q * P + a * D - beta * e * P
        f2 = q * P - (a + dd) * D
        f3 = mu * (1 - e) + rho * g * e - eta * N * e
        out = np.empty(3 + 3 * n)
        out[0], out[1], out[2] = f1, f2, f3
        if n:
            S = y[3:].reshape(3, n)
            A = np.array([[r * (1 - N / K) - r * P / K - q - beta * e, -r * P / K + a, -beta * P],
                          [q, -(a + dd), 0.0],
                          [rho * e * h / hp ** 2 - eta * e, -eta * e, -mu + rho * g - eta * N]])
            B = np.zeros((3, n))
            cols = {0: (P * (1 - N / K), 0.0, 0.0), 1: (-P, P, 0.0), 2: (D, -D, 0.0), 3: (0.0, -D, 0.0), 4: (-e * P, 0.0, 0.0),
                    5: (0.0, 0.0, g * e), 6: (0.0, 0.0, -N * e), 7: (0.0, 0.0, -rho * P * e / hp ** 2), 8: (0.0, 0.0, 1 - e),
                    9: (r * P * N / K ** 2, 0.0, 0.0)}
            for j in range(n_model):
                B[:, j] = cols[int(idx[j])]
            out[3:] = (A @ S + B).ravel()
        return out

    return rhs


def forward_sensitivities(theta: dict, x0, sens_params: list, times, rtol: float = 1e-10, atol: float = 1e-13,
                          ic_sens: dict | None = None, method: str = "DOP853"):
    """States and forward sensitivities of x(t_j) with respect to ``sens_params`` (natural scale).

    ``sens_params`` entries are model-parameter names (in PARAMS) or names listed in ``ic_sens`` (3-vector dx0/dp for
    initial-state quantities such as ln N0). Returns (X (T,3), S (T,3,n)). Analytic Jacobians (see _aug_rhs_factory).
    """
    from scipy.integrate import solve_ivp
    ic_sens = ic_sens or {}
    model_names = [p for p in sens_params if p in PARAMS]
    extra = [p for p in sens_params if p not in PARAMS]
    order = model_names + extra                      # columns: model params first, then initial-state quantities
    idx = [PARAMS.index(p) for p in model_names]
    th = np.array([theta[k] for k in PARAMS], float)
    rhs = _aug_rhs_factory(th, idx, len(model_names), len(extra))
    n = len(order)
    S0 = np.zeros((3, n))
    for j, nm in enumerate(extra):
        S0[:, len(model_names) + j] = ic_sens[nm]
    y0 = np.concatenate([np.asarray(x0, float), S0.ravel()])
    times = np.asarray(times, float)
    tmax = float(times.max()) if len(times) else 1.0
    sol = solve_ivp(rhs, (0.0, tmax), y0, method=method, rtol=rtol, atol=atol, t_eval=times if len(times) else None)
    X = sol.y[:3].T
    S = sol.y[3:].T.reshape(-1, 3, n) if n else np.zeros((len(sol.t), 3, 0))
    # return columns in the order requested by the caller
    perm = [order.index(p) for p in sens_params]
    return X, S[:, :, perm]


def numerical_spectrum(J: np.ndarray, sigma: np.ndarray | float = 1.0):
    """Singular values of the noise-scaled sensitivity matrix J (rows = observations, cols = log-parameters).

    Returns the full spectrum, the right singular vectors (parameter combinations; rows of Vt) and the 'effective' numbers of
    estimable directions at several relative thresholds. The spectrum, not a single threshold, is the result.
    """
    Jw = J / (np.asarray(sigma)[:, None] if np.ndim(sigma) else sigma)
    U, s, Vt = np.linalg.svd(Jw, full_matrices=False)
    return {"singular_values": s, "Vt": Vt,
            "n_above_1": int(np.sum(s > 1.0)), "n_above_3": int(np.sum(s > 3.0)),
            "condition": float(s[0] / s[-1]) if s[-1] > 0 else float("inf")}


# ----------------------------------------------------------------------------------------------- analytic pieces
def linear_regime_symbolic():
    """Linear (low-burden) sub-model, N-only observation: symbolic generating-series rank and invariants.

    x' = M x, M = [[m, a], [q, -s]] with m = r - q - beta e0 (e = e0 frozen at low burden), s = a + delta_D; y = N = P + D.
    The output series has only 4 independent numbers (lambda_+, lambda_-, c_+, c_-) whereas M and x0 carry 6 unknowns
    (m, a, q, s, P0, D0): at least a 2-parameter family of (m, a, q, s, P0, D0) gives the same N(t). Returns the symbolic
    Jacobian rank and an exact null basis.
    """
    import sympy as sp
    m, a, q, s, P0, D0 = sp.symbols("m a q s P0 D0")
    M = sp.Matrix([[m, a], [q, -s]])
    x0 = sp.Matrix([P0, D0])
    one = sp.Matrix([[1, 1]])
    y = [sp.simplify((one * (M ** k) * x0)[0]) for k in range(8)]
    J = sp.Matrix([[sp.diff(yk, v) for v in (m, a, q, s, P0, D0)] for yk in y])
    pt = {m: sp.Rational(3, 29), a: sp.Rational(5, 41), q: sp.Rational(7, 53), s: sp.Rational(11, 37), P0: sp.Rational(2, 17), D0: sp.Rational(9, 23)}
    Jn = J.subs(pt)
    return {"rank": int(Jn.rank()), "n_unknown": 6, "null_dim": 6 - int(Jn.rank()), "series": y, "J_symbolic": J, "point": pt,
            "null_basis": Jn.nullspace()}


# ----------------------------------------------------------------------------------------------- event-only identifiability
def event_time_gradient(theta: dict, N_detect: float = 0.1, horizon: float = 5000.0, keys=None, rel_step: float = 1e-4):
    """Gradient of the first-crossing time T with respect to ln(parameters) (central differences of the event time).

    ``theta`` is a normalised-parameter dict (r, q, a, delta_D, beta, rho, eta, N0, fD, e0, N_detect optional). Only valid where
    the event exists and crosses transversally; returns (T, grad dict) or (None, None).
    """
    from .normalized import NormalizedParams
    from .solver import simulate

    keys = keys or ["r", "q", "a", "delta_D", "beta", "rho", "eta", "N0", "fD", "N_detect"]

    def T_of(th, Nd):
        n = NormalizedParams(r=th["r"], q=th["q"], a=th["a"], delta_D=th["delta_D"], beta=th["beta"], rho=th["rho"], eta=th["eta"])
        p = n.to_params()
        y0 = [th["N0"] * (1 - th["fD"]), th["N0"] * th["fD"], th.get("e0", 1.0) * p.E_free]
        s = simulate(p, y0, Nd, horizon, rtol=1e-11, atol=1e-14, n_points=10)
        return s.time if s.event else None

    T0 = T_of(theta, N_detect)
    if T0 is None:
        return None, None
    grad = {}
    for k in keys:
        th_p, th_m = dict(theta), dict(theta)
        Nd_p = Nd_m = N_detect
        if k == "N_detect":
            Nd_p, Nd_m = N_detect * (1 + rel_step), N_detect * (1 - rel_step)
        elif k == "fD":
            f = theta["fD"]
            th_p["fD"], th_m["fD"] = min(f * (1 + rel_step), 1 - 1e-9), f * (1 - rel_step)
        else:
            th_p[k], th_m[k] = theta[k] * (1 + rel_step), theta[k] * (1 - rel_step)
        Tp, Tm = T_of(th_p, Nd_p), T_of(th_m, Nd_m)
        if Tp is None or Tm is None:
            return None, None
        grad[k] = (Tp - Tm) / (2 * np.log1p(rel_step) if k != "fD" else (np.log(th_p["fD"]) - np.log(th_m["fD"])))
    return T0, grad
