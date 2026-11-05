import numpy as np
import pytest

from dormancy import NormalizedParams, Params, critical_k, jacobian, simulate
from dormancy.normalized import (gauge_transform, jacobian_normalized, rhs_normalized, simulate_normalized,
                                 to_normalized_state, to_original_state)
from dormancy.model import rhs

Y0 = [1e-6, 1e-4, 0.4]


def test_roundtrip_and_beta_definition():
    p = Params(k=0.14)
    n = NormalizedParams.from_params(p)
    assert n.beta == pytest.approx(p.k * p.s / p.mu)
    q = n.to_params(s=p.s)
    assert q.k == pytest.approx(p.k) and q.s == p.s
    assert n.to_params().s == n.mu and n.to_params().k == pytest.approx(n.beta)    # default gauge s = mu
    assert NormalizedParams.from_params(n.to_params(s=0.3)).beta == pytest.approx(n.beta)


def test_beta_c_equals_r_eff_and_matches_phase1_boundary():
    p = Params()
    n = NormalizedParams.from_params(p)
    assert n.beta_c == pytest.approx(p.r_eff)
    assert n.beta_c == pytest.approx(critical_k(p) * p.s / p.mu)


def test_rhs_transformation_is_exact():
    p = Params(k=0.2, s=0.07)
    n = NormalizedParams.from_params(p)
    rng = np.random.default_rng(0)
    for _ in range(20):
        y = np.array([rng.uniform(0, .5), rng.uniform(0, .5), rng.uniform(0, 2)])
        z = to_normalized_state(y, p)
        dy = np.array(rhs(0, y, p))
        dz = np.array(rhs_normalized(0, z, n))
        assert np.allclose(dz[:2], dy[:2], rtol=1e-13, atol=1e-15)
        assert dz[2] == pytest.approx(dy[2] * p.mu / p.s, rel=1e-12)       # de/dt = (mu/s) dE/dt
        # Jacobians conjugate by diag(1,1,mu/s)
        S = np.diag([1, 1, p.mu / p.s])
        assert np.allclose(jacobian_normalized(z, n), S @ jacobian(y, p) @ np.linalg.inv(S), atol=1e-13)


@pytest.mark.parametrize("f", [0.25, 0.5, 0.75])
@pytest.mark.parametrize("s", [0.002, 0.02, 2.0])
def test_trajectories_equivalent_across_gauges(f, s):
    base = Params()
    p = gauge_transform(base.with_(k=f * critical_k(base)), s / base.s)
    n = NormalizedParams.from_params(p)
    y0 = [1e-6, 1e-4, p.E_free]
    so = simulate(p, y0, 0.1, 5000, n_points=500)
    t, z, ev, T = simulate_normalized(n, to_normalized_state(y0, p), 0.1, 5000, n_points=500)
    assert so.event == ev == 1
    assert so.time == pytest.approx(T, abs=1e-5)
    # state at the event: P, D agree; E = (s/mu) e
    assert np.allclose(so.y[:2, -1], z[:2, -1], atol=1e-9)
    assert so.y[2, -1] == pytest.approx(to_original_state(z[:, -1], p)[2], rel=1e-8)


def test_gauge_invariance_of_burden_dynamics_in_original_model():
    p = Params(k=0.14)
    times = []
    for c in (0.1, 1.0, 10.0):
        q = gauge_transform(p, c)
        sim = simulate(q, [1e-6, 1e-4, 0.4 * c], 0.1, 3000)
        times.append(sim.time)
    assert max(times) - min(times) < 1e-5
    assert critical_k(gauge_transform(p, 10.0)) == pytest.approx(critical_k(p) / 10.0)    # k_c scales as 1/c; beta_c fixed


def test_nonequivalent_parameter_change_is_detected():
    # changing beta (not just the gauge) must change the dynamics: the test is not vacuous
    p = Params(k=0.14)
    a = simulate(p, Y0, 0.1, 3000).time
    b = simulate(p.with_(k=0.15), Y0, 0.1, 3000).time
    assert abs(a - b) > 1.0


def test_normalized_validation_and_defaults():
    with pytest.raises(ValueError):
        NormalizedParams(mu=0.0)
    with pytest.raises(ValueError):
        NormalizedParams(beta=-1.0)
    from pathlib import Path
    import yaml
    cfg = yaml.safe_load(open(Path(__file__).parents[1] / "configs" / "normalized.yaml"))
    n = NormalizedParams(**cfg["parameters"])
    assert n.beta == pytest.approx(NormalizedParams.from_params(Params()).beta)      # baseline.yaml <-> normalized.yaml
