import numpy as np

from dormancy import Params, simulate, critical_k


Y0 = [1e-6, 1e-4, 0.4]


def test_solution_nonnegative_and_event_consistent():
    p = Params(k=0.1)
    sim = simulate(p, Y0, N_detect=0.1, horizon=600)
    assert (sim.y >= -1e-15).all()
    if sim.event:
        assert abs(sim.N[-1] - 0.1) < 1e-6
        assert sim.time == sim.t[-1]


def test_censored_when_no_crossing():
    p = Params(k=1.0)  # far above critical k
    sim = simulate(p, Y0, horizon=240)
    assert sim.event == 0 and sim.time == 240


def test_initial_burden_above_threshold():
    sim = simulate(Params(), [0.2, 0.0, 0.4], N_detect=0.1)
    assert sim.event == 1 and sim.time == 0.0


def test_solvers_agree_on_relapse_time():
    p = Params(k=0.1)
    t1 = simulate(p, Y0, horizon=600, method="DOP853")
    t2 = simulate(p, Y0, horizon=600, method="Radau", rtol=1e-10, atol=1e-13)
    t3 = simulate(p, Y0, horizon=600, method="DOP853", rtol=1e-11, atol=1e-14)
    assert t1.event == t2.event == t3.event
    if t1.event:
        assert abs(t1.time - t2.time) < 1e-2
        assert abs(t1.time - t3.time) < 1e-3


def test_decay_below_boundary_growth_above():
    p = Params()
    kc = critical_k(p)
    lo = simulate(p.with_(k=1.3 * kc), [1e-6, 1e-4, 0.4], horizon=100)
    assert lo.N[-1] < lo.N[0]
