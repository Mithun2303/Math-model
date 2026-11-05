import numpy as np
import pytest

from dormancy import Params, critical_k, simulate
from dormancy.relapse import (T_manifold, T_matched, T_original, T_projected, T_two_mode, k_detection_limit,
                              unstable_mode)
from dormancy.model import jacobian
from dormancy.equilibria import tumour_free

P0 = Params()
KC = critical_k(P0)
Y0 = [1e-6, 1e-4, 0.4]


def num(p, Nd, y0=Y0):
    s = simulate(p, y0, Nd, 1e5, n_points=300)
    return s.time if s.event else np.inf


def test_original_formula_is_ln_ratio_over_lambda():
    p = P0.with_(k=0.5 * KC)
    lam = np.max(np.linalg.eigvals(jacobian(tumour_free(p).state, p)).real)
    assert T_original(p, Y0, 0.1) == pytest.approx(np.log(0.1 / 1.01e-4) / lam, rel=1e-12)
    assert 237 < T_original(p, Y0, 0.1) < 238         # Phase-1 value, about 240 months


def test_unstable_mode_eigenpair_and_normalisation():
    p = P0.with_(k=0.5 * KC)
    m = unstable_mode(p)
    J = jacobian(tumour_free(p).state, p)
    assert np.allclose(J @ m.w, m.lam * m.w, atol=1e-13)
    assert np.allclose(m.l @ J, m.lam * m.l, atol=1e-13)
    assert m.l @ m.w == pytest.approx(1.0)
    assert not unstable_mode(P0.with_(k=1.2 * KC)).valid      # tumour-free stable: no growing mode


def test_projected_estimate_equals_two_mode_for_small_threshold():
    p = P0.with_(k=0.5 * KC)
    assert T_projected(p, Y0, 1e-5) == pytest.approx(T_two_mode(p, Y0, 1e-5), rel=1e-3)


def test_projection_correction_moves_estimate_towards_numeric():
    p = P0.with_(k=0.5 * KC)
    t_num = num(p, 0.1)
    assert abs(T_projected(p, Y0, 0.1) - t_num) < abs(T_original(p, Y0, 0.1) - t_num)


def test_linear_estimate_is_accurate_at_low_threshold():
    p = P0.with_(k=0.5 * KC)
    assert T_projected(p, Y0, 1e-4 * 1.0 + 1.01e-4 * 5) == pytest.approx(num(p, 6.05e-4), rel=0.03)


def test_manifold_estimate_exact_without_immune_feedback():
    # with rho = eta = 0 the immune level is constant and only the logistic nonlinearity remains
    for f in (0.25, 0.5, 0.75):
        p = P0.with_(k=f * KC, rho=0.0, eta=0.0)
        assert T_manifold(p, Y0, 0.1) == pytest.approx(num(p, 0.1), rel=0.02)


def test_manifold_estimate_good_when_Pdet_much_smaller_than_h():
    p = P0.with_(k=0.5 * KC)
    assert 0.01 / p.phi < 0.1 * p.h
    assert T_manifold(p, Y0, 0.01) == pytest.approx(num(p, 0.01), rel=0.01)


def test_matched_estimate_accurate_across_k_and_composition():
    for f in (0.1, 0.5, 0.75):
        for fD in (0.0, 0.5, 0.99):
            p = P0.with_(k=f * KC)
            y0 = [1.01e-4 * (1 - fD), 1.01e-4 * fD, p.E_free]
            assert T_matched(p, y0, 0.1) == pytest.approx(num(p, 0.1, y0), rel=0.01)


def test_detection_limit_divergence_and_no_event_beyond():
    kdet = k_detection_limit(P0, 0.1)
    assert 0 < kdet < KC
    # beyond k_det the stable nonzero burden is below N_det: no event; below it an event occurs
    assert num(P0.with_(k=1.02 * kdet), 0.1, [1e-6, 1e-4, 0.4]) == np.inf
    t1, t2 = num(P0.with_(k=0.9 * kdet), 0.1), num(P0.with_(k=0.99 * kdet), 0.1)
    assert np.isfinite(t1) and np.isfinite(t2) and t2 > 2 * t1
    assert T_matched(P0.with_(k=1.02 * kdet), Y0, 0.1) == np.inf


def test_initial_burden_above_threshold_gives_zero():
    assert T_original(P0.with_(k=0.1), [0.2, 0, 0.4], 0.1) == 0.0
    assert T_two_mode(P0.with_(k=0.1), [0.2, 0, 0.4], 0.1) == 0.0


def test_tumour_free_stable_predicts_no_relapse():
    p = P0.with_(k=1.2 * KC)
    for fn in (T_original, T_projected, T_manifold, T_matched):
        assert fn(p, Y0, 0.1) == np.inf


def test_numerical_convergence_of_event_time_in_tolerance():
    p = P0.with_(k=0.5 * KC)
    t1 = simulate(p, Y0, 0.1, 1e4, rtol=1e-8, atol=1e-12).time
    t2 = simulate(p, Y0, 0.1, 1e4, rtol=1e-11, atol=1e-15).time
    assert abs(t1 - t2) < 1e-4


def test_matched_estimate_tail_is_composition_independent_and_dominant_at_large_threshold():
    from dormancy.relapse import T_matched_parts
    p = P0.with_(k=0.5 * KC)
    parts = []
    for fD in (0.0, 0.5, 0.99):
        y0 = [1.01e-4 * (1 - fD), 1.01e-4 * fD, p.E_free]
        parts.append(T_matched_parts(p, y0, 0.1))
    taus = [x[2] for x in parts]
    assert max(taus) - min(taus) < 1e-6                              # tau depends on parameters only
    assert all(tau / tot > 0.6 for (tot, _, tau) in parts)           # numerical tail is most of T at N_det = 0.1
    tot, closed, tau = T_matched_parts(p, Y0, 1e-3)
    assert tau == 0.0 and tot == closed                              # below N_match the estimate is purely closed form


def test_no_event_just_above_kdet_at_baseline():
    kd = k_detection_limit(P0, 0.1)
    for f in (1.005, 1.03):
        assert num(P0.with_(k=f * kd), 0.1) == np.inf
