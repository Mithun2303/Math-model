"""Phase 5 (Agent 2) tests: model-vs-METABRIC descriptive comparison utilities (no fitting)."""
import shutil
from pathlib import Path

import numpy as np
import pytest

from dormancy.clinical_comparison import (conditional_late_risk, ecdf_at, evaluate_scenario, inside_ci, km_cuminc_grid, km_reference,
                                          late_increment, load_prereg, load_prior_times, map_times, prereg_precedes, ratio_240_120,
                                          scenario_grid, shifted_cdf, signed_difference, sup_norm, verify_prereg, wilson)

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / "configs" / "phase5.yaml"
HASH = ROOT / "results" / "tables" / "p5_prereg_hash.txt"
TAB = ROOT / "results" / "tables"
rng0 = np.random.default_rng(7)
T = np.where(rng0.random(400) < 0.6, rng0.uniform(40, 240, 400), np.inf)


def test_ecdf_equals_kaplan_meier_without_censoring():
    fin = np.isfinite(T)
    t_km = np.where(fin, T, 240.0)            # administrative censoring at the horizon = event indicator 0
    F_km = km_cuminc_grid(t_km, fin.astype(int), np.arange(0, 241.0, 5.0))
    F_e = ecdf_at(T, np.arange(0, 241.0, 5.0))
    assert np.allclose(F_km, F_e, atol=1e-12)


def test_time_origin_transformation_exact():
    for tau, d in [(0, 0), (6, 12), (24, 6)]:
        Tc = map_times(T, tau, d)
        assert np.all(np.isposinf(Tc[~np.isfinite(T)]))
        assert np.allclose(Tc[np.isfinite(T)], T[np.isfinite(T)] + tau + d)
    h = np.array([60.0, 120.0, 240.0])
    assert np.allclose(ecdf_at(map_times(T, 12, 6), h), shifted_cdf(T, h, 18.0))
    assert np.allclose(ecdf_at(map_times(T, 12, 6), h), ecdf_at(map_times(T, 18, 0), h))     # only tau + d matters
    assert np.all(shifted_cdf(T, np.array([0.0, 10.0]), 50.0) == 0)                          # t - shift < 0 -> 0


def test_no_negative_times_and_nonnegative_cdf():
    Tc = map_times(T, 24, 12, death_first=0.272, seed=3)
    assert np.all(Tc[np.isfinite(Tc)] > 0)
    F = ecdf_at(Tc, np.arange(0, 300.0))
    assert F[0] == 0 and np.all(np.diff(F) >= 0) and F[-1] <= 1


def test_death_first_reproducible_and_seeded():
    a = map_times(T, 0, 0, 0.272, 24, seed=11); b = map_times(T, 0, 0, 0.272, 24, seed=11); c = map_times(T, 0, 0, 0.272, 24, seed=12)
    assert np.array_equal(a, b) and not np.array_equal(a, c)


def test_comparison_arithmetic():
    assert signed_difference(0.30, 0.45) == pytest.approx(-0.15)
    assert inside_ci(0.5, 0.48, 0.53) and not inside_ci(0.47, 0.48, 0.53) and inside_ci(0.48, 0.48, 0.53)
    assert sup_norm([0, 0.1, 0.5], [0, 0.3, 0.4]) == pytest.approx(0.2)
    assert late_increment(0.2, 0.5) == pytest.approx(0.3)
    assert conditional_late_risk(0.2, 0.5) == pytest.approx(0.3 / 0.8)
    assert ratio_240_120(0.2, 0.5) == pytest.approx(2.5) and np.isinf(ratio_240_120(0.0, 0.1))
    lo, hi = wilson(50, 100); assert lo < 0.5 < hi and 0 <= lo and hi <= 1


def test_scenario_grid_size_and_unique():
    cfg = load_prereg(CFG)
    g = scenario_grid(cfg)
    assert len(g) == cfg["grid"]["size"] == 1764
    assert len({(s["prior"], s["threshold"], s["tau"], s["d"], s["death_first"]) for s in g}) == 1764


def _fake_km():
    return {h: {"cum_inc": 0.1 * (i + 1), "ci_lo": 0.1 * (i + 1) - 0.02, "ci_hi": 0.1 * (i + 1) + 0.02, "n_at_risk": 100, "n_events_by_horizon": 10}
            for i, h in enumerate([60, 120, 180, 240])}


def test_scenario_reproducibility_and_shift_monotonicity():
    cfg = load_prereg(CFG)
    sc = {"scenario_id": 5, "prior": "P1_base", "threshold": 0.1, "tau": 6, "d": 6, "death_first": 0.105}
    grid_t = np.arange(0.0, 241.0)
    kmF = np.linspace(0, 0.5, 241)
    a = evaluate_scenario(sc, T, _fake_km(), kmF, cfg, grid_t); b = evaluate_scenario(sc, T, _fake_km(), kmF, cfg, grid_t)
    assert np.array_equal(a.pop("_Fg"), b.pop("_Fg")) and a == b
    # increasing tau (death_first = 0) can only lower the simulated cumulative incidence at every horizon
    F = [evaluate_scenario({**sc, "tau": tau, "death_first": 0.0}, T, _fake_km(), kmF, cfg, grid_t) for tau in (0, 6, 12, 24)]
    for h in (60, 120, 180, 240):
        v = [r[f"F_sim_{h}"] for r in F]; assert all(x >= y - 1e-12 for x, y in zip(v, v[1:]))


@pytest.mark.skipif(not (TAB / "p3_mc_P1_samples.csv").exists(), reason="Phase 3 tables absent")
def test_threshold_monotonicity_real_samples():
    cfg = load_prereg(CFG)
    files = cfg["simulation_source"]["priors"]["P1_base"]["files"][:1]
    prev = None
    for th in cfg["grid"]["thresholds"]:
        Tm = load_prior_times(TAB, files, th)
        assert np.all(Tm[np.isfinite(Tm)] > 0)
        F = ecdf_at(Tm, [60, 120, 240])
        if prev is not None:
            assert np.all(F <= prev + 1e-12)          # higher threshold -> fewer/later crossings (same samples)
        prev = F


def test_prereg_hash_verified_and_detects_change(tmp_path):
    assert len(verify_prereg(CFG, HASH)) == 64
    bad = tmp_path / "phase5.yaml"; shutil.copy(CFG, bad); bad.write_text(bad.read_text() + "\n# tweak\n")
    with pytest.raises(RuntimeError):
        verify_prereg(bad, HASH)


def test_prereg_recorded_before_comparison_tables():
    outs = [TAB / "p5_consistency_grid.csv", TAB / "p5_horizon_comparison.csv", TAB / "p5_consistency_robust_findings.csv"]
    if not all(o.exists() for o in outs):
        pytest.skip("comparison tables not generated yet")
    assert prereg_precedes(HASH, *outs)
    assert "recorded_utc" in HASH.read_text()
