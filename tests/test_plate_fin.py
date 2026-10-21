import numpy as np
import pytest
from hx.plate_fin import (OffsetStripFin, j_manglik_bergles, f_manglik_bergles,
                          fin_efficiency, surface_efficiency)
from hx.double_pipe import Stream
from hx.plate_fin import (PlateFinCore, side_geometry, side_performance,
                          core_mass, core_volume, rate_plate_fin)
from hx.methods import effectiveness_crossflow_unmixed

# Typical aero-engine offset strip fin
fin = OffsetStripFin(s=1.5e-3, h=6.0e-3, t=0.15e-3, l=3.0e-3)


def test_ratios():
    assert np.isclose(fin.alpha, 0.25)
    assert np.isclose(fin.delta, 0.05)
    assert np.isclose(fin.gamma, 0.1)


def test_hydraulic_diameter():
    assert np.isclose(fin.D_h, 2.2967e-3, rtol=1e-4)


def test_j_and_f_reference_values():
    """Hand-evaluated Manglik-Bergles values for this geometry."""
    assert np.isclose(j_manglik_bergles(1000, fin), 0.01529, rtol=1e-3)
    assert np.isclose(f_manglik_bergles(1000, fin), 0.06276, rtol=1e-3)


def test_j_and_f_fall_with_reynolds():
    for corr in (j_manglik_bergles, f_manglik_bergles):
        assert corr(300, fin) > corr(1000, fin) > corr(3000, fin)


def test_j_over_f_is_physical():
    """Offset strip fins typically give j/f ~ 0.2-0.3 (Reynolds analogy territory)."""
    for Re in (300, 1000, 3000):
        assert 0.15 < j_manglik_bergles(Re, fin) / f_manglik_bergles(Re, fin) < 0.35


def test_outside_fitted_range_raises():
    with pytest.raises(ValueError):
        j_manglik_bergles(50, fin)
    with pytest.raises(ValueError):
        f_manglik_bergles(20000, fin)


def test_fin_efficiency_aluminium():
    """h = 500 W/m2K, aluminium k = 170, t = 0.15 mm, L = h_fin/2 = 3 mm -> mL = 0.594."""
    eta = fin_efficiency(500, 170, 0.15e-3, 3.0e-3)
    assert np.isclose(eta, 0.8969, rtol=1e-3)


def test_fin_efficiency_limits():
    assert fin_efficiency(1e-6, 170, 0.15e-3, 3e-3) > 0.999     # tiny h -> ~100%
    assert fin_efficiency(500, 170, 0.15e-3, 3e-3) > fin_efficiency(500, 20, 0.15e-3, 3e-3)  # Al beats Ti-ish k


def test_surface_efficiency():
    assert np.isclose(surface_efficiency(1.0, 0.8), 1.0)
    assert np.isclose(surface_efficiency(0.9, 0.8), 0.92)

oil_like = Stream(m_dot=0.4, cp=2100, mu=2e-3, k=0.13, T_in=120, rho=850)
fuel_like = Stream(m_dot=0.6, cp=2000, mu=1.2e-3, k=0.11, T_in=40, rho=800)
core = PlateFinCore(L_hot=0.2, L_cold=0.15, N_hot=3, fin_hot=fin, fin_cold=fin)


def test_side_geometry():
    g = side_geometry(core, "hot")
    assert g["n_ch"] == 90
    assert np.isclose(g["A_ff"], 2.43e-3)
    assert np.isclose(g["A"], 0.84645, rtol=1e-4)
    assert side_geometry(core, "cold")["n_ch"] == 121


def test_side_performance_hot():
    p = side_performance(core, "hot", oil_like)
    assert np.isclose(p["Re"], 189.0, rtol=1e-3)
    assert np.isclose(p["h"], 1218.9, rtol=1e-3)
    assert np.isclose(p["eta_o"], 0.8291, rtol=1e-3)
    assert np.isclose(p["dP"], 1044.0, rtol=1e-3)


def test_mass_and_volume_hand_calc():
    assert np.isclose(core_mass(core), 0.7106, rtol=1e-3)
    assert np.isclose(core_volume(core), 1.4115e-3, rtol=1e-4)


def test_rating_reference_case():
    r = rate_plate_fin(core, oil_like, fuel_like)
    assert np.isclose(r["UA"], 465.24, rtol=1e-3)
    assert np.isclose(r["q"], 24255, rtol=1e-3)
    assert np.isclose(r["Th_out"], 91.12, atol=0.02)
    assert np.isclose(r["kW_per_kg"], 34.13, rtol=1e-3)


def test_rating_energy_balance_and_consistency():
    r = rate_plate_fin(core, oil_like, fuel_like)
    q_hot = oil_like.m_dot * oil_like.cp * (oil_like.T_in - r["Th_out"])
    q_cold = fuel_like.m_dot * fuel_like.cp * (r["Tc_out"] - fuel_like.T_in)
    assert np.isclose(q_hot, q_cold)
    assert np.isclose(r["eps"], effectiveness_crossflow_unmixed(r["NTU"], 0.4 * 2100 / (0.6 * 2000)))


def test_too_many_layers_leaves_correlation_range():
    """10 hot layers slows the oil to Re ~ 57: must refuse to extrapolate."""
    big = PlateFinCore(L_hot=0.2, L_cold=0.15, N_hot=10, fin_hot=fin, fin_cold=fin)
    with pytest.raises(ValueError):
        rate_plate_fin(big, oil_like, fuel_like)