import pytest
from hx.properties import water, engine_oil, props
import numpy as np


def test_oil_matches_incropera_ex_11_1():
    """Book used oil at 80 C (353 K): cp 2131, mu 3.25e-2, k 0.138."""
    p = engine_oil(353 - 273.15)
    assert abs(p["cp"] - 2131) < 1
    assert abs(p["mu"] - 3.25e-2) < 0.01e-2
    assert abs(p["k"] - 0.138) < 0.001


def test_water_close_to_incropera_table():
    """Book used water at 35 C: cp 4178, mu 725e-6, k 0.625."""
    p = water(35)
    assert abs(p["cp"] - 4178) / 4178 < 0.005
    assert abs(p["mu"] - 725e-6) / 725e-6 < 0.02
    assert abs(p["k"] - 0.625) / 0.625 < 0.01


def test_oil_out_of_range_raises():
    with pytest.raises(ValueError):
        engine_oil(200)


def test_dispatch():
    assert props("water", 35) == water(35)
    with pytest.raises(ValueError):
        props("mercury", 35)

from hx.properties import turbine_oil, jet_fuel


def test_oil_viscosity_hits_datasheet_points():
    """Walther fit must reproduce both Mobil Jet Oil II datasheet values exactly."""
    for T, nu_cSt in [(40, 27.6), (100, 5.1)]:
        p = turbine_oil(T)
        assert np.isclose(p["mu"] / p["rho"] * 1e6, nu_cSt, rtol=1e-6)


def test_oil_viscosity_interpolated():
    p = turbine_oil(80)
    assert np.isclose(p["mu"] / p["rho"] * 1e6, 7.896, rtol=1e-3)
    assert np.isclose(p["mu"], 7.563e-3, rtol=1e-3)


def test_oil_conductivity_matches_nist():
    assert abs(turbine_oil(28.4)["k"] - 0.1451) < 0.0015
    assert abs(turbine_oil(126.8)["k"] - 0.1307) < 0.0015


def test_oil_out_of_range():
    with pytest.raises(ValueError):
        turbine_oil(250)


def test_jet_fuel_surrogate_values():
    p = jet_fuel(80)
    assert np.isclose(p["rho"], 705.4, rtol=2e-3)
    assert np.isclose(p["cp"], 2408, rtol=5e-3)
    assert np.isclose(p["mu"], 6.317e-4, rtol=5e-3)


def test_jet_fuel_rejects_frozen_surrogate():
    with pytest.raises(ValueError):
        jet_fuel(-20)


def test_new_fluids_dispatch():
    assert props("turbine_oil", 80) == turbine_oil(80)
    assert props("jet_fuel", 80) == jet_fuel(80)