import pytest
from hx.properties import water, engine_oil, props


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