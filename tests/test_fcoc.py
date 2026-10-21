import numpy as np
from dataclasses import replace
from hx.plate_fin import OffsetStripFin, PlateFinCore, rate_plate_fin_auto
from hx.fcoc import evaluate_mission, is_feasible

oil_fin = OffsetStripFin(s=1.5e-3, h=6.0e-3, t=0.15e-3, l=3.0e-3)
fuel_fin = OffsetStripFin(s=1.0e-3, h=2.5e-3, t=0.10e-3, l=3.0e-3)
best = PlateFinCore(L_hot=0.4, L_cold=0.2, N_hot=1, fin_hot=oil_fin, fin_cold=fuel_fin)


def test_takeoff_reference_point():
    r = rate_plate_fin_auto(best, "turbine_oil", "jet_fuel", 0.5, 1.10, 130, 40)
    assert np.isclose(r["q"], 45355, rtol=2e-3)
    assert np.isclose(r["Tc_out"], 57.98, atol=0.1)
    assert np.isclose(r["Re_hot"], 265.9, rtol=2e-3)
    assert np.isclose(r["mass"], 0.7167, rtol=1e-3)
    assert r["iterations"] < 20


def test_best_design_is_feasible_across_mission():
    res = evaluate_mission(best)
    assert is_feasible(res)
    assert np.isclose(res["idle_descent"]["result"]["Tc_out"], 109.2, atol=0.3)


def test_takeoff_is_the_sizing_point():
    """Smallest margin (achieved / required) should be at take-off."""
    res = evaluate_mission(best)
    from hx.fcoc import MISSION
    margins = {k: v["result"]["q"] / MISSION[k][2] for k, v in res.items()}
    assert min(margins, key=margins.get) == "takeoff"


def test_too_many_layers_fails_gracefully():
    """3 oil layers drop oil Re below 120: infeasible, but no crash."""
    res = evaluate_mission(replace(best, N_hot=3))
    assert not is_feasible(res)


def test_titanium_heavier_and_worse():
    ti = replace(best, k_metal=6.7, rho_metal=4430)
    r_al = rate_plate_fin_auto(best, "turbine_oil", "jet_fuel", 0.5, 1.10, 130, 40)
    r_ti = rate_plate_fin_auto(ti, "turbine_oil", "jet_fuel", 0.5, 1.10, 130, 40)
    assert r_ti["mass"] > 1.5 * r_al["mass"]
    assert r_ti["q"] < 0.65 * r_al["q"]