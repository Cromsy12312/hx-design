"""Fuel-cooled oil cooler (FCOC) case study: CFM56-class turbofan.

Mission data are ILLUSTRATIVE ASSUMPTIONS (manufacturer FCOC loads are not
published). Fuel flows are in the range of ICAO emissions databank LTO values.
"""
from hx.plate_fin import rate_plate_fin_auto

FUEL_LIMIT_C = 150.0    # Jet A thermal stability: ASTM coker test 149 C; NASA CR-72951 cites 177 C
M_OIL = 0.5             # kg/s, assumed
T_OIL_IN = 130.0        # C, oil returning from bearing chambers, assumed

# name: (fuel flow kg/s, fuel inlet C, required oil heat rejection W)
MISSION = {
    "takeoff":      (1.10, 40.0, 45e3),
    "climb":        (0.90, 40.0, 40e3),
    "cruise":       (0.35, 20.0, 30e3),
    "idle_descent": (0.12, 30.0, 18e3),
}


def evaluate_mission(core, mission=MISSION, m_oil=M_OIL, T_oil_in=T_OIL_IN,
                     fuel_limit=FUEL_LIMIT_C):
    """Rate the core at every mission point.

    For each point:
      - try rate_plate_fin_auto('turbine_oil' hot, 'jet_fuel' cold)
      - except ValueError -> the point is outside the j/f correlation range:
            record valid=False, ok=False (don't crash: the sweep needs this)
      - meets_duty = q >= required ;  fuel_ok = Tc_out <= fuel_limit
      - ok = valid and meets_duty and fuel_ok

    Returns dict: point name -> dict(result or None, valid, meets_duty, fuel_ok, ok).
    """
    # TODO
    results = {}
    for name, (m_fuel, T_fuel_in, Q_required) in mission.items():
        try:
            r = rate_plate_fin_auto(core, "turbine_oil", "jet_fuel", m_oil, m_fuel, T_oil_in, T_fuel_in)
            valid = True
        except ValueError:
            r = None
            valid = False

        meets_duty = r is not None and r["q"] >= Q_required
        fuel_ok = r is not None and r["Tc_out"] <= fuel_limit
        ok = valid and meets_duty and fuel_ok

        results[name] = dict(result=r, valid=valid, meets_duty=meets_duty,
                             fuel_ok=fuel_ok, ok=ok)
    return results

    


def is_feasible(mission_results):
    """True only if every point is ok."""
    return all(v["ok"] for v in mission_results.values())