"""Fluid property lookup.

- Water: CoolProp (IAPWS-IF97 formulation)
- Engine oil (unused): Incropera 7th ed., Table A.5, linear interpolation
All functions take temperature in degC and return SI units.
"""
import numpy as np
from CoolProp.CoolProp import PropsSI

# Incropera Table A.5, engine oil (unused), 273-430 K
_OIL_T_K = np.array([273, 280, 290, 300, 310, 320, 330, 340, 350,
                     360, 370, 380, 390, 400, 410, 420, 430])
_OIL_RHO = np.array([899.1, 895.3, 890.0, 884.1, 877.9, 871.8, 865.8, 859.9, 853.9,
                     847.8, 841.8, 836.0, 830.6, 825.1, 818.9, 812.1, 806.5])
_OIL_CP = 1e3 * np.array([1.796, 1.827, 1.868, 1.909, 1.951, 1.993, 2.035, 2.076, 2.118,
                          2.161, 2.206, 2.250, 2.294, 2.337, 2.381, 2.427, 2.471])
_OIL_MU = 1e-2 * np.array([385, 217, 99.9, 48.6, 25.3, 14.1, 8.36, 5.31, 3.56,
                           2.52, 1.86, 1.41, 1.10, 0.874, 0.698, 0.564, 0.470])
_OIL_K = 1e-3 * np.array([147, 144, 145, 145, 145, 143, 141, 139, 138,
                          138, 137, 136, 135, 134, 133, 133, 132])
# Mobil Jet Oil II: Walther (ASTM D341) fit to datasheet viscosities
_x40, _x100 = np.log10(40 + 273.15), np.log10(100 + 273.15)
_y40 = np.log10(np.log10(27.6 + 0.7))
_y100 = np.log10(np.log10(5.1 + 0.7))
_WALTHER_B = -(_y100 - _y40) / (_x100 - _x40)     # minus the slope -> positive
_WALTHER_A = _y40 + _WALTHER_B * _x40             # intercept


def water(T_C, P=101325.0):
    """Liquid water properties from CoolProp at T [degC], P [Pa].

    CoolProp needs T in kelvin. Keys: 'C' (cp), 'V' (viscosity),
    'L' (conductivity), 'D' (density).
    e.g. PropsSI('C', 'T', T_K, 'P', P, 'Water')

    Returns dict(rho=..., cp=..., mu=..., k=...).
    """
    T_K = T_C + 273.15
    rho = PropsSI('D', 'T', T_K, 'P', P, 'Water')
    cp = PropsSI('C', 'T', T_K, 'P', P, 'Water')
    mu = PropsSI('V', 'T', T_K, 'P', P, 'Water')
    k = PropsSI('L', 'T', T_K, 'P', P, 'Water')

    return dict(rho=rho, cp=cp, mu=mu, k=k)



def engine_oil(T_C):
    """Unused engine oil properties at T [degC], from Incropera Table A.5.

    Use np.interp on each property array. Raise ValueError if T is outside
    the table range (273-430 K); np.interp would silently clamp to the end values.

    Returns dict(rho=..., cp=..., mu=..., k=...).
    """
    T_K = T_C + 273.15
    if not (273 <= T_K <= 430):
        raise ValueError("Temperature out of range for engine oil properties.")

    rho = np.interp(T_K, _OIL_T_K, _OIL_RHO)
    cp = np.interp(T_K, _OIL_T_K, _OIL_CP)
    mu = np.interp(T_K, _OIL_T_K, _OIL_MU)
    k = np.interp(T_K, _OIL_T_K, _OIL_K)

    return dict(rho=rho, cp=cp, mu=mu, k=k)


def props(fluid, T_C):
    """Dispatch by name: 'water' or 'engine_oil'. Raise ValueError otherwise."""
    if fluid == "water":
        return water(T_C)
    elif fluid == "engine_oil":
        return engine_oil(T_C)
    elif fluid == "turbine_oil":
        return turbine_oil(T_C)
    elif fluid == "jet_fuel":
        return jet_fuel(T_C)
    else:
        raise ValueError("Unknown fluid. Use 'water' or 'engine_oil'.")

def turbine_oil(T_C):
    """Mobil Jet Oil II (MIL-PRF-23699) at T [degC].

    Viscosity, ASTM D341 (Walther):
        log10(log10(nu + 0.7)) = A - B*log10(T_K),   nu in cSt (mm2/s)
        Fit A and B to the two datasheet points (40 C: 27.6 cSt, 100 C: 5.1 cSt).
        -> compute A, B once at module level, not on every call
        Invert to get nu, then mu = nu * 1e-6 * rho.
    Density:  rho = 1003.5 * (1 - 7e-4 * (T_C - 15))
    k:        0.1496 - 1.58e-4 * T_C              (NIST IR 8263)
    cp:       reuse engine_oil(T_C)["cp"]          (stated assumption)

    Valid range: raise ValueError outside 0-200 C.
    Returns dict(rho, cp, mu, k).
    """
    # TODO
    if not (0 <= T_C <= 200):
        raise ValueError("Temperature out of range for turbine oil properties.")

    T_K = T_C + 273.15
    y = _WALTHER_A - _WALTHER_B * np.log10(T_K)
    nu_cSt = 10 ** (10 ** y) - 0.7

    rho = 1003.5 * (1 - 7e-4 * (T_C - 15))
    mu = nu_cSt * 1e-6 * rho
    cp = engine_oil(T_C)["cp"]
    k = 0.1496 - 1.58e-4 * T_C
    return dict(rho=rho, cp=cp, mu=mu, k=k)



def jet_fuel(T_C, P=5e5):
    """Jet A-1 surrogate: n-dodecane from CoolProp ('n-Dodecane').
    P defaults to 5 bar (fuel system pressure keeps it liquid).
    Same structure as water(). Raise ValueError below 0 C (surrogate freezes at -9.6 C).
    """
    # TODO
    if T_C < 0:
        raise ValueError("Temperature out of range for jet fuel properties.")
    else:
        T_K = T_C + 273.15
        rho = PropsSI('D', 'T', T_K, 'P', P, 'n-Dodecane')
        cp = PropsSI('C', 'T', T_K, 'P', P, 'n-Dodecane')
        mu = PropsSI('V', 'T', T_K, 'P', P, 'n-Dodecane')
        k = PropsSI('L', 'T', T_K, 'P', P, 'n-Dodecane')
    return dict(rho=rho, cp=cp, mu=mu, k=k)