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
    else:
        raise ValueError("Unknown fluid. Use 'water' or 'engine_oil'.")