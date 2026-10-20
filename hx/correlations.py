"""Convection and friction correlations for internal flow.

References:
    Incropera et al., Fundamentals of Heat and Mass Transfer, 7th ed., Ch. 8
"""
import numpy as np

RE_LAMINAR = 2300  # transition Reynolds number for tubes


def reynolds_from_mdot(m_dot, D_h, mu, A_flow):
    """Re = rho*u*D_h/mu = m_dot*D_h / (A_flow*mu)."""
    Re = m_dot * D_h / (A_flow * mu)
    return Re



def dittus_boelter(Re, Pr, heating=True):
    """Nu = 0.023 Re^0.8 Pr^n   (Incropera Eq. 8.60)
    n = 0.4 if the fluid is being heated, 0.3 if it is being cooled.
    Valid: Re > 10,000, 0.6 < Pr < 160.
    """
    if heating == True:
        n = 0.4
    else:
        n = 0.3
    Nu = 0.023 * Re**0.8 * Pr**n
    return Nu


def petukhov_friction(Re):
    """Darcy friction factor for smooth tubes (Incropera Eq. 8.21):
    f = (0.790 ln(Re) - 1.64)^-2,  valid 3000 < Re < 5e6.
    """
    if Re < 3000 or Re > 5e6:
        raise ValueError(f"Re = {Re} is out of range for petukhov_friction")
    else:
        f = (0.790 * np.log(Re) - 1.64) ** -2
        return f



def gnielinski(Re, Pr):
    """Nu = (f/8)(Re - 1000)Pr / (1 + 12.7 (f/8)^0.5 (Pr^(2/3) - 1))
    (Incropera Eq. 8.62), with f from petukhov_friction.
    Valid: 3000 < Re < 5e6, 0.5 < Pr < 2000.
    """
    if Re < 3000 or Re > 5e6 or Pr < 0.5 or Pr > 2000:
        raise ValueError(f"Re = {Re}, Pr = {Pr} are out of range for gnielinski")
    else:
        f = petukhov_friction(Re)
        Nu = (f / 8) * (Re - 1000) * Pr / (1 + 12.7 * (f / 8) ** 0.5 * (Pr ** (2 / 3) - 1))
        return Nu



# Incropera Table 8.2: laminar annulus, inner surface heated, outer insulated
_ANNULUS_RATIO = np.array([0.05, 0.10, 0.25, 0.50, 1.00])
_ANNULUS_NU_I = np.array([17.46, 11.56, 7.37, 5.74, 4.86])


def nu_laminar_annulus_inner(Di_over_Do):
    """Nu at the inner surface of an annulus, fully developed laminar flow.
    Linear interpolation of Incropera Table 8.2 (use np.interp).
    """
    if Di_over_Do < 0.05 or Di_over_Do > 1.0:
        raise ValueError(f"Di/Do = {Di_over_Do} is out of range for nu_laminar_annulus_inner")
    else:
        Nu = np.interp(Di_over_Do, _ANNULUS_RATIO, _ANNULUS_NU_I)
        return Nu


def nu_tube(Re, Pr, method="gnielinski", heating=True):
    """Pick the right correlation for a circular tube:
    - Re < RE_LAMINAR: laminar, Nu = 3.66 (constant wall temperature)
    - otherwise: Gnielinski or Dittus-Boelter, depending on `method`
    Raise ValueError for an unknown method.
    """
    if Re < RE_LAMINAR:
        Nu = 3.66
    elif method == "gnielinski":
        Nu = gnielinski(Re, Pr)
    elif method == "dittus_boelter":
        Nu = dittus_boelter(Re, Pr, heating)
    else:
        raise ValueError(f"Unknown method: {method}")
    return Nu   

def haaland_friction(Re, rel_roughness=0.0):
    X = -1.8 * np.log10((rel_roughness / 3.7) ** 1.11 + 6.9 / Re)
    return 1 / X**2
