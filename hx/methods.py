import numpy as np

def heat_capacity_rate(m_dot, cp):
    """C = m_dot * cp  [W/K]"""
    C = m_dot * cp
    return C

def lmtd(Th_in, Th_out, Tc_in, Tc_out, arrangement="counterflow"):
    """Log mean temperature difference [K].

    counterflow: dT1 = Th_in - Tc_out,  dT2 = Th_out - Tc_in
    parallel:    dT1 = Th_in - Tc_in,   dT2 = Th_out - Tc_out
    LMTD = (dT1 - dT2) / ln(dT1 / dT2)

    Edge case: if dT1 == dT2 (balanced counterflow), LMTD = dT1.
    """
    if arrangement == "counterflow":
        dT1, dT2 = Th_in - Tc_out, Th_out - Tc_in
    elif arrangement == "parallel":
        dT1, dT2 = Th_in - Tc_in, Th_out - Tc_out
    else:
        raise ValueError(f"Unknown arrangement: {arrangement!r}")

    if np.isclose(dT1, dT2):
        return dT1
    return (dT1 - dT2) / np.log(dT1 / dT2)


def effectiveness(NTU, Cr, arrangement="counterflow"):
    """Effectiveness eps from NTU and Cr = Cmin/Cmax  (Incropera Table 11.3).

    counterflow:  eps = (1 - exp(-NTU(1-Cr))) / (1 - Cr*exp(-NTU(1-Cr)))
                  if Cr == 1: eps = NTU / (1 + NTU)
    parallel:     eps = (1 - exp(-NTU(1+Cr))) / (1 + Cr)
    Any arrangement with Cr == 0: eps = 1 - exp(-NTU)
    """
    if arrangement == "parallel":
        eps = (1 - np.exp(-NTU * (1 + Cr))) / (1 + Cr)
    elif arrangement == "counterflow":
        if np.isclose(Cr, 1):
            eps = NTU / (1 + NTU)
        else:
            eps = (1 - np.exp(-NTU * (1 - Cr))) / (1 - Cr * np.exp(-NTU * (1 - Cr)))
    else:
        raise ValueError(f"Unknown arrangement: {arrangement!r}")
    return eps

def ntu_from_effectiveness(eps, Cr, arrangement="counterflow"):
    """Inverse of the above: NTU needed to reach a given eps (Incropera Table 11.4).

    counterflow:  NTU = 1/(Cr-1) * ln((eps-1)/(eps*Cr-1))
                  if Cr == 1: NTU = eps / (1 - eps)
    parallel:     NTU = -ln(1 - eps(1+Cr)) / (1+Cr)
    """

    if arrangement == "counterflow":
        if np.isclose(Cr, 1):
            NTU = eps / (1 - eps)
        else:
            NTU = 1/(Cr - 1) * np.log((eps - 1) / (eps * Cr - 1)) 
    elif arrangement == "parallel":
        NTU = -np.log(1 - eps * (1 + Cr)) / (1 + Cr) #for parallel arrangement
    else:
        raise ValueError(f"Unknown arrangement: {arrangement!r}")
    return NTU

def f_factor_1_2(R, S):
    """LMTD correction factor F for 1 shell pass, 2 (or 2n) tube passes.

    R = (Th_in - Th_out) / (Tc_out - Tc_in)
    S = (Tc_out - Tc_in) / (Th_in - Tc_in)

    F = sqrt(R^2+1) * ln[(1-S)/(1-R*S)]
        / ( (R-1) * ln[ (2 - S(R+1-sqrt(R^2+1))) / (2 - S(R+1+sqrt(R^2+1))) ] )

    Special case R == 1 (formula is 0/0):
    F = (sqrt(2) * S/(1-S)) / ln[ (2 - S(2-sqrt(2))) / (2 - S(2+sqrt(2))) ]

    Raise ValueError if the logs' arguments are <= 0: that combination of
    temperatures is impossible in a single 1-2 shell (needs more shells).
    """
    R = float(R)
    S = float(S)
    if R == 1:
        if S >= 1:
            raise ValueError(f"Impossible combination of temperatures: R={R}, S={S}")
        else: 
            F = (np.sqrt(2) * S / (1 - S)) / np.log((2 - S * (2 - np.sqrt(2))) / (2 - S * (2 + np.sqrt(2))))
    else:
        if S >= 1 or (1 - R * S) <= 0 or (2 - S * (R + 1 - np.sqrt(R**2 + 1))) <= 0 or (2 - S * (R + 1 + np.sqrt(R**2 + 1))) <= 0:
            raise ValueError(f"Impossible combination of temperatures: R={R}, S={S}")
        else:
            F = np.sqrt(R**2 + 1) * np.log((1 - S) / (1 - R * S)) / ((R - 1) * np.log((2 - S * (R + 1 - np.sqrt(R**2 + 1))) / (2 - S * (R + 1 + np.sqrt(R**2 + 1)))))
    
    return F


def effectiveness_shell_1_2(NTU, Cr):
    """eps for one shell pass, 2n tube passes (Incropera Eq. 11.30a):

    eps = 2 / ( 1 + Cr + sqrt(1+Cr^2) * (1 + exp(-NTU*sqrt(1+Cr^2)))
                                       / (1 - exp(-NTU*sqrt(1+Cr^2))) )
    """
    eps = 2 / (1 + Cr + np.sqrt(1 + Cr**2) * (1 + np.exp(-NTU * np.sqrt(1 + Cr**2))) / (1 - np.exp(-NTU * np.sqrt(1 + Cr**2))))
    return eps


def ntu_shell_1_2(eps, Cr):
    """Inverse of the above (Incropera Eq. 11.31a/b):

    E   = (2/eps - (1+Cr)) / sqrt(1+Cr^2)
    NTU = -ln((E-1)/(E+1)) / sqrt(1+Cr^2)

    Raise ValueError if E <= 1 (eps too high for a single 1-2 shell).
    """
    E = (2 / eps - (1 + Cr)) / np.sqrt(1 + Cr**2)
    if E <= 1:
        raise ValueError(f"eps={eps} is too high for a single 1-2 shell with Cr={Cr}")
    NTU = -np.log((E - 1) / (E + 1)) / np.sqrt(1 + Cr**2)
    return NTU


# TODO TOMMOROW
#
#
#
def effectiveness_crossflow_unmixed(NTU, Cr):
    """Crossflow, both fluids unmixed (Incropera Eq. 11.32)."""
    if np.isclose(Cr, 0):
        eps = 1 - np.exp(-NTU)
    else:
        inner = np.exp(-Cr * NTU**0.78) - 1         
        eps = 1 - np.exp((1 / Cr) * NTU**0.22 * inner) 
    return eps


def ntu_crossflow_unmixed(eps, Cr):
    """Inverse of the above, found numerically.

    Hint: scipy.optimize.brentq(f, a, b) finds x in [a, b] where f(x) = 0,
    as long as f(a) and f(b) have opposite signs.
      - write a small function of NTU that is zero at the answer
      - choose a bracket, e.g. a = 1e-6, b = 50
      - raise ValueError if eps is unreachable (f(b) still < 0)
    """
    # TODO
    def f(NTU):
        return effectiveness_crossflow_unmixed(NTU, Cr) - eps
    try:
        from scipy.optimize import brentq
        NTU = brentq(f, 1e-6, 50)
    except ValueError:
        raise ValueError(f"eps={eps} is unreachable for Cr={Cr}")
    return NTU