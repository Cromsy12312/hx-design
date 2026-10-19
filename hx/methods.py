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