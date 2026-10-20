"""Double-pipe (concentric tube) heat exchanger sizing.

Thin-walled inner tube assumed (wall resistance neglected), as in Incropera Ex 11.1.
"""
from dataclasses import dataclass
import numpy as np

from hx.methods import heat_capacity_rate, lmtd
from hx.correlations import (
    reynolds_from_mdot, nu_tube, nu_laminar_annulus_inner, RE_LAMINAR,
)


@dataclass
class Stream:
    """One fluid stream, with properties at its mean temperature."""
    m_dot: float  # kg/s
    cp: float     # J/kg.K
    mu: float     # Pa.s
    k: float      # W/m.K
    T_in: float   # degC

    @property
    def Pr(self):
        return self.cp * self.mu / self.k


def size_double_pipe(hot, cold, Th_out, Di, Do, hot_in="annulus",
                     arrangement="counterflow", method="gnielinski",
                     Rf_i=0.0, Rf_o=0.0):
    """Required tube length to cool `hot` to Th_out.

    Steps:
      1. q from the hot-side energy balance; Tc_out from the cold side
      2. LMTD
      3. Re, Nu, h for the tube side and the annulus side
         - tube: area = pi Di^2/4,  D_h = Di,  use nu_tube
         - annulus: area = pi(Do^2 - Di^2)/4,  D_h = Do - Di
           laminar -> nu_laminar_annulus_inner(Di/Do); turbulent -> nu_tube
         - the tube-side fluid is heated if it is the cold stream
      4. U = 1 / (1/h_i + Rf_i + Rf_o + 1/h_o)   (thin wall)
      5. L = q / (U * pi * Di * LMTD)

    Returns a dict with q, Tc_out, LMTD, Re_i, Re_o, h_i, h_o, U, L.
    """
    q = hot.m_dot * hot.cp * (hot.T_in - Th_out)
    tc_out = cold.T_in + q / (cold.m_dot * cold.cp)
    LMTD = lmtd(Th_in=hot.T_in, Th_out=Th_out, Tc_in=cold.T_in, Tc_out=tc_out, arrangement=arrangement)
    tube_area = np.pi * Di**2 / 4
    annulus_area = np.pi * (Do**2 - Di**2) / 4
    Dh_tube = Di
    Dh_annulus = Do - Di
    if hot_in == "annulus":
        tube_fluid, ann_fluid = cold, hot
    elif hot_in == "tube":
        tube_fluid, ann_fluid = hot, cold
    else:
        raise ValueError(f"hot_in must be 'tube' or 'annulus', got {hot_in!r}")
    tube_heated = tube_fluid is cold   # the cold stream is the one being heated

    Re_i = reynolds_from_mdot(tube_fluid.m_dot, Dh_tube, tube_fluid.mu, tube_area)
    Re_o = reynolds_from_mdot(ann_fluid.m_dot, Dh_annulus, ann_fluid.mu, annulus_area)

    Nu_i = nu_tube(Re_i, tube_fluid.Pr, method=method, heating=tube_heated)
    if Re_o < RE_LAMINAR:
        Nu_o = nu_laminar_annulus_inner(Di / Do)
    else:
        Nu_o = nu_tube(Re_o, ann_fluid.Pr, method=method, heating=not tube_heated)

    h_i = Nu_i * tube_fluid.k / Dh_tube
    h_o = Nu_o * ann_fluid.k / Dh_annulus
    U = 1 / (1 / h_i + Rf_i + Rf_o + 1 / h_o)
    L = q / (U * np.pi * Di * LMTD)
    return dict(q=q, Tc_out=tc_out, LMTD=LMTD, Re_i=Re_i, Re_o=Re_o,
                h_i=h_i, h_o=h_o, U=U, L=L)
