"""Double-pipe (concentric tube) heat exchanger sizing.

Thin-walled inner tube assumed (wall resistance neglected), as in Incropera Ex 11.1.
"""
from dataclasses import dataclass
import numpy as np

from hx.methods import effectiveness, heat_capacity_rate, lmtd
from hx.correlations import (
    reynolds_from_mdot, nu_tube, nu_laminar_annulus_inner, RE_LAMINAR,
)
from hx.properties import props


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

def _overall_U(hot, cold, Di, Do, hot_in, method, Rf_i, Rf_o):
    """Tube/annulus film coefficients and thin-wall overall U.
    Returns dict with Re_i, Re_o, h_i, h_o, U."""
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
    return dict(Re_i=Re_i, Re_o=Re_o, h_i=h_i, h_o=h_o, U=U)


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
    film = _overall_U(hot, cold, Di, Do, hot_in, method, Rf_i, Rf_o)
    L = q / (film["U"] * np.pi * Di * LMTD)
    return dict(q=q, Tc_out=tc_out, LMTD=LMTD, L=L, **film)

def rate_double_pipe(hot, cold, L, Di, Do, hot_in="annulus",
                     arrangement="counterflow", method="gnielinski",
                     Rf_i=0.0, Rf_o=0.0):
    """Outlet temperatures and duty for an exchanger of known length L.

    Steps:
      1. U from _overall_U  (properties are fixed, so U doesn't need the outlets)
      2. A = pi * Di * L
      3. Ch, Cc; Cmin, Cmax; Cr = Cmin/Cmax
      4. NTU = U*A/Cmin  ->  eps = effectiveness(NTU, Cr, arrangement)
      5. q = eps * Cmin * (Th_in - Tc_in)
      6. Th_out = Th_in - q/Ch ;  Tc_out = Tc_in + q/Cc

    Returns dict with q, Th_out, Tc_out, NTU, eps, U.
    """
    film = _overall_U(hot, cold, Di, Do, hot_in, method, Rf_i, Rf_o)
    A = np.pi * Di * L
    Ch = heat_capacity_rate(hot.m_dot, hot.cp)
    Cc = heat_capacity_rate(cold.m_dot, cold.cp)
    C_min, C_max = min(Ch, Cc), max(Ch, Cc)
    Cr = C_min / C_max
    NTU = film["U"] * A / C_min
    eps = effectiveness(NTU, Cr, arrangement)
    q = eps * C_min * (hot.T_in - cold.T_in)
    Th_out = hot.T_in - q / Ch
    Tc_out = cold.T_in + q / Cc
    return dict(q=q, Th_out=Th_out, Tc_out=Tc_out, NTU=NTU, eps=eps, U=film["U"])

def rate_double_pipe_auto(hot_fluid, cold_fluid, m_hot, m_cold, Th_in, Tc_in,
                          L, Di, Do, tol=0.01, max_iter=50, **kwargs):
    """Rating with properties evaluated at each stream's mean temperature.

    hot_fluid / cold_fluid are names for props(), e.g. 'engine_oil', 'water'.
    **kwargs are passed straight on to rate_double_pipe (hot_in, method, ...).

    Loop:
      - build Stream objects from props(fluid, T_mean)
      - r = rate_double_pipe(...)
      - new means from the inlet and r's outlets
      - stop when both outlets change by < tol; raise RuntimeError after max_iter

    Returns r plus 'iterations'.
    """
    Th_mean, Tc_mean = Th_in, Tc_in        # first guess: inlet temperatures
    Th_out_old, Tc_out_old = None, None

    for i in range(1, max_iter + 1):
        ph = props(hot_fluid, Th_mean)
        pc = props(cold_fluid, Tc_mean)
        hot = Stream(m_dot=m_hot, cp=ph["cp"], mu=ph["mu"], k=ph["k"], T_in=Th_in)
        cold = Stream(m_dot=m_cold, cp=pc["cp"], mu=pc["mu"], k=pc["k"], T_in=Tc_in)

        r = rate_double_pipe(hot, cold, L, Di, Do, **kwargs)

        if Th_out_old is not None and Tc_out_old is not None:
            if abs(r["Th_out"] - Th_out_old) < tol and abs(r["Tc_out"] - Tc_out_old) < tol:
                r["iterations"] = i
                return r

        Th_mean = (Th_in + r["Th_out"]) / 2
        Tc_mean = (Tc_in + r["Tc_out"]) / 2
        Th_out_old, Tc_out_old = r["Th_out"], r["Tc_out"]


    raise RuntimeError(f"Property iteration did not converge in {max_iter} iterations")