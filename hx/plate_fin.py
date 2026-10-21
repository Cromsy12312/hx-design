"""Compact plate-fin heat exchanger with offset strip fins.

References:
    Manglik & Bergles (1995), Exp. Thermal Fluid Sci. 10, 171-180
    Incropera et al., 7th ed., Sec. 3.6 (fin efficiency)
"""
from dataclasses import dataclass
import numpy as np

from hx.methods import effectiveness_crossflow_unmixed, heat_capacity_rate


@dataclass
class OffsetStripFin:
    s: float   # clear fin spacing [m]
    h: float   # clear fin height [m]
    t: float   # fin thickness [m]
    l: float   # strip length in flow direction [m]


    @property
    def alpha(self):
        return self.s / self.h

    @property
    def delta(self):
        return self.t / self.l

    @property
    def gamma(self):
        return self.t / self.s

    @property
    def D_h(self):
        """Hydraulic diameter (Manglik & Bergles):
        D_h = 4 s h l / ( 2(s l + h l + t h) + t s )
        """
        D_h = 4 * self.s * self.h * self.l / (2 * (self.s * self.l + self.h * self.l + self.t * self.h) + self.t * self.s)
        return D_h


def j_manglik_bergles(Re, fin):
    """Colburn j factor:
    j = 0.6522 Re^-0.5403 alpha^-0.1541 delta^0.1499 gamma^-0.0678
        * [1 + 5.269e-5 Re^1.340 alpha^0.504 delta^0.456 gamma^-1.055]^0.1

    Raise ValueError outside the fitted range 120 <= Re <= 10,000.
    """
    # TODO
    if not (120 <= Re <= 10000):
        raise ValueError(f"Re = {Re} is out of range for j_manglik_bergles")
    else:
        j = 0.6522 * Re**-0.5403 * fin.alpha**-0.1541 * fin.delta**0.1499 * fin.gamma**-0.0678 * (1 + 5.269e-5 * Re**1.340 * fin.alpha**0.504 * fin.delta**0.456 * fin.gamma**-1.055)**0.1
        return j
    


def f_manglik_bergles(Re, fin):
    """Fanning friction factor:
    f = 9.6243 Re^-0.7422 alpha^-0.1856 delta^0.3053 gamma^-0.2659
        * [1 + 7.669e-8 Re^4.429 alpha^0.920 delta^3.767 gamma^0.236]^0.1

    Same validity check as above.
    """
    if not (120 <= Re <= 10000):
        raise ValueError(f"Re = {Re} is out of range for f_manglik_bergles")
    else:
        f = 9.6243 * Re**-0.7422 * fin.alpha**-0.1856 * fin.delta**0.3053 * fin.gamma**-0.2659 * (1 + 7.669e-8 * Re**4.429 * fin.alpha**0.920 * fin.delta**3.767 * fin.gamma**0.236)**0.1
        return f

def fin_efficiency(h_conv, k_fin, t, L_fin):
    """Straight fin, adiabatic tip (Incropera Table 3.4):
    m = sqrt(2 h / (k t))
    eta_f = tanh(m L) / (m L)

    In a plate-fin core each fin spans between two plates and is heated from
    both ends, so the effective fin length is HALF the fin height: L_fin = h/2.
    """
    m = np.sqrt(2 * h_conv / (k_fin * t))
    eta_f = np.tanh(m * L_fin) / (m * L_fin)
    return eta_f


def surface_efficiency(eta_f, fin_area_fraction):
    """Overall surface efficiency (Incropera Eq. 11.3):
    eta_o = 1 - (A_fin / A_total) * (1 - eta_f)

    For a plate-fin channel, A_fin/A_total ~ h / (h + s)
    (two fin walls of height h vs two plate faces of width s per channel).
    """
    eta_o = 1 - fin_area_fraction * (1 - eta_f)
    return eta_o


@dataclass
class PlateFinCore:
    L_hot: float              # hot-side flow length [m]
    L_cold: float             # cold-side flow length [m]
    N_hot: int                # number of hot layers
    fin_hot: OffsetStripFin
    fin_cold: OffsetStripFin
    t_plate: float = 0.5e-3   # parting plate thickness [m]
    k_metal: float = 170.0    # aluminium alloy [W/m.K]
    rho_metal: float = 2700.0 # [kg/m3]

    # TODO: @property N_cold = N_hot + 1
    @property
    def N_cold(self):
        return self.N_hot + 1


def side_geometry(core, side):
    """Geometry of one side ('hot' or 'cold').

    For the hot side: N = N_hot, fin = fin_hot, L_flow = L_hot, width = L_cold
    For the cold side: swap all of those.

    channels per layer  n_ch = floor(width / (s + t))
    free-flow area      A_ff = N * n_ch * s * h
    heat transfer area  A    = 4 * A_ff * L_flow / D_h     (rearranged definition of D_h)
    fin area fraction        = h / (h + s)

    Raise ValueError for an unknown side.
    Returns dict(N, n_ch, A_ff, A, L_flow, fin_area_fraction, fin).
    """
    if side == 'hot':
        N = core.N_hot
        fin = core.fin_hot
        L_flow = core.L_hot
        width = core.L_cold
    elif side == 'cold':
        N = core.N_cold
        fin = core.fin_cold
        L_flow = core.L_cold
        width= core.L_hot
    else:
        raise ValueError(f"Unknown side: {side!r}")
    n_ch = int(np.floor(width / (fin.s + fin.t)))
    A_ff = N * n_ch * fin.s * fin.h
    A = 4 * A_ff * L_flow / fin.D_h
    fin_area_fraction = fin.h / (fin.h + fin.s)
    return dict(N=N, n_ch=n_ch, A_ff=A_ff, A=A, L_flow=L_flow, fin_area_fraction=fin_area_fraction, fin=fin)



def side_performance(core, side, fluid):
    """Film coefficient, surface efficiency and pressure drop for one side.

    G   = m_dot / A_ff                 (mass velocity in the fin passages)
    Re  = G D_h / mu
    h   = j G cp Pr^(-2/3)             (definition of the Colburn j factor)
    eta_f from fin_efficiency(h, k_metal, t, L_fin = fin.h / 2)
    eta_o from surface_efficiency
    dP  = 4 f (L_flow / D_h) G^2 / (2 rho)      (Fanning f, core friction only)

    Returns dict(G, Re, h, eta_f, eta_o, dP) plus everything from side_geometry.
    """
    G = fluid.m_dot / side_geometry(core, side)['A_ff']
    Re = G * side_geometry(core, side)['fin'].D_h / fluid.mu
    h = j_manglik_bergles(Re, side_geometry(core, side)['fin']) * G * fluid.cp * fluid.Pr**(-2/3)
    eta_f = fin_efficiency(h, core.k_metal, side_geometry(core, side)['fin'].t, side_geometry(core, side)['fin'].h / 2)
    eta_o = surface_efficiency(eta_f, side_geometry(core, side)['fin_area_fraction'])
    dP = 4 * f_manglik_bergles(Re, side_geometry(core, side)['fin']) * (side_geometry(core, side)['L_flow'] / side_geometry(core, side)['fin'].D_h) * G**2 / (2 * fluid.rho)
    return dict(G=G, Re=Re, h=h, eta_f=eta_f, eta_o=eta_o, dP=dP)


def core_mass(core):
    """Metal mass [kg] (side bars and headers neglected).
    plates: (N_hot + N_cold + 1) * t_plate * L_hot * L_cold
    fins, per side: N * t * (h + s)/(s + t) * width * L_flow
        (the folded fin sheet is ~ (h + s)/(s + t) times longer than the width it covers)
    mass = rho_metal * (plates + fins)
    """
    plates = (core.N_hot + core.N_cold + 1) * core.t_plate * core.L_hot * core.L_cold
    fins_hot = core.N_hot * core.fin_hot.t * (core.fin_hot.h + core.fin_hot.s) / (core.fin_hot.s + core.fin_hot.t) * core.L_cold * core.L_hot
    fins_cold = core.N_cold * core.fin_cold.t * (core.fin_cold.h + core.fin_cold.s) / (core.fin_cold.s + core.fin_cold.t) * core.L_hot * core.L_cold
    mass = core.rho_metal * (plates + fins_hot + fins_cold)
    return mass


def core_volume(core):
    """Envelope volume [m3]:
    stack height = N_hot (h_hot + t_hot) + N_cold (h_cold + t_cold) + (N_hot + N_cold + 1) t_plate
    volume = L_hot * L_cold * stack height
    """
    stack_height = core.N_hot * (core.fin_hot.h + core.fin_hot.t) + core.N_cold * (core.fin_cold.h + core.fin_cold.t) + (core.N_hot + core.N_cold + 1) * core.t_plate
    volume = core.L_hot * core.L_cold * stack_height
    return volume


def rate_plate_fin(core, hot, cold):
    """Rating of a given core (properties fixed).

    1. side_performance for both sides
    2. 1/UA = 1/(eta_o,h h_h A_h) + 1/(eta_o,c h_c A_c)    (thin, highly conductive plates: wall neglected)
    3. Cmin, Cmax, NTU = UA/Cmin, eps = effectiveness_crossflow_unmixed
    4. q, Th_out, Tc_out (same as rate_double_pipe)
    5. mass, volume, and q per kg

    Returns dict(q, Th_out, Tc_out, UA, NTU, eps, dP_hot, dP_cold,
                 Re_hot, Re_cold, mass, volume, kW_per_kg).
    """
    hot_perf = side_performance(core, 'hot', hot)
    cold_perf = side_performance(core, 'cold', cold)
    RA = 1/((hot_perf['eta_o'] * hot_perf['h'] * side_geometry(core, 'hot')['A'])) + 1/(cold_perf['eta_o'] * cold_perf['h'] * side_geometry(core, 'cold')['A'])
    UA = 1 / RA
    Ch = heat_capacity_rate(hot.m_dot, hot.cp) 
    Cc = heat_capacity_rate(cold.m_dot, cold.cp)
    Cmin, Cmax = min(Ch, Cc), max(Ch, Cc)
    NTU = UA / Cmin
    eps = effectiveness_crossflow_unmixed(NTU, Cmin / Cmax)
    q = eps * Cmin * (hot.T_in - cold.T_in)
    Th_out = hot.T_in - q / Ch                      
    Tc_out = cold.T_in + q / Cc
    mass = core_mass(core)
    volume = core_volume(core)
    kW_per_kg = q / mass * 1e-3
    return dict(q=q, Th_out=Th_out, Tc_out=Tc_out, UA=UA, NTU=NTU, eps=eps, dP_hot=hot_perf['dP'], dP_cold=cold_perf['dP'],
                Re_hot=hot_perf['Re'], Re_cold=cold_perf['Re'], mass=mass, volume=volume, kW_per_kg=kW_per_kg)
