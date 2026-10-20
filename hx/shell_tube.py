"""Shell-and-tube heat exchanger: Kern method (1 shell pass, 2n tube passes).

References:
    Sinnott, Coulson & Richardson's Chemical Engineering Vol. 6, 4th ed., Sec. 12.9.3
"""
from dataclasses import dataclass
import numpy as np

from hx.correlations import haaland_friction, nu_tube

# Sinnott Table 12.4: bundle diameter constants, pitch = 1.25 do
# {layout: {tube passes: (K1, n1)}}
_BUNDLE_K = {
    "triangular": {1: (0.319, 2.142), 2: (0.249, 2.207), 4: (0.175, 2.285),
                   6: (0.0743, 2.499), 8: (0.0365, 2.675)},
    "square":     {1: (0.215, 2.207), 2: (0.156, 2.291), 4: (0.158, 2.263),
                   6: (0.0402, 2.617), 8: (0.0331, 2.643)},
}


@dataclass
class ShellTubeGeometry:
    do: float              # tube outside diameter [m]
    di: float              # tube inside diameter [m]
    L: float               # effective tube length [m]
    Nt: int                # number of tubes
    n_passes: int          # tube passes
    Ds: float              # shell inside diameter [m]
    baffle_spacing: float  # [m]
    pitch_ratio: float = 1.25
    layout: str = "triangular"
    k_wall: float = 50.0   # tube wall conductivity [W/m.K]

    @property
    def pitch(self):
        return self.pitch_ratio * self.do


def bundle_diameter(Nt, do, n_passes, layout="triangular"):
    """Db = do * (Nt / K1)^(1/n1)   (Sinnott Eq. 12.3b), constants from _BUNDLE_K."""
    Db = do * (Nt / _BUNDLE_K[layout][n_passes][0]) ** (1 / _BUNDLE_K[layout][n_passes][1])
    return Db



def tube_side(geom, fluid, method="gnielinski", heating=True):
    """Tube-side film coefficient.
    flow area = (Nt / n_passes) * pi di^2 / 4 ;  u = m_dot / (rho * area)
    Re = rho u di / mu ;  h = Nu k / di  (Nu from nu_tube)
    Returns dict(A_flow, u, Re, Pr, h_i).
    """
    flow_area = (geom.Nt / geom.n_passes) * np.pi * geom.di**2 / 4
    u = fluid.m_dot / (fluid.rho * flow_area)
    Re = fluid.rho * u * geom.di / fluid.mu
    h_i = nu_tube(Re, fluid.Pr, method=method, heating=heating) * fluid.k / geom.di
    return dict(A_flow=flow_area, u=u, Re=Re, Pr=fluid.Pr, h_i=h_i)


def shell_side_kern(geom, fluid, mu_wall=None):
    """Kern method shell-side coefficient (Sinnott Sec. 12.9.3).
    As = (pt - do) * Ds * lB / pt                         (12.21)
    Gs = m_dot / As ;  us = Gs / rho
    de: triangular  1.10/do * (pt^2 - 0.917 do^2)          (12.23)
        square      1.27/do * (pt^2 - 0.785 do^2)          (12.22)
    Re = Gs de / mu ;  Pr = cp mu / k
    Nu = 0.36 Re^0.55 Pr^(1/3) (mu/mu_wall)^0.14   (Kern's correlation;
         replaces reading jh off Fig. 12.29; valid ~2e3 < Re < 1e6)
    If mu_wall is None, take the viscosity correction as 1.
    Returns dict(As, Gs, us, de, Re, Pr, h_s).
    """
    As = (geom.pitch - geom.do) * geom.Ds * geom.baffle_spacing / geom.pitch
    Gs = fluid.m_dot / As
    us = Gs / fluid.rho
    if geom.layout == "triangular":
        de = 1.10 / geom.do * (geom.pitch**2 - 0.917 * geom.do**2)
    elif geom.layout == "square":
        de = 1.27 / geom.do * (geom.pitch**2 - 0.785 * geom.do**2)
    else:
        raise ValueError(f"layout must be 'triangular' or 'square', got {geom.layout!r}")
    Re = Gs * de / fluid.mu
    Pr = fluid.cp * fluid.mu / fluid.k
    if mu_wall is None:
        mu_ratio = 1.0
    else:
        mu_ratio = fluid.mu / mu_wall
    Nu = 0.36 * Re**0.55 * Pr**(1 / 3) * mu_ratio**0.14
    h_s = Nu * fluid.k / de
    return dict(As=As, Gs=Gs, us=us, de=de, Re=Re, Pr=Pr, h_s=h_s)



def overall_U_shell_tube(h_s, h_i, geom, R_fo=0.0, R_fi=0.0):
    """U based on tube OUTSIDE area (Sinnott Eq. 12.2), with wall resistance:
    1/Uo = 1/h_s + R_fo + do*ln(do/di)/(2 k_w) + (do/di)*R_fi + (do/di)/h_i
    R_fo, R_fi are fouling RESISTANCES [m2.K/W] (= 1 / fouling coefficient).
    """
    R_wall = geom.do * np.log(geom.do / geom.di) / (2 * geom.k_wall)
    Uo = 1 / (1 / h_s + R_fo + R_wall+ (geom.do / geom.di) * R_fi + (geom.do / geom.di) / h_i)
    return Uo

def jf_shell_kern(Re, baffle_cut=25):
    """Shell-side friction factor, Sinnott Fig. 12.30, 25% baffle cut.
    Power-law fit digitised from the chart: jf = 0.183 Re^-0.142 (1e4 < Re < 1e6).
    Raise ValueError for any other baffle cut (not digitised).
    """
    if baffle_cut != 25:
        raise ValueError(f"Only 25% baffle cut is implemented, got {baffle_cut}")
    else:
        if Re < 1e4 or Re > 1e6:
            raise ValueError(f"Re = {Re} is out of range for jf_shell_kern")
        else:
            jf = 0.183 * Re**-0.142
    return jf



def tube_side_dp(geom, fluid, tube, roughness=0.0):
    """Tube-side pressure drop (Sinnott Eq. 12.20, viscosity correction neglected):
    dP = n_passes * ( f * L/di + 2.5 ) * rho * u^2 / 2
    f = Darcy friction factor from haaland_friction(Re, roughness/di)
    (Sinnott writes 8*jf; 8*jf is the Darcy f.) The 2.5 is the return losses per pass.
    `tube` is the dict returned by tube_side().
    """
    dP = geom.n_passes * (haaland_friction(tube["Re"], roughness / geom.di) * geom.L / geom.di + 2.5) * fluid.rho * tube["u"]**2 / 2
    return dP
  


def shell_side_dp(geom, fluid, shell, baffle_cut=25):
    """Shell-side pressure drop (Sinnott Eq. 12.26, viscosity correction neglected):
    dP = 8 * jf * (Ds/de) * (L/lB) * rho * us^2 / 2
    `shell` is the dict returned by shell_side_kern().
    """
    dP = 8 * jf_shell_kern(shell["Re"], baffle_cut) * (geom.Ds / shell["de"]) * (geom.L / geom.baffle_spacing) * fluid.rho * shell["us"]**2 / 2
    return dP
