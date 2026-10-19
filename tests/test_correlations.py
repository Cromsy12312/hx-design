import numpy as np
from hx.correlations import (
    reynolds_from_mdot, dittus_boelter, petukhov_friction,
    gnielinski, nu_laminar_annulus_inner, nu_tube,
)

# Incropera Ex 11.1 water side: inner tube
m_c, Di, mu_c, k_c, Pr_c = 0.2, 0.025, 725e-6, 0.625, 4.85
# Oil side: annulus
m_h, Do, mu_h, k_h = 0.1, 0.045, 3.25e-2, 0.138


def test_water_reynolds():
    Re = reynolds_from_mdot(m_c, Di, mu_c, np.pi * Di**2 / 4)
    assert abs(Re - 14050) < 10


def test_water_dittus_boelter_matches_incropera():
    Re = reynolds_from_mdot(m_c, Di, mu_c, np.pi * Di**2 / 4)
    Nu = dittus_boelter(Re, Pr_c, heating=True)
    assert abs(Nu - 90) < 0.5
    h = Nu * k_c / Di
    assert abs(h - 2250) < 15


def test_gnielinski_close_to_dittus_boelter():
    """Gnielinski is more accurate; for this case it is about 4% higher."""
    Nu_g = gnielinski(14050, Pr_c)
    assert abs(Nu_g - 93.8) < 1.0


def test_petukhov_value():
    assert np.isclose(petukhov_friction(14050), 0.0287, atol=2e-4)


def test_oil_annulus_reynolds_and_nu():
    A = np.pi * (Do**2 - Di**2) / 4
    Re = reynolds_from_mdot(m_h, Do - Di, mu_h, A)
    assert abs(Re - 56.0) < 0.5
    Nu = nu_laminar_annulus_inner(Di / Do)
    assert abs(Nu - 5.63) < 0.05      # book interpolates to 5.63
    h = Nu * k_h / (Do - Di)
    assert abs(h - 38.8) < 0.5


def test_nu_tube_switches_to_laminar():
    assert nu_tube(1000, 5.0) == 3.66
    assert nu_tube(20000, 5.0) > 3.66