import numpy as np
from hx.double_pipe import Stream
from hx.shell_tube import (ShellTubeGeometry, bundle_diameter, tube_side,
                           shell_side_kern, overall_U_shell_tube)

# Sinnott Ex 12.1, book property values
methanol = Stream(m_dot=100000 / 3600, cp=2840, mu=0.34e-3, k=0.19, T_in=95, rho=750)
water = Stream(m_dot=68.9, cp=4200, mu=0.8e-3, k=0.59, T_in=25, rho=995)

geom = ShellTubeGeometry(do=0.020, di=0.016, L=4.83, Nt=918, n_passes=2,
                         Ds=0.894, baffle_spacing=0.178)


def test_bundle_diameter():
    assert abs(bundle_diameter(918, 0.020, 2) - 0.826) < 0.002


def test_tube_side():
    t = tube_side(geom, water, heating=True)
    assert abs(t["A_flow"] - 0.092) < 0.001
    assert abs(t["u"] - 0.75) < 0.01
    assert abs(t["Re"] - 14925) / 14925 < 0.005
    assert abs(t["h_i"] - 3812) / 3812 < 0.03      # Gnielinski vs book chart


def test_shell_side():
    s = shell_side_kern(geom, methanol)
    assert abs(s["As"] - 0.032) < 0.0005
    assert abs(s["de"] - 0.0142) < 0.0001          # book prints 14.4 mm (arithmetic slip)
    assert abs(s["h_s"] - 2740) / 2740 < 0.03      # Kern correlation vs book chart


def test_overall_U_matches_sinnott():
    t = tube_side(geom, water, heating=True)
    s = shell_side_kern(geom, methanol)
    U = overall_U_shell_tube(s["h_s"], t["h_i"], geom, R_fo=1 / 5000, R_fi=1 / 3000)
    assert abs(U - 738) / 738 < 0.01


def test_book_h_values_reproduce_book_U():
    """Feeding the book's own h values in must give exactly the book's U."""
    U = overall_U_shell_tube(2740, 3812, geom, R_fo=1 / 5000, R_fi=1 / 3000)
    assert abs(U - 738) < 1