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

from dataclasses import replace
from hx.shell_tube import jf_shell_kern, tube_side_dp, shell_side_dp


def test_jf_shell_matches_book_reading():
    assert abs(jf_shell_kern(36762) - 0.04) / 0.04 < 0.05


def test_tube_dp_rough_close_to_book():
    t = tube_side(geom, water)
    dp_rough = tube_side_dp(geom, water, t, roughness=45e-6)
    dp_smooth = tube_side_dp(geom, water, t, roughness=0.0)
    assert abs(dp_rough - 7211) / 7211 < 0.06
    assert dp_smooth < dp_rough


def test_shell_dp_close_to_book():
    s = shell_side_kern(geom, methanol)
    assert abs(shell_side_dp(geom, methanol, s) - 272e3) / 272e3 < 0.07


def test_doubling_baffle_spacing():
    """Book estimates 68 kPa; full recalculation gives ~39 kPa (L/lB also halves)."""
    wide = replace(geom, baffle_spacing=2 * geom.baffle_spacing)
    s1, s2 = shell_side_kern(geom, methanol), shell_side_kern(wide, methanol)
    dp1, dp2 = shell_side_dp(geom, methanol, s1), shell_side_dp(wide, methanol, s2)
    assert 6 < dp1 / dp2 < 8.5                       # ~1/7, not the book's 1/4
    assert abs(s2["h_s"] / s1["h_s"] - 2**-0.55) < 0.01


def test_jf_rejects_other_baffle_cuts():
    import pytest
    with pytest.raises(ValueError):
        jf_shell_kern(36762, baffle_cut=35)