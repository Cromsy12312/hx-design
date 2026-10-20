import numpy as np
import pytest
from hx.methods import heat_capacity_rate, lmtd, effectiveness, ntu_from_effectiveness, effectiveness_crossflow_unmixed, ntu_crossflow_unmixed

# Incropera Ex 11.1 data
m_h, cp_h = 0.1, 2131      # oil
m_c, cp_c = 0.2, 4178      # water
Th_in, Th_out, Tc_in = 100.0, 60.0, 30.0


def test_energy_balance():
    q = heat_capacity_rate(m_h, cp_h) * (Th_in - Th_out)
    assert abs(q - 8524) < 1
    Tc_out = Tc_in + q / heat_capacity_rate(m_c, cp_c)
    assert abs(Tc_out - 40.2) < 0.1


def test_lmtd_matches_incropera():
    assert abs(lmtd(Th_in, Th_out, Tc_in, 40.2) - 43.2) < 0.1


def test_lmtd_and_ntu_give_same_area():
    """Core consistency check: both methods must size the same exchanger."""
    U = 38.1
    Ch, Cc = heat_capacity_rate(m_h, cp_h), heat_capacity_rate(m_c, cp_c)
    Cmin, Cmax = min(Ch, Cc), max(Ch, Cc)
    q = Ch * (Th_in - Th_out)
    Tc_out = Tc_in + q / Cc

    A_lmtd = q / (U * lmtd(Th_in, Th_out, Tc_in, Tc_out))

    eps = q / (Cmin * (Th_in - Tc_in))
    A_ntu = ntu_from_effectiveness(eps, Cmin / Cmax) * Cmin / U

    assert np.isclose(A_lmtd, A_ntu, rtol=1e-6)


def test_effectiveness_roundtrip():
    for arr in ["counterflow", "parallel"]:
        eps = effectiveness(1.5, 0.5, arr)
        assert np.isclose(ntu_from_effectiveness(eps, 0.5, arr), 1.5)


def test_balanced_counterflow_edge_case():
    assert np.isclose(lmtd(80, 50, 20, 50), 30)
    assert np.isclose(effectiveness(2.0, 1.0), 2 / 3)


from hx.methods import f_factor_1_2, effectiveness_shell_1_2, ntu_shell_1_2


def test_f_factor_sinnott_12_1():
    """Sinnott reads F = 0.85 from a chart; the exact formula gives 0.812."""
    R = (95 - 40) / (40 - 25)
    S = (40 - 25) / (95 - 25)
    assert abs(f_factor_1_2(R, S) - 0.812) < 0.002


def test_f_factor_limits():
    assert np.isclose(f_factor_1_2(2.0, 1e-6), 1.0, atol=1e-3)   # tiny S -> F ~ 1
    assert 0 < f_factor_1_2(1.0, 0.5) < 1                          # R = 1 branch works


def test_shell_lmtd_and_ntu_agree():
    """F*LMTD sizing and shell eps-NTU sizing must give the same UA."""
    Ch, Cc = 1.0, 55 / 15            # Sinnott 12.1 ratio: hot changes 55 K, cold 15 K
    q = Ch * (95 - 40)
    R, S = (95 - 40) / (40 - 25), (40 - 25) / (95 - 25)
    UA_lmtd = q / (f_factor_1_2(R, S) * lmtd(95, 40, 25, 40))
    eps = q / (Ch * (95 - 25))
    UA_ntu = ntu_shell_1_2(eps, Ch / Cc) * Ch
    assert np.isclose(UA_lmtd, UA_ntu, rtol=1e-9)


def test_shell_roundtrip():
    eps = effectiveness_shell_1_2(1.2, 0.6)
    assert np.isclose(ntu_shell_1_2(eps, 0.6), 1.2)


def test_crossflow_incropera_11_3_sizing():
    """Book reads NTU ~ 2.0 from Fig 11.14 -> A = 37.8 m2; exact equation gives ~38.3 m2."""
    Cmin, Cmax, U = 1889, 4197, 100
    eps = (4197 * (125 - 35)) / (Cmin * (300 - 35))      # = 0.755
    NTU = ntu_crossflow_unmixed(eps, Cmin / Cmax)
    A = NTU * Cmin / U
    assert abs(A - 37.8) / 37.8 < 0.02


def test_crossflow_incropera_11_4_rating():
    """Book reads eps = 0.82 from the chart; exact equation gives ~0.845."""
    NTU = 100 * 40 / 1500
    eps = effectiveness_crossflow_unmixed(NTU, 1500 / 4197)
    assert abs(eps - 0.82) / 0.82 < 0.04


def test_crossflow_roundtrip():
    for NTU, Cr in [(0.5, 0.2), (1.5, 0.6), (3.0, 1.0)]:
        eps = effectiveness_crossflow_unmixed(NTU, Cr)
        assert np.isclose(ntu_crossflow_unmixed(eps, Cr), NTU, rtol=1e-6)


def test_crossflow_cr_zero_limit():
    """With one fluid of infinite capacity, every arrangement gives eps = 1 - exp(-NTU)."""
    assert np.isclose(effectiveness_crossflow_unmixed(2.0, 0.0), 1 - np.exp(-2.0))


def test_crossflow_unreachable_eps_raises():
    with pytest.raises(ValueError):
        ntu_crossflow_unmixed(1.0, 0.5)