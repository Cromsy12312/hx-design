from hx.double_pipe import Stream, size_double_pipe, rate_double_pipe

# Incropera Ex 11.1, book property values
oil = Stream(m_dot=0.1, cp=2131, mu=3.25e-2, k=0.138, T_in=100)
water = Stream(m_dot=0.2, cp=4178, mu=725e-6, k=0.625, T_in=30)


def run(method):
    return size_double_pipe(oil, water, Th_out=60, Di=0.025, Do=0.045,
                            hot_in="annulus", method=method)


def test_incropera_11_1_intermediates():
    r = run("dittus_boelter")
    assert abs(r["q"] - 8524) < 1
    assert abs(r["Tc_out"] - 40.2) < 0.05
    assert abs(r["LMTD"] - 43.2) < 0.05
    assert abs(r["Re_o"] - 56.0) < 0.5
    assert abs(r["h_i"] - 2250) / 2250 < 0.01
    assert abs(r["h_o"] - 38.8) / 38.8 < 0.01


def test_incropera_11_1_length():
    """Book answer: L = 65.9 m."""
    r = run("dittus_boelter")
    assert abs(r["L"] - 65.9) / 65.9 < 0.01


def test_gnielinski_barely_changes_length():
    """The oil side controls U, so the water-side correlation hardly matters."""
    L_db = run("dittus_boelter")["L"]
    L_gn = run("gnielinski")["L"]
    assert abs(L_gn - L_db) / L_db < 0.01


def test_rating_inverts_sizing():
    """Rate the exchanger we just sized -> must recover Th_out = 60 C."""
    L = run("dittus_boelter")["L"]
    r = rate_double_pipe(oil, water, L=L, Di=0.025, Do=0.045,
                         hot_in="annulus", method="dittus_boelter")
    assert abs(r["Th_out"] - 60.0) < 1e-6
    assert abs(r["Tc_out"] - 40.2) < 0.01


def test_longer_exchanger_cools_more():
    r1 = rate_double_pipe(oil, water, L=40, Di=0.025, Do=0.045)
    r2 = rate_double_pipe(oil, water, L=80, Di=0.025, Do=0.045)
    assert r2["Th_out"] < r1["Th_out"]