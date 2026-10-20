"""Plot 1: effect of fouling on required area (Sinnott Ex 12.1 duty)."""
import numpy as np
import matplotlib.pyplot as plt

from hx.double_pipe import Stream
from hx.methods import lmtd, f_factor_1_2
from hx.shell_tube import (ShellTubeGeometry, tube_side, shell_side_kern,
                           overall_U_shell_tube)

# --- Sinnott Ex 12.1 case ---
methanol = Stream(m_dot=100000 / 3600, cp=2840, mu=0.34e-3, k=0.19, T_in=95, rho=750)
water = Stream(m_dot=68.9, cp=4200, mu=0.8e-3, k=0.59, T_in=25, rho=995)
geom = ShellTubeGeometry(do=0.020, di=0.016, L=4.83, Nt=918, n_passes=2,
                         Ds=0.894, baffle_spacing=0.178)
Th_in, Th_out, Tc_in, Tc_out = 95, 40, 25, 40

# Duty and corrected mean temperature difference (fixed by the process)
q = methanol.m_dot * methanol.cp * (Th_in - Th_out)
R = (Th_in - Th_out) / (Tc_out - Tc_in)
S = (Tc_out - Tc_in) / (Th_in - Tc_in)
dT_m = f_factor_1_2(R, S) * lmtd(Th_in, Th_out, Tc_in, Tc_out)

# Film coefficients (don't depend on fouling)
h_i = tube_side(geom, water, heating=True)["h_i"]
h_s = shell_side_kern(geom, methanol)["h_s"]

# --- Sweep: scale Sinnott's fouling resistances from 0 (clean) to 2x ---
R_fo_book, R_fi_book = 1 / 5000, 1 / 3000
scale = np.linspace(0, 2, 81)
A_req = []
for s in scale:
    U = overall_U_shell_tube(h_s, h_i, geom, R_fo=s * R_fo_book, R_fi=s * R_fi_book)
    A_req.append(q / (U * dT_m))
A_req = np.array(A_req)
A_clean = A_req[0]

# Total fouling resistance referred to the tube outside area
R_total = scale * (R_fo_book + (geom.do / geom.di) * R_fi_book)

# --- Plot ---
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(R_total * 1e4, A_req, lw=2)
i_book = np.argmin(abs(scale - 1))
ax.plot(R_total[i_book] * 1e4, A_req[i_book], "o", ms=8)
ax.annotate(f"Sinnott design fouling\n{A_req[i_book]:.0f} m² (+{100*(A_req[i_book]/A_clean-1):.0f}% vs clean)",
            (R_total[i_book] * 1e4, A_req[i_book]),
            xytext=(10, -40), textcoords="offset points")
ax.axhline(A_clean, ls="--", lw=1, color="grey")
ax.text(0.2, A_clean + 4, f"clean: {A_clean:.0f} m²", color="grey")

ax.set_xlabel("Total fouling resistance  [×10⁻⁴ m²K/W, outside-area basis]")
ax.set_ylabel("Required heat transfer area  [m²]")
ax.set_title("Fouling allowance vs required area\nShell-and-tube methanol cooler (Sinnott Ex 12.1 duty)")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig("figures/fouling_vs_area.png", dpi=200)
print(f"q = {q/1e3:.0f} kW, dT_m = {dT_m:.1f} K")
print(f"Clean area: {A_clean:.0f} m2, design area: {A_req[i_book]:.0f} m2")
plt.show()