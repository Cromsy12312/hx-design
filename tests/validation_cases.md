# Validation cases

## Case 1 — Double-pipe counterflow (Incropera 7th ed., Example 11.1, p. 716–718)

### Inputs
- Hot fluid: unused engine oil in annulus, m_dot = 0.1 kg/s, T_in = 100 °C, T_out = 60 °C (required)
- Cold fluid: water in inner tube, m_dot = 0.2 kg/s, T_in = 30 °C
- Inner tube D_i = 0.025 m (thin wall), annulus outer D_o = 0.045 m
- Arrangement: counterflow
- Assumptions: no wall resistance, no fouling, fully developed flow, constant properties
- Book properties:
  - Oil @ 80 °C (Table A.5): cp = 2131 J/kg·K, mu = 3.25e-2 N·s/m², k = 0.138 W/m·K
  - Water @ 35 °C (Table A.6): cp = 4178 J/kg·K, mu = 725e-6 N·s/m², k = 0.625 W/m·K, Pr = 4.85

### Published results
- q = 8524 W
- T_c,out = 40.2 °C
- ΔT_lm = 43.2 °C
- Water: Re = 14,050 (turbulent), Nu = 90 (Dittus–Boelter), h_i = 2250 W/m²K
- Oil: D_h = 0.020 m, Re = 56.0 (laminar), Nu_i = 5.63 (Table 8.2, D_i/D_o = 0.56), h_o = 38.8 W/m²K
- U = 38.1 W/m²K
- Required length L = 65.9 m

## Case 2 — Shell-and-tube, Kern method (Sinnott, Coulson & Richardson Vol. 6, 4th ed., Example 12.1, pp. 675–679)

### Inputs
- Shell side (hot): methanol, 100,000 kg/h (27.78 kg/s), 95 °C → 40 °C
- Tube side (cold): brackish water, 25 °C → 40 °C (corrosive, so tube side)
- 1 shell pass, 2 tube passes
- Tubes: 20 mm OD, 16 mm ID, effective length 4.83 m, cupro-nickel (k = 50 W/m·K)
- Triangular pitch 1.25 × OD = 25 mm
- Split-ring floating head, bundle clearance 68 mm
- Baffle cut 25%; baffle spacing = Ds/5 = 178 mm (initial design)
- Fouling: methanol 5000 W/m²K, brackish water 3000 W/m²K (as coefficients)
- Book properties:
  - Methanol @ 68 °C: rho = 750, mu = 0.34e-3, cp = 2840, k = 0.19, mu_wall = 0.37e-3
  - Water @ 33 °C: rho = 995, mu = 0.8e-3, cp = 4200, k = 0.59

### Published results
- Q = 4340 kW, water flow = 68.9 kg/s
- ΔT_lm = 31 °C, R = 3.67, S = 0.21, Ft = 0.85, ΔT_m = 26 °C
- Assumed U = 600 → A = 278 m², 918 tubes, Db = 826 mm, Ds = 894 mm
- Tube side: 459 tubes/pass, flow area 0.092 m², u = 0.75 m/s, Re = 14,925, Pr = 5.7
  - h_i = 3852 W/m²K (eq. 12.17, water-specific) / 3812 W/m²K (jh method, used)
- Shell side: As = 0.032 m², Gs = 868 kg/s·m², de = 14.4 mm, Re = 36,762, Pr = 5.1
  - jh = 3.3e-3 → h_s = 2740 W/m²K, viscosity correction 0.99 (negligible)
- Overall: U_o = 738 W/m²K (above the assumed 600, so design OK)
- Pressure drop:
  - Tube side: jf = 4.3e-3 → ΔP_t = 7.2 kPa
  - Shell side: jf = 4e-2 → ΔP_s = 272 kPa (too high)
  - Redesign with baffle spacing doubled: ΔP_s ≈ 68 kPa, h_s = 1573 W/m²K, U = 615 W/m²K