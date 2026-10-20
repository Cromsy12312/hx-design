# hx-design

Heat exchanger sizing and rating tool for industrial and aero-engine applications.
*(intro: write Friday)*

## Features
- LMTD and effectiveness-NTU methods (counterflow, parallel)
- Double-pipe sizing and rating, with temperature-dependent properties (CoolProp, Incropera tables)
- *(shell-and-tube, plate-fin: coming)*

## Validation

### Case 1: Incropera Ex 11.1 (double-pipe lube-oil cooler)
| Quantity | Incropera | This tool | Error |
|---|---|---|---|
| q | 8524 W | 8524 W | 0.0% |
| h_i | 2250 | 2249 | 0.05% |
| h_o | 38.8 | 38.93 | +0.3% |
| U | 38.1 | 38.27 | +0.4% |
| L | 65.9 m | 65.65 m | −0.4% |

Notes:
- Gap comes from the book rounding Nu_o to 5.63 (interpolated value 5.64)
- Oil side ≈ 98% of thermal resistance, so Gnielinski vs Dittus–Boelter changes L by <1%
- With temperature-dependent properties: Th_out = 59.91 °C vs book 60 °C (4 iterations)
- Using inlet-temperature properties instead of the mean gives a 1.2 °C error (oil μ halves from 80 to 100 °C)

### Case 2: Sinnott Ex 12.1 (shell-and-tube, Kern)
*(Wednesday)*

## Results
*(plots: Wed/Thu)*

## Limitations
- Transition region (Re 2300–3000) uses Gnielinski without correction
- *(add as you go)*

## Running
pip install -r requirements.txt
python -m pytest
