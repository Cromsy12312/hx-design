"""FCOC design study: sweep -> feasibility -> Plot 2 (mission) + Plot 3 (Pareto).

Run from the repo root:  python -m examples.fcoc_study
"""
from itertools import product

import numpy as np
import matplotlib.pyplot as plt

from hx.plate_fin import OffsetStripFin, PlateFinCore
from hx.fcoc import evaluate_mission, is_feasible, pareto_front, MISSION, FUEL_LIMIT_C

# --- 1. Design space ---
OIL_FINS = [OffsetStripFin(s=1.5e-3, h=6e-3, t=0.15e-3, l=3e-3),
            OffsetStripFin(s=2.0e-3, h=6e-3, t=0.15e-3, l=3e-3)]
FUEL_FINS = [OffsetStripFin(s=1.0e-3, h=2.5e-3, t=0.10e-3, l=3e-3),
             OffsetStripFin(s=1.0e-3, h=3.5e-3, t=0.10e-3, l=3e-3)]
N_HOTS = [1, 2, 3]
LENGTHS = np.arange(0.10, 0.50, 0.05)
MATERIALS = {"aluminium": (170.0, 2700.0), "titanium": (6.7, 4430.0)}


def sweep():
    """Evaluate every design in the space; return feasible records + counts."""
    feasible = []
    tried = {m: 0 for m in MATERIALS}
    n_feasible = {m: 0 for m in MATERIALS}

    for oil_fin, fuel_fin, N, L_hot, L_cold in product(OIL_FINS, FUEL_FINS, N_HOTS,
                                                      LENGTHS, LENGTHS):
        for material, (k_metal, rho_metal) in MATERIALS.items():
            core = PlateFinCore(L_hot=L_hot, L_cold=L_cold, N_hot=N,
                                fin_hot=oil_fin, fin_cold=fuel_fin,
                                k_metal=k_metal, rho_metal=rho_metal)
            tried[material] += 1
            res = evaluate_mission(core)
            if is_feasible(res):
                n_feasible[material] += 1
                dP_fuel = max(res[p]["result"]["dP_cold"] for p in res) / 1e3   # kPa
                dP_oil = max(res[p]["result"]["dP_hot"] for p in res) / 1e3     # kPa
                mass = res["takeoff"]["result"]["mass"]
                feasible.append(dict(material=material, mass=mass, dP_fuel=dP_fuel, dP_oil=dP_oil,
                                     core=core, res=res))
    return feasible, tried, n_feasible


def plot_mission(design):
    """Plot 2: required vs achieved duty, and fuel outlet temperature, per mission point."""
    res = design["res"]
    names = list(MISSION.keys())
    labels = [n.replace("_", "\n") for n in names]
    q_req = [MISSION[n][2] / 1e3 for n in names]
    q_ach = [res[n]["result"]["q"] / 1e3 for n in names]
    T_fuel = [res[n]["result"]["Tc_out"] for n in names]

    x = np.arange(len(names))
    w = 0.38
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))

    ax1.bar(x - w / 2, q_req, w, label="Required")
    ax1.bar(x + w / 2, q_ach, w, label="Achieved")
    for xi, req, ach in zip(x, q_req, q_ach):
        ax1.text(xi + w / 2, ach + 0.8, f"+{100 * (ach / req - 1):.0f}%",
                 ha="center", fontsize=9)
    ax1.set_xticks(x, labels)
    ax1.set_ylabel("Oil heat rejection [kW]")
    ax1.set_title("Duty across the mission")
    ax1.legend()
    ax1.grid(axis="y", alpha=0.3)

    ax2.plot(x, T_fuel, "o-", lw=2)
    ax2.axhline(FUEL_LIMIT_C, ls="--", color="grey")
    ax2.text(x[0], FUEL_LIMIT_C + 2, f"fuel limit {FUEL_LIMIT_C:.0f} °C", color="grey")
    for xi, T in zip(x, T_fuel):
        ax2.annotate(f"{T:.0f} °C", (xi, T), xytext=(0, 8),
                     textcoords="offset points", ha="center", fontsize=9)
    ax2.set_xticks(x, labels)
    ax2.set_ylabel("Fuel outlet temperature [°C]")
    ax2.set_ylim(0, FUEL_LIMIT_C + 20)
    ax2.set_title("Fuel heating across the mission")
    ax2.grid(alpha=0.3)

    c = design["core"]
    fig.suptitle(f"Lightest feasible FCOC: aluminium, {design['mass']:.2f} kg, "
                 f"{c.N_hot} oil layer(s), {c.L_hot:.2f} × {c.L_cold:.2f} m")
    fig.tight_layout()
    fig.savefig("figures/fcoc_mission.png", dpi=200)


def plot_pareto(feasible, design, tried, n_feasible):
    """Plot 3: mass vs oil-side pressure drop, with the Pareto front."""
    al = [d for d in feasible if d["material"] == "aluminium"]
    pts = [(d["mass"], d["dP_oil"]) for d in al]
    front = pareto_front(pts)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter([p[0] for p in pts], [p[1] for p in pts], s=22, alpha=0.3,
               color="grey", label=f"Feasible aluminium designs ({len(pts)})")
    ax.plot([p[0] for p in front], [p[1] for p in front], "o-", lw=2, ms=7,
            label="Pareto front (no lighter design has lower ΔP)")
    for m, dp in front:
        ax.annotate(f"{m:.2f} kg\n{dp:.1f} kPa", (m, dp), xytext=(8, 4),
                    textcoords="offset points", fontsize=8)
    ax.plot(design["mass"], design["dP_oil"], "*", ms=18, zorder=5,
            label=f"Lightest: {design['mass']:.2f} kg")

    ax.text(0.98, 0.60,
            f"Titanium: {n_feasible['titanium']} feasible of {tried['titanium']} tried\n"
            f"(k = 6.7 W/m·K collapses fin efficiency)",
            transform=ax.transAxes, ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round", fc="white", alpha=0.9))

    ax.set_xlabel("Core mass [kg]")
    ax.set_ylabel("Max oil-side pressure drop over mission [kPa]")
    ax.set_ylim(bottom=0)
    ax.set_title("FCOC trade-off: core mass vs oil pressure drop\n"
                 "(all feasible designs meet duty, fuel limit and correlation range)")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig("figures/fcoc_pareto.png", dpi=200)


def main():
    feasible, tried, n_feasible = sweep()
    for m in MATERIALS:
        print(f"{m}: {n_feasible[m]} feasible of {tried[m]} tried")

    al_feasible = [d for d in feasible if d["material"] == "aluminium"]
    if not al_feasible:
        raise RuntimeError("No feasible aluminium designs: widen the design space")
    lightest = min(al_feasible, key=lambda d: d["mass"])

    c = lightest["core"]
    print(f"\nLightest: {lightest['mass']:.3f} kg, N_hot={c.N_hot}, "
          f"L_hot={c.L_hot:.2f} m, L_cold={c.L_cold:.2f} m, "
          f"oil fin s={c.fin_hot.s*1e3:.1f} mm, fuel fin h={c.fin_cold.h*1e3:.1f} mm, "
          f"max fuel dP={lightest['dP_fuel']:.1f} kPa")
    for name, point in lightest["res"].items():
        r = point["result"]
        print(f"  {name:13s} q={r['q']/1e3:5.1f} kW (need {MISSION[name][2]/1e3:.0f}), "
              f"fuel out={r['Tc_out']:5.1f} C, Re oil={r['Re_hot']:.0f}, Re fuel={r['Re_cold']:.0f}")

    plot_mission(lightest)
    plot_pareto(feasible, lightest, tried, n_feasible)
    plt.show()


if __name__ == "__main__":
    main()