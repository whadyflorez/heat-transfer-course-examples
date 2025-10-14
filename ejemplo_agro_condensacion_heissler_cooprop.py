#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cálculo con CoolProp de:
  - h por condensación laminar en superficie vertical (Nusselt)
  - Número de Biot
  - Tiempo transitorio por Heisler (1 término) para placa (Bi→∞ y opcional Bi finito)

Caso: pouch agroindustrial con puré de tomate; autoclave con vapor saturado a 121°C.
"""

from math import pi, sin, tan, log, exp
import numpy as np

# --- CoolProp ---
from CoolProp.CoolProp import PropsSI

# -----------------------------------------------
# (1) Correlación de Nusselt: condensación laminar en placa vertical
#     Propiedades:
#       - ρ_l, μ_l, k_l evaluadas a T_film = (T_sat + T_s)/2
#       - ρ_v y h_fg a saturación a T_sat
# -----------------------------------------------
def h_condensacion_nusselt_vertical_coolprop(
    T_sat_C, T_s_C, L_cond, fluid="Water", g=9.81
):
    """
    Retorna h [W/m^2-K] por Nusselt para condensación laminar en superficie vertical,
    evaluando propiedades con CoolProp.

    Convención usual:
      - Propiedades del condensado (ρ_l, μ_l, k_l) a T_film = (T_sat + T_s)/2
      - ρ_v y h_fg a saturación a T_sat
    """
    T_sat = T_sat_C + 273.15
    T_s = T_s_C + 273.15
    T_film = 0.5 * (T_sat + T_s)

    # Propiedades en la película (líquido saturado a T_film)
    rho_l = PropsSI('D', 'T', T_film, 'Q', 0, fluid)              # kg/m^3
    mu_l  = PropsSI('VISCOSITY', 'T', T_film, 'Q', 0, fluid)      # Pa·s
    k_l   = PropsSI('CONDUCTIVITY', 'T', T_film, 'Q', 0, fluid)   # W/m-K

    # Propiedades a saturación a T_sat
    rho_v = PropsSI('D', 'T', T_sat, 'Q', 1, fluid)               # kg/m^3
    h_g   = PropsSI('H', 'T', T_sat, 'Q', 1, fluid)               # J/kg
    h_l   = PropsSI('H', 'T', T_sat, 'Q', 0, fluid)               # J/kg
    h_fg  = h_g - h_l                                             # J/kg

    dT = T_sat - T_s  # K

    # Nusselt vertical laminar:
    # h = 0.943 * [ rho_l (rho_l - rho_v) g h_fg k_l^3 / ( mu_l L_cond ΔT ) ]^(1/4)
    num = rho_l * (rho_l - rho_v) * g * h_fg * (k_l**3)
    den = mu_l * L_cond * dT
    h = 0.943 * (num / den) ** 0.25

    props = {
        "T_film_K": T_film, "rho_l": rho_l, "mu_l": mu_l, "k_l": k_l,
        "rho_v": rho_v, "h_fg": h_fg, "DeltaT_K": dT
    }
    return h, props

# -----------------------------------------------
# (2) Difusividad térmica del alimento
# -----------------------------------------------
def difusividad_termica(k, rho, cp):
    return k / (rho * cp)

# -----------------------------------------------
# (3) Heisler 1 término
#   - Bi → ∞ (superficie a T_surf fija): θ/θ_i = (4/π) exp(-π² Fo / 4)
#   - Bi finito (placa): θ/θ_i = A1 exp(-μ1² Fo), con μ1 tan μ1 = Bi
# -----------------------------------------------
def tiempo_heisler_infinite_bi(T_surf, T_i, T_center_obj, L, alpha):
    R = (T_surf - T_center_obj) / (T_surf - T_i)
    if not (0 < R < 1):
        raise ValueError("El objetivo de temperatura produce R fuera de (0,1).")
    t = -(4 * L**2) / (pi**2 * alpha) * log(R * pi / 4.0)
    return t

def _mu1_from_bi(Bi):
    # Bisección en (0, π/2) para resolver μ1 tan μ1 = Bi
    a, b = 1e-12, pi/2 - 1e-12
    fa = a * tan(a) - Bi
    fb = b * tan(b) - Bi
    for _ in range(200):
        m = 0.5 * (a + b)
        fm = m * tan(m) - Bi
        if fa * fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
        if abs(b - a) < 1e-12:
            break
    return 0.5 * (a + b)

def tiempo_heisler_finite_bi(T_inf, T_i, T_center_obj, L, alpha, Bi):
    if Bi <= 0:
        raise ValueError("Bi debe ser > 0.")
    R = (T_inf - T_center_obj) / (T_inf - T_i)
    if not (0 < R < 1):
        raise ValueError("El objetivo de temperatura produce R fuera de (0,1).")
    mu1 = _mu1_from_bi(Bi)
    A1 = 4.0 * sin(mu1) / (2.0 * mu1 + sin(2.0 * mu1))
    if R >= A1:
        raise ValueError(
            f"Para Bi={Bi:.2f}, R={R:.4f} ≥ A1={A1:.4f}. "
            "El objetivo es demasiado exigente para 1 término."
        )
    Fo = - (1.0 / mu1**2) * log(R / A1)
    t = Fo * L**2 / alpha
    return t, mu1, A1, Fo

# -----------------------------------------------
# (4) T_centro(t) para graficar (opcional)
# -----------------------------------------------
def T_centro_infinite_bi(t, T_surf, T_i, L, alpha):
    theta_ratio = (4.0 / pi) * exp(-(pi**2) * alpha * t / (4.0 * L**2))
    return T_surf - theta_ratio * (T_surf - T_i)

def T_centro_finite_bi(t, T_inf, T_i, L, alpha, Bi):
    mu1 = _mu1_from_bi(Bi)
    A1 = 4.0 * sin(mu1) / (2.0 * mu1 + sin(2.0 * mu1))
    Fo = alpha * t / (L**2)
    theta_ratio = A1 * exp(-mu1**2 * Fo)
    return T_inf - theta_ratio * (T_inf - T_i)

# ===============================================
# Ejemplo numérico (datos del problema)
# ===============================================
if __name__ == "__main__":
    # Condiciones térmicas
    T_sat_C = 121.0      # °C (vapor saturado)
    T_s_C   = 100.0      # °C (aprox. superficie del pouch)
    T_i     = 25.0       # °C (inicial del alimento)
    T_c_obj = 100.0      # °C (objetivo en el centro)

    # Geometría
    L_cond = 0.25        # m (altura característica para condensación)
    L      = 0.015       # m (semiespesor; 2L = 0.03 m)

    # Propiedades del alimento (puré tipo agua espesa)
    k_s   = 0.60         # W/m-K
    rho_s = 1050.0       # kg/m^3
    cp_s  = 4000.0       # J/kg-K
    alpha = difusividad_termica(k_s, rho_s, cp_s)

    # (1) h por Nusselt (CoolProp)
    h, props = h_condensacion_nusselt_vertical_coolprop(
        T_sat_C=T_sat_C, T_s_C=T_s_C, L_cond=L_cond, fluid="Water"
    )

    # (2) Biot
    Bi = h * L / k_s

    # (3) Heisler (Bi→∞)
    t_inf = tiempo_heisler_infinite_bi(
        T_surf=T_sat_C, T_i=T_i, T_center_obj=T_c_obj, L=L, alpha=alpha
    )
    Fo_inf = alpha * t_inf / (L**2)

    # (3b) Heisler con Bi finito (1 término, opcional)
    try:
        t_bi, mu1, A1, Fo_bi = tiempo_heisler_finite_bi(
            T_inf=T_sat_C, T_i=T_i, T_center_obj=T_c_obj, L=L, alpha=alpha, Bi=Bi
        )
        usado_bi_finito = True
    except Exception as e:
        usado_bi_finito = False
        msg_bi = str(e)

    # --------- Salidas ----------
    print("\n=== PROPIEDADES (CoolProp) ===")
    print(f"T_film        = {props['T_film_K']-273.15:8.2f} °C")
    print(f"rho_l (film)  = {props['rho_l']:8.2f} kg/m^3")
    print(f"mu_l (film)   = {props['mu_l']:8.3e} Pa·s")
    print(f"k_l  (film)   = {props['k_l']:8.3f} W/m-K")
    print(f"rho_v (Tsat)  = {props['rho_v']:8.5f} kg/m^3")
    print(f"h_fg (Tsat)   = {props['h_fg']/1e6:8.3f} MJ/kg")
    print(f"ΔT            = {props['DeltaT_K']:8.2f} K")

    print("\n=== RESULTADOS ===")
    print(f"h (Nusselt)   = {h:8.2f} W/m^2-K")
    print(f"Biot          = {Bi:8.2f} (-)")
    print(f"alpha         = {alpha:8.3e} m^2/s")

    print("\n-- Heisler (1 término), Bi→∞ (superficie a T_surf = T_sat) --")
    print(f"t_infinite_Bi = {t_inf:8.2f} s  ({t_inf/60:5.2f} min)")
    print(f"Fo (∞ Bi)     = {Fo_inf:8.3f}  (válido si Fo ≳ 0.2)")

    if usado_bi_finito:
        print("\n-- Heisler (1 término) con Bi finito (placa convectiva) --")
        print(f"μ1            = {mu1:8.5f} rad")
        print(f"A1            = {A1:8.5f} (-)")
        print(f"t_finite_Bi   = {t_bi:8.2f} s  ({t_bi/60:5.2f} min)")
        print(f"Fo (Bi finito)= {Fo_bi:8.3f}")
    else:
        print("\n[Aviso] No se calculó el caso Bi finito (1 término):")
        print("       ", msg_bi)

    # --------- (Opcional) Gráfica T_centro(t) ---------
    try:
        import matplotlib.pyplot as plt

        t_max = t_inf * 1.3
        t_vec = np.linspace(0.0, t_max, 200)

        Tc_inf_vec = [T_centro_infinite_bi(t, T_surf=T_sat_C, T_i=T_i, L=L, alpha=alpha)
                      for t in t_vec]

        plt.figure()
        plt.plot(t_vec/60.0, Tc_inf_vec, label="Heisler 1T (Bi→∞)")
        if usado_bi_finito:
            Tc_bi_vec = [T_centro_finite_bi(t, T_inf=T_sat_C, T_i=T_i, L=L, alpha=alpha, Bi=Bi)
                         for t in t_vec]
            plt.plot(t_vec/60.0, Tc_bi_vec, linestyle="--", label=f"Heisler 1T (Bi={Bi:.0f})")

        plt.axhline(T_c_obj, linestyle=":", label=f"Objetivo centro = {T_c_obj} °C")
        plt.xlabel("Tiempo [min]")
        plt.ylabel("T centro [°C]")
        plt.title("Calentamiento transitorio del centro (placa)")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()
    except Exception as e:
        print("\n[Nota] No se generó la gráfica (matplotlib no disponible).")
        print("      Error:", e)
