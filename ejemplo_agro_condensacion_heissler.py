#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cálculo de h (condensación por película en superficie vertical, Nusselt),
número de Biot y tiempo transitorio con Heisler (1 término) para una placa.

Caso de estudio: pouch agroindustrial con puré de tomate en autoclave a 121°C.
"""

from math import pi, sin, tan, log, exp
import numpy as np

# -----------------------------------------------
# (1) Correlación de Nusselt: condensación laminar en placa vertical
# -----------------------------------------------
def h_condensacion_nusselt_vertical(
    T_sat, T_s, rho_l, k_l, mu_l, h_fg, L_cond, rho_v=1.0, g=9.81
):
    """
    Retorna h [W/m^2-K] por Nusselt para condensación laminar en superficie vertical.
    """
    dT = T_sat - T_s
    num = rho_l * (rho_l - rho_v) * g * h_fg * (k_l**3)
    den = mu_l * L_cond * dT
    h = 0.943 * (num / den) ** 0.25
    return h

# -----------------------------------------------
# (2) Propiedades efectivas y difusividad
# -----------------------------------------------
def difusividad_termica(k, rho, cp):
    return k / (rho * cp)

# -----------------------------------------------
# (3) Heisler 1 término
#   - Caso Bi -> ∞ (superficie a T_surf fija): θ/θ_i = (4/π) exp(-π² Fo / 4)
#   - Caso Bi finito (placa): θ/θ_i = A1 * exp(-μ1² Fo), con μ1 tan μ1 = Bi
# -----------------------------------------------
def tiempo_heisler_infinite_bi(T_surf, T_i, T_center_obj, L, alpha):
    """
    Tiempo (s) para que T(centro)=T_center_obj con superficie a T_surf (Bi→∞).
    Usa la forma 1-término: (T_surf - T_c)/(T_surf - T_i) = (4/π) exp(-π² α t / (4 L²))
    """
    R = (T_surf - T_center_obj) / (T_surf - T_i)
    if not (0 < R < 1):
        raise ValueError("El objetivo de temperatura produce R fuera de (0,1). Revisa datos.")
    t = -(4 * L**2) / (pi**2 * alpha) * log(R * pi / 4.0)
    return t

def _mu1_from_bi(Bi):
    """
    Resuelve μ1 tan(μ1) = Bi para μ1 in (0, π/2) usando bisección simple.
    (Sin SciPy para mantener el script auto-contenido)
    """
    # Para Bi>0, la primera raíz está en (0, π/2)
    a, b = 1e-9, pi/2 - 1e-9
    fa = a * tan(a) - Bi
    fb = b * tan(b) - Bi
    # Bisección
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
    """
    Tiempo (s) para placa con convección en la superficie (Bi finito, Heisler 1 término).
    θ/θ_i (centro) = A1 * exp(-μ1² Fo), donde:
      μ1 tan(μ1) = Bi,
      A1 = 4 sin(μ1) / (2 μ1 + sin(2 μ1))
    """
    if Bi <= 0:
        raise ValueError("Bi debe ser > 0 para el caso convectivo.")
    R = (T_inf - T_center_obj) / (T_inf - T_i)
    if not (0 < R < 1):
        raise ValueError("El objetivo de temperatura produce R fuera de (0,1). Revisa datos.")
    mu1 = _mu1_from_bi(Bi)
    A1 = 4.0 * sin(mu1) / (2.0 * mu1 + sin(2.0 * mu1))
    if R >= A1:  # si R>=A1, el 1er término no puede alcanzar ese R (se necesitarían más términos)
        raise ValueError(
            f"Para Bi={Bi:.2f}, R={R:.4f} ≥ A1={A1:.4f}. "
            "El objetivo es demasiado exigente para 1 término; usa la serie completa o revisa datos."
        )
    Fo = - (1.0 / mu1**2) * log(R / A1)
    t = Fo * L**2 / alpha
    return t, mu1, A1, Fo

# -----------------------------------------------
# (4) Centro T(t) para graficar (infinite Bi y finite Bi)
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
    T_sat = 121.0  # °C (vapor saturado)
    T_s_aprox = 100.0  # °C (superficie del pouch ~ isoterma por Bi alto)
    T_i = 25.0     # °C (inicial del alimento)
    T_c_obj = 100.0  # °C (objetivo en el centro)

    # Geometría
    L_cond = 0.25     # m (altura característica para condensación)
    L = 0.015         # m (semiespesor de la placa; 2L = 0.03 m)

    # Propiedades del condensado (agua) a T de película ~ (T_sat + T_s)/2
    rho_l = 958.0     # kg/m3
    rho_v = 1.0       # kg/m3
    k_l = 0.68        # W/m-K
    mu_l = 2.8e-4     # Pa·s
    h_fg = 2.2e6      # J/kg
    g = 9.81          # m/s2

    # Propiedades del alimento
    k_s = 0.60        # W/m-K
    rho_s = 1050.0    # kg/m3
    cp_s = 4000.0     # J/kg-K
    alpha = difusividad_termica(k_s, rho_s, cp_s)

    # (1) h por Nusselt (condensación)
    h = h_condensacion_nusselt_vertical(
        T_sat=T_sat, T_s=T_s_aprox, rho_l=rho_l, k_l=k_l, mu_l=mu_l,
        h_fg=h_fg, L_cond=L_cond, rho_v=rho_v, g=g
    )

    # (2) Biot
    Bi = h * L / k_s

    # (3) Tiempo con Heisler 1 término
    t_inf = tiempo_heisler_infinite_bi(
        T_surf=T_sat, T_i=T_i, T_center_obj=T_c_obj, L=L, alpha=alpha
    )
    Fo_inf = alpha * t_inf / (L**2)

    # (3b) (Opcional) Tiempo con Bi finito (placa convectiva)
    #      Útil si quieres comparar con el caso exacto convectivo 1-término
    try:
        t_bi, mu1, A1, Fo_bi = tiempo_heisler_finite_bi(
            T_inf=T_sat, T_i=T_i, T_center_obj=T_c_obj, L=L, alpha=alpha, Bi=Bi
        )
        usado_bi_finito = True
    except Exception as e:
        t_bi = None
        usado_bi_finito = False
        msg_bi = str(e)

    # --------- Salida de resultados ----------
    print("\n=== RESULTADOS ===")
    print(f"h (Nusselt, condensación)        = {h:8.2f} W/m^2-K")
    print(f"Biot (h L / k_s)                 = {Bi:8.2f}  (-)")
    print(f"alpha (difusividad)              = {alpha:.3e} m^2/s")

    print("\n-- Heisler (1 término), superficie a T_surf = T_sat (Bi→∞) --")
    print(f"t_infinite_Bi                    = {t_inf:8.2f} s  ({t_inf/60:5.2f} min)")
    print(f"Fo (infinite Bi)                 = {Fo_inf:8.3f}  (válido si Fo ≳ 0.2)")

    if usado_bi_finito:
        print("\n-- Heisler (1 término) con Bi finito (placa convectiva) --")
        print(f"μ1 (raíz μ1 tan μ1 = Bi)         = {mu1:8.5f} rad")
        print(f"A1                               = {A1:8.5f} (-)")
        print(f"t_finite_Bi                      = {t_bi:8.2f} s  ({t_bi/60:5.2f} min)")
        print(f"Fo (finite Bi)                   = {Fo_bi:8.3f}")
        print("(Con Bi ≈ alto, la diferencia vs. Bi→∞ es pequeña.)")
    else:
        print("\n[Aviso] No se calculó el caso Bi finito (1 término):")
        print("       ", msg_bi)

    # --------- (Opcional) Curva T_centro(t) para visualizar ---------
    try:
        import matplotlib.pyplot as plt

        t_max = t_inf * 1.3
        t_vec = np.linspace(0.0, t_max, 200)

        Tc_inf_vec = [T_centro_infinite_bi(t, T_surf=T_sat, T_i=T_i, L=L, alpha=alpha)
                      for t in t_vec]

        plt.figure()
        plt.plot(t_vec/60.0, Tc_inf_vec, label="Heisler 1T (Bi→∞)")
        if usado_bi_finito:
            Tc_bi_vec = [T_centro_finite_bi(t, T_inf=T_sat, T_i=T_i, L=L, alpha=alpha, Bi=Bi)
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
        print("\n[Nota] No se generó gráfica (matplotlib no disponible).")
        print("      Error:", e)
