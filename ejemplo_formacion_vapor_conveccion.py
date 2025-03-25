import CoolProp.CoolProp as CP
import numpy as np

def boiling_heat_flux(T_s, fluid, surface_material, A):
    # Propiedades a la temperatura de saturación
    T_sat = CP.PropsSI('T','P',101325,'Q',0,fluid)
    rho_l = CP.PropsSI('D','T',T_sat,'Q',0,fluid)
    rho_v = CP.PropsSI('D','T',T_sat,'Q',1,fluid)
    mu_l = CP.PropsSI('V','T',T_sat,'Q',0,fluid)
    h_fg = CP.PropsSI('H','T',T_sat,'Q',1,fluid) - CP.PropsSI('H','T',T_sat,'Q',0,fluid)
    cp_l = CP.PropsSI('C','T',T_sat,'Q',0,fluid)
    Pr_l = CP.PropsSI('PRANDTL','T',T_sat,'Q',0,fluid)
    sigma = CP.PropsSI('SURFACE_TENSION','T',T_sat,'Q',0,fluid)
    g = 9.81  # m/s^2
    
    # Coeficiente empírico para la combinación fluido-superficie (suposición)
    C_sf = 0.013  # Depende del material de la superficie
    n = 0.33
    
    # Exceso de temperatura
    Delta_Te = T_s - T_sat
    
    # Cálculo del flujo de calor
    q_s = mu_l * h_fg * ( (g * (rho_l - rho_v) / sigma )**0.5 ) * ((cp_l * Delta_Te) / (C_sf * h_fg * Pr_l**n))**3
    
    # Cálculo de la cantidad de vapor formado
    m_vapor = (q_s * A) / h_fg  # kg/s
    
    return q_s, m_vapor

# Datos del problema
T_s = 380  # Temperatura de la superficie en Kelvin
fluid = 'Water'
surface_material = 'Copper'  # Material de la superficie
A = 0.1  # Área de la superficie en m²

q_s, m_vapor = boiling_heat_flux(T_s, fluid, surface_material, A)
print(f"Flujo de calor por ebullición nucleada: {q_s:.2f} W/m²")
print(f"Masa de vapor formada: {m_vapor:.6f} kg/s")
