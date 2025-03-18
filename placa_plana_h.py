import CoolProp.CoolProp as CP
import numpy as np

# Propiedades del aire a 25°C
T_inf = 298.15  # K (25°C)
P = 101325  # Pa (presión atmosférica)

rho = CP.PropsSI('D', 'T', T_inf, 'P', P, 'Air')  # Densidad (kg/m³)
mu = CP.PropsSI('V', 'T', T_inf, 'P', P, 'Air')  # Viscosidad dinámica (Pa.s)
nu = mu / rho  # Viscosidad cinemática (m²/s)
k = CP.PropsSI('L', 'T', T_inf, 'P', P, 'Air')  # Conductividad térmica (W/m.K)
Pr = CP.PropsSI('PRANDTL', 'T', T_inf, 'P', P, 'Air')  # Número de Prandtl

# Condiciones del problema
U_inf = 10  # m/s (velocidad del aire)
L = 1  # m (longitud de la placa)

# Cálculo del número de Reynolds
Re_x = (rho * U_inf * L) / mu

# Cálculo del número de Nusselt según régimen de flujo
if Re_x < 5e5:
    Nu_x = 0.332 * Re_x**0.5 * Pr**(1/3)  # Flujo laminar
else:
    Nu_x = 0.0296 * Re_x**(4/5) * Pr**(1/3)  # Flujo turbulento

# Coeficiente de convección
h_x = (Nu_x * k) / L

# Resultados
print(f"Número de Reynolds: {Re_x:.2f}")
print(f"Número de Nusselt: {Nu_x:.2f}")
print(f"Coeficiente de convección h: {h_x:.2f} W/m²K")
