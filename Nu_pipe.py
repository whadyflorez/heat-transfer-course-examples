import CoolProp.CoolProp as CP
import numpy as np

# Propiedades del agua a 50°C
T_inf = 323.15  # K (50°C)
P = 101325  # Pa (presión atmosférica)

rho = CP.PropsSI('D', 'T', T_inf, 'P', P, 'Water')  # Densidad (kg/m³)
mu = CP.PropsSI('V', 'T', T_inf, 'P', P, 'Water')  # Viscosidad dinámica (Pa.s)
k = CP.PropsSI('L', 'T', T_inf, 'P', P, 'Water')  # Conductividad térmica (W/m.K)
Pr = CP.PropsSI('PRANDTL', 'T', T_inf, 'P', P, 'Water')  # Número de Prandtl

# Condiciones del problema
U = 0.12  # m/s (velocidad del agua)
D = 0.02  # m (diámetro del tubo)
L = 5.0  # m (longitud del tubo)
mu_w = mu  # Se asume temperatura de pared similar a la del fluido

# Cálculo del número de Reynolds
Re = (rho * U * D) / mu

# Cálculo del número de Nusselt
if Re < 2300:  # Flujo laminar
    Nu = 1.86 * (Re * Pr * (D / L))**(1/3) * (mu / mu_w)**0.14
elif Re > 4000:  # Flujo turbulento (Dittus-Boelter)
    Nu = 0.023 * Re**0.8 * Pr**0.3
else:
    Nu = np.nan  # Zona de transición, no se aplica correlación directa

# Coeficiente de convección
h = (Nu * k) / D

# Resultados
print(f"Número de Reynolds: {Re:.2f}")
print(f"Número de Nusselt: {Nu:.2f}")
print(f"Coeficiente de convección h: {h:.2f} W/m²K")