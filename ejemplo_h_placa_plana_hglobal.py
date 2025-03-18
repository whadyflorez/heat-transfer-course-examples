import CoolProp.CoolProp as CP
import numpy as np

# Datos de entrada
T_s = 323.15  # K (temperatura de la superficie de la placa, 50°C)
T_inf = 298.15  # K (25°C)
P = 201325  # Pa (presión atmosférica)
U_inf = 10  # m/s (velocidad del aire)
L = 1  # m (longitud de la placa)

# Temperatura de película
T_film = (T_s + T_inf) / 2  # Promedio de temperatura de pared y fluido

# Propiedades del aire evaluadas a la temperatura de película
rho = CP.PropsSI('D', 'T', T_film, 'P', P, 'Air')  # Densidad (kg/m³)
mu = CP.PropsSI('V', 'T', T_film, 'P', P, 'Air')  # Viscosidad dinámica (Pa.s)
nu = mu / rho  # Viscosidad cinemática (m²/s)
k = CP.PropsSI('L', 'T', T_film, 'P', P, 'Air')  # Conductividad térmica (W/m.K)
Pr = CP.PropsSI('PRANDTL', 'T', T_film, 'P', P, 'Air')

# Cálculo del número de Reynolds
Re_L = (rho * U_inf * L) / mu

# Cálculo del número de Nusselt promedio según régimen de flujo
if Re_L < 5e5:
    Nu_avg = 0.664 * Re_L**0.5 * Pr**(1/3)  # Flujo laminar (promedio)
elif Re_L > 5e5:
    Nu_avg = 0.037 * Re_L**(4/5) * Pr**(1/3)  # Flujo turbulento (promedio)
else:
    Nu_avg = np.nan  # Zona de transición

# Coeficiente de convección promedio sobre la placa
h_avg = (Nu_avg * k) / L

# Resultados
print(f"Temperatura de película: {T_film:.2f} K")
print(f"Número de Reynolds (Re_L): {Re_L:.2f}")
print(f"Número de Nusselt promedio (Nu_avg): {Nu_avg:.2f}")
print(f"Coeficiente de convección promedio h: {h_avg:.2f} W/m²K")
