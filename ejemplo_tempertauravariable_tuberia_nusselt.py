import CoolProp.CoolProp as CP
import numpy as np

def calcular_h_tubo(T_in, T_out, P, D, L, U, fluid):
    """
    Calcula el coeficiente de transferencia de calor promedio para flujo interno
    en un tubo de diámetro D y longitud L con un fluido en régimen laminar o turbulento.
    """
    # Temperatura de película (promedio entre entrada y salida)
    T_avg = (T_in + T_out) / 2
    
    # Propiedades del fluido a la temperatura de película
    rho = CP.PropsSI('D', 'T', T_avg, 'P', P, 'Air')
    mu = CP.PropsSI('V', 'T', T_avg, 'P', P, 'Air')
    k = CP.PropsSI('CONDUCTIVITY', 'T', T_avg, 'P', P, 'Air')
    Pr = CP.PropsSI('PRANDTL', 'T', T_avg, 'P', P, 'Air')

    # Cálculo del número de Reynolds basado en la longitud total del tubo
    Re_L = (rho * U * D) / mu

    # Cálculo del Número de Nusselt promedio según régimen de flujo
    if Re_L < 2300:  # Flujo laminar
        Nu_L = 3.66 + (0.065 * (Re_L * Pr * (D/L)) / (1 + 0.04 * (Re_L * Pr * (D/L))**(2/3)))
    elif Re_L > 4000:  # Flujo turbulento (Dittus-Boelter)
        if Pr < 0.7:
            Nu_L = 0.023 * (Re_L**0.8) * (Pr**0.4)  # Para calentamiento
        else:
            Nu_L = 0.023 * (Re_L**0.8) * (Pr**(1/3))  # Para enfriamiento o neutro
    else:
        Nu_L = np.nan  # Caso intermedio sin correlación específica
    
    # Cálculo del coeficiente de convección promedio
    h_L = (Nu_L * k) / D
    
    print(f"Temperatura promedio: {T_avg:.2f} K")
    print(f"Número de Reynolds (Re_L): {Re_L:.2f}")
    print(f"Número de Nusselt promedio (Nu_L): {Nu_L:.2f}")
    print(f"Coeficiente de transferencia de calor (h): {h_L:.2f} W/m²K")
    
    return h_L

# Parámetros de entrada
T_in = 300.15  # K (27°C)
T_out = 330.15  # K (57°C)
T_wall = 340.15  # K (67°C)
P = 101325  # Pa
D = 0.02  # m (diámetro del tubo)
L = 1.0  # Longitud en m
U = 1.5  # m/s
fluid = 'Air'

# Cálculo del coeficiente de transferencia de calor promedio
h_prom = calcular_h_tubo(T_in, T_out, P, D, L, U, fluid)

