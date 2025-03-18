import CoolProp.CoolProp as CP
import numpy as np

def propiedades_aire(T):
    """ Obtiene las propiedades del aire a la temperatura T usando CoolProp."""
    P = 101325  # Presión atmosférica en Pa
    rho = CP.PropsSI('D', 'T', T, 'P', P, 'Air')
    mu = CP.PropsSI('V', 'T', T, 'P', P, 'Air')
    k = CP.PropsSI('CONDUCTIVITY', 'T', T, 'P', P, 'Air')
    cp = CP.PropsSI('C', 'T', T, 'P', P, 'Air')
    nu = mu / rho  # Viscosidad cinemática
    alpha = k / (rho * cp)  # Difusividad térmica
    beta = 1 / T  # Aproximación para gases ideales
    Pr = CP.PropsSI('PRANDTL', 'T', T, 'P', P, 'Air')
    return rho, mu, k, cp, nu, alpha, beta, Pr

def coeficiente_conveccion_cilindro(T_s, T_inf, D, L, orientacion):
    """ Calcula el coeficiente de convección para un cilindro en aire más frío."""
    g = 9.81  # m/s²
    T_film = (T_s + T_inf) / 2  # Temperatura de película
    rho, mu, k, cp, nu, alpha, beta, Pr = propiedades_aire(T_film)
    
    if orientacion == "vertical":
        L_char = L  # Longitud característica: altura del cilindro
    elif orientacion == "horizontal":
        L_char = D  # Longitud característica: diámetro del cilindro
    else:
        raise ValueError("Orientación no válida. Usa 'vertical' o 'horizontal'")
    
    Ra = g * beta * (T_s - T_inf) * L_char**3 / (nu * alpha)  # Número de Rayleigh
    
    if orientacion == "vertical":
        Nu = (0.825 + (0.387 * Ra**(1/6)) / (1 + (0.492 / Pr)**(9/16))**(8/27))**2
    else:  # horizontal
        Nu = 0.36 + (0.518 * Ra**(1/4)) / (1 + (0.559 / Pr)**(9/16))**(4/9)
    
    h = (Nu * k) / L_char  # Coeficiente de convección
    
    print(f"Cilindro {orientacion.capitalize()} - Ra: {Ra:.2e}, Nu: {Nu:.2f}, h: {h:.2f} W/m²K")
    return h

# Parámetros de entrada
T_s = 350.15  # K (77°C)
T_inf = 300.15  # K (27°C)
D = 0.05  # m (diámetro del cilindro)
L = 0.5  # m (longitud del cilindro)

# Cálculo del coeficiente de convección
h_vertical = coeficiente_conveccion_cilindro(T_s, T_inf, D, L, "vertical")
h_horizontal = coeficiente_conveccion_cilindro(T_s, T_inf, D, L, "horizontal")
