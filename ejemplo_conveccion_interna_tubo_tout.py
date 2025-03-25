import CoolProp.CoolProp as CP
import numpy as np

def calcular_flujo_interno(T_in, T_wall, P, D, L, U, fluid):
    """
    Calcula la temperatura de salida del fluido y el coeficiente de convección promedio
    en un tubo de pared a temperatura fija.
    """
    A_s = np.pi * D * L  # Área de la superficie interna del tubo
    A_c = np.pi * (D/2)**2  # Área de la sección transversal
    
    # Iteración para encontrar T_out
    T_out = T_in  # Inicializamos con la temperatura de entrada
    error = 1e-6  # Criterio de convergencia
    max_iter = 100
    iter_count = 0
    
    while iter_count < max_iter:
        T_prom = (T_in + T_out) / 2  # Temperatura promedio
        
        # Propiedades del fluido a T_prom
        rho = CP.PropsSI('D', 'T', T_prom, 'P', P, fluid)
        mu = CP.PropsSI('V', 'T', T_prom, 'P', P, fluid)
        k = CP.PropsSI('CONDUCTIVITY', 'T', T_prom, 'P', P, fluid)
        Pr = CP.PropsSI('PRANDTL', 'T', T_prom, 'P', P, fluid)
        cp = CP.PropsSI('C', 'T', T_prom, 'P', P, fluid)
        
        # Flujo másico
        m_dot = rho * U * A_c
        
        # Número de Reynolds
        Re = (rho * U * D) / mu
        
        # Cálculo del número de Nusselt
        if Re < 2300:  # Flujo laminar
            Nu = 3.66  # Condición de pared a temperatura constante
        elif Re > 4000:  # Flujo turbulento (Dittus-Boelter)
            Nu = 0.023 * (Re**0.8) * (Pr**0.4)  
        else:
            Nu = np.nan  # Transición no modelada
        
        # Coeficiente de convección
        h = (Nu * k) / D
        
        # Cálculo de la nueva T_out a partir del balance de energía
        Q = h * A_s * (T_wall - T_prom)
#        T_out_new = T_in + Q / (m_dot * cp) 
        perim=np.pi*D
        T_out_new=(T_wall-T_in)*np.exp(-perim*L*h/(m_dot*cp))
        
        # Criterio de convergencia
        if abs(T_out_new - T_out) < error:
            break
        
        T_out = T_out_new
        iter_count += 1
    
    print(f"Temperatura de salida del fluido: {T_out:.2f} K")
    print(f"Coeficiente de convección promedio (h): {h:.2f} W/m²K")
    
    return T_out, h

# Parámetros de entrada
T_in = 300.15  # K (27°C)
T_wall = 350.15  # K (77°C)
P = 101325  # Pa
D = 0.02  # m (diámetro del tubo)
L = 1.0  # Longitud del tubo
U = 1.5  # m/s (velocidad promedio del fluido)
fluid = 'Air'

# Cálculo de la temperatura de salida y h
T_out, h_prom = calcular_flujo_interno(T_in, T_wall, P, D, L, U, fluid)
