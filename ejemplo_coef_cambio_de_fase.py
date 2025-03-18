import CoolProp.CoolProp as CP
import numpy as np

def propiedades_agua(T):
    """ Obtiene las propiedades del agua a la temperatura T."""
    P = 101325  # Presión atmosférica en Pa
    rho_l = CP.PropsSI('D', 'T', T, 'Q', 0, 'Water')
    rho_v = CP.PropsSI('D', 'T', T, 'Q', 1, 'Water')
    mu_l = CP.PropsSI('V', 'T', T, 'Q', 0, 'Water')
    k_l = CP.PropsSI('CONDUCTIVITY', 'T', T, 'Q', 0, 'Water')
    cp_l = CP.PropsSI('C', 'T', T, 'Q', 0, 'Water')
    h_fg = CP.PropsSI('H', 'T', T, 'Q', 1, 'Water') - CP.PropsSI('H', 'T', T, 'Q', 0, 'Water')
    nu_l = mu_l / rho_l  # Viscosidad cinemática
    return rho_l, rho_v, mu_l, k_l, cp_l, h_fg, nu_l

def ebullicion_rohsenow(T_s, T_sat, A):
    """ Calcula el coeficiente de convección en ebullición nucleada y la masa de vapor formado."""
    C_sf = 0.013  # Coeficiente para superficies metálicas
    g = 9.81  # Gravedad
    rho_l, rho_v, mu_l, k_l, cp_l, h_fg, nu_l = propiedades_agua(T_sat)
    Pr = cp_l * mu_l / k_l  # Número de Prandtl
    
    q = C_sf * (cp_l * (T_s - T_sat) / (h_fg * Pr**0.33))**3 * h_fg * rho_v * g  # Flujo de calor
    h = q / (T_s - T_sat)  # Coeficiente de convección
    m_vapor = (q * A) / h_fg  # Masa de vapor generado
    
    print(f"Ebullición: h = {h:.2f} W/m²K, q = {q:.2f} W/m², Vapor generado = {m_vapor:.4f} kg/s")
    return h, m_vapor

def condensacion_nusselt(T_s, T_sat, L, W):
    """ Calcula el coeficiente de convección en condensación en película y la masa de líquido condensado."""
    g = 9.81  # Gravedad
    rho_l, rho_v, mu_l, k_l, cp_l, h_fg, nu_l = propiedades_agua(T_sat)
    
    Gr = (rho_l * (rho_l - rho_v) * g * L**3) / (mu_l**2)  # Número de Grashof
    Pr = cp_l * mu_l / k_l  # Número de Prandtl
    Ra = Gr * Pr  # Número de Rayleigh
    Nu = 0.943 * (Ra**(1/4))  # Número de Nusselt
    h = (Nu * k_l) / L  # Coeficiente de convección
    
    q = h * (T_sat - T_s)  # Flujo de calor
    A = L * W  # Área de la placa
    m_liquido = (q * A) / h_fg  # Masa de líquido condensado
    
    print(f"Condensación: h = {h:.2f} W/m²K, q = {q:.2f} W/m², Líquido condensado = {m_liquido:.4f} kg/s")
    return h, m_liquido

# Parámetros de entrada para ebullición
T_s_eb = 380.15  # K (107°C)
T_sat_eb = 373.15  # K (100°C, agua a presión atmosférica)
A_eb = 0.1  # m² (área de la superficie)

# Parámetros de entrada para condensación
T_s_cond = 300.15  # K (27°C)
T_sat_cond = 373.15  # K (100°C)
L_cond = 0.5  # m (longitud de la placa)
W_cond = 0.3  # m (ancho de la placa)

# Cálculo de ebullición
h_eb, m_vapor = ebullicion_rohsenow(T_s_eb, T_sat_eb, A_eb)

# Cálculo de condensación
h_cond, m_liquido = condensacion_nusselt(T_s_cond, T_sat_cond, L_cond, W_cond)
