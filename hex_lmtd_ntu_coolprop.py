
# Intercambiador doble-tubo: LMTD y NTU con propiedades CoolProp y Dittus-Boelter
from math import pi, log, exp
try:
    from CoolProp.CoolProp import PropsSI
except Exception as e:
    raise SystemExit("Se requiere CoolProp. Instale con 'pip install CoolProp'. Error: %s" % e)

# --- Datos de proceso ---
m_dot_h = 0.50           # kg/s (agua caliente) tubo interior
m_dot_c = 0.60           # kg/s (agua fría) anular
Th_in = 80.0 + 273.15    # K
Th_out_target = 50.0 + 273.15  # K
Tc_in = 20.0 + 273.15    # K
p = 1e5                  # Pa
# Geometría
Di = 0.020               # m
Do = 0.024               # m
Dan_i = 0.030            # m
k_tube = 385.0           # W/m-K
Rf_i = 0.0               # m^2 K/W
Rf_o = 0.0               # m^2 K/W

def water_props(T, p=1e5):
    rho = PropsSI("D", "T", T, "P", p, "Water")
    cp  = PropsSI("C", "T", T, "P", p, "Water")
    k   = PropsSI("L", "T", T, "P", p, "Water")
    mu  = PropsSI("V", "T", T, "P", p, "Water")
    return rho, cp, k, mu

def dittus_boelter(Re, Pr, k, Dh, heating=False):
    n = 0.4 if not heating else 0.3
    Nu = 0.023 * (Re**0.8) * (Pr**n)
    h = Nu * k / Dh
    return Nu, h

# --- Paso 1: Q y Tc_out por objetivo Th_out_target ---
# Propiedades a temperaturas de película estimadas
Tfilm_h = 0.5*(Th_in + Th_out_target)
rho_h, cp_h, k_h, mu_h = water_props(Tfilm_h, p)
Qdot = m_dot_h * cp_h * (Th_in - Th_out_target)
Tc_out = Tc_in + Qdot/(m_dot_c*cp_h)  # aproximar cp_c ~ cp_h la primera vez
Tfilm_c = 0.5*(Tc_in + Tc_out)
rho_c, cp_c, k_c, mu_c = water_props(Tfilm_c, p)
Tc_out = Tc_in + Qdot/(m_dot_c*cp_c)  # corregido

Ch = m_dot_h * cp_h
Cc = m_dot_c * cp_c
Cmin = min(Ch, Cc); Cmax = max(Ch, Cc)

# --- Paso 2: h_i, h_o ---
Ai = pi*(Di**2)/4.0
Ao = pi*(Dan_i**2 - Do**2)/4.0
ui = m_dot_h / (rho_h * Ai)
uo = m_dot_c / (rho_c * Ao)

Re_i = rho_h * ui * Di / mu_h
Pr_i = cp_h * mu_h / k_h
Re_o = rho_c * uo * (Dan_i - Do) / mu_c
Pr_o = cp_c * mu_c / k_c

Nu_i, h_i = dittus_boelter(Re_i, Pr_i, k_h, Di, heating=False)   # caliente se enfría -> n=0.4
Nu_o, h_o = dittus_boelter(Re_o, Pr_o, k_c, (Dan_i - Do), heating=True)  # frío se calienta -> n=0.3

# --- Paso 3: U y LMTD ---
U_inv = (1.0/h_i) + log(Do/Di)/(2.0*k_tube) + (1.0/h_o)*(Di/Do) + Rf_i + Rf_o*(Di/Do)
U = 1.0 / U_inv

dT1 = (Th_in - Tc_out)
dT2 = (Th_out_target - Tc_in)
dTlm = (dT1 - dT2) / (log(dT1/dT2))
A_req = Qdot / (U * dTlm)
L_req = A_req / (pi * Di)

# --- Verificación NTU ---
NTU = U * A_req / Cmin
Cr = Cmin / Cmax
eps = (1.0 - exp(-NTU*(1.0 - Cr))) / (1.0 - Cr*exp(-NTU*(1.0 - Cr)))
Q_NTU = eps * Cmin * (Th_in - Tc_in)
Th_out = Th_in - Q_NTU / Ch
Tc_out_chk = Tc_in + Q_NTU / Cc

def K2C(T): return T - 273.15

print("=== Intercambiador doble-tubo (agua-agua) ===")
print(f"Re_i={Re_i:.0f}, Re_o={Re_o:.0f}, Pr_i={Pr_i:.2f}, Pr_o={Pr_o:.2f}")
print(f"h_i={h_i:.1f} W/m2K, h_o={h_o:.1f} W/m2K, U={U:.1f} W/m2K")
print(f"Q={Qdot/1e3:.3f} kW, dTlm={dTlm:.2f} K")
print(f"A requerida={A_req:.4f} m2  ->  L={L_req:.3f} m (basado en A_i)")
print(f"NTU={NTU:.3f}, Cr={Cr:.3f}, eps={eps:.3f}")
print(f"Salidas: Th_out={K2C(Th_out):.2f} °C, Tc_out={K2C(Tc_out_chk):.2f} °C")
