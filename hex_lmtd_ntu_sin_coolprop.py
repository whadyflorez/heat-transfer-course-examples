
# Intercambiador doble-tubo: LMTD y NTU sin CoolProp (propiedades constantes)
from math import pi, log, exp

# --- Datos de proceso ---
m_dot_h = 0.50
m_dot_c = 0.60
Th_in = 80.0 + 273.15
Th_out_target = 50.0 + 273.15
Tc_in = 20.0 + 273.15

# Propiedades constantes (aprox agua 20-70°C)
rho_h = 985.0; cp_h = 4180.0; k_h = 0.64; mu_h = 0.0005
rho_c = 998.0; cp_c = 4180.0; k_c = 0.60; mu_c = 0.0010

Di = 0.020; Do = 0.024; Dan_i = 0.030
k_tube = 385.0; Rf_i = 0.0; Rf_o = 0.0

def dittus_boelter(Re, Pr, k, Dh, heating=False):
    n = 0.4 if not heating else 0.3
    Nu = 0.023 * (Re**0.8) * (Pr**n)
    h = Nu * k / Dh
    return Nu, h

Qdot = m_dot_h * cp_h * (Th_in - Th_out_target)
Tc_out = Tc_in + Qdot/(m_dot_c*cp_c)

Ch = m_dot_h * cp_h; Cc = m_dot_c * cp_c
Cmin = min(Ch, Cc); Cmax = max(Ch, Cc)

Ai = pi*(Di**2)/4.0
Ao = pi*(Dan_i**2 - Do**2)/4.0
ui = m_dot_h / (rho_h * Ai)
uo = m_dot_c / (rho_c * Ao)

Re_i = rho_h * ui * Di / mu_h
Pr_i = cp_h * mu_h / k_h
Re_o = rho_c * uo * (Dan_i - Do) / mu_c
Pr_o = cp_c * mu_c / k_c

Nu_i, h_i = dittus_boelter(Re_i, Pr_i, k_h, Di, heating=False)
Nu_o, h_o = dittus_boelter(Re_o, Pr_o, k_c, (Dan_i - Do), heating=True)

U_inv = (1.0/h_i) + log(Do/Di)/(2.0*k_tube) + (1.0/h_o)*(Di/Do) + Rf_i + Rf_o*(Di/Do)
U = 1.0 / U_inv

dT1 = (Th_in - Tc_out)
dT2 = (Th_out_target - Tc_in)
dTlm = (dT1 - dT2) / (log(dT1/dT2))

A_req = Qdot / (U * dTlm)
L_req = A_req / (pi * Di)

NTU = U * A_req / Cmin
Cr = Cmin / Cmax
eps = (1.0 - exp(-NTU*(1.0 - Cr))) / (1.0 - Cr*exp(-NTU*(1.0 - Cr)))
Q_NTU = eps * Cmin * (Th_in - Tc_in)
Th_out = Th_in - Q_NTU / Ch
Tc_out_chk = Tc_in + Q_NTU / Cc

def K2C(T): return T - 273.15
print("=== Intercambiador doble-tubo [sin CoolProp] ===")
print(f"Re_i={Re_i:.0f}, Re_o={Re_o:.0f}, Pr_i={Pr_i:.2f}, Pr_o={Pr_o:.2f}")
print(f"h_i={h_i:.1f} W/m2K, h_o={h_o:.1f} W/m2K, U={U:.1f} W/m2K")
print(f"Q={Qdot/1e3:.3f} kW, dTlm={dTlm:.2f} K")
print(f"A requerida={A_req:.4f} m2  ->  L={L_req:.3f} m (basado en A_i)")
print(f"NTU={NTU:.3f}, Cr={Cr:.3f}, eps={eps:.3f}")
print(f"Salidas: Th_out={K2C(Th_out):.2f} °C, Tc_out={K2C(Tc_out_chk):.2f} °C")
