
# shell_tube_lmtd_ntu_coolprop.py
# Ejemplo didáctico de intercambiador coraza-tubos: cálculos LMTD (con F) y NTU
# Requiere: pip install CoolProp
from math import pi, log, exp
try:
    from CoolProp.CoolProp import PropsSI
except Exception as e:
    raise SystemExit("Instale CoolProp: pip install CoolProp. Error: %s" % e)

# --------------------- Datos de proceso (ejemplo) ----------------------
m_dot_s = 1.20        # kg/s (shell - fluido caliente)
m_dot_t = 2.00        # kg/s (tubes - fluido frío)
Ts_in = 120.0 + 273.15  # K (shell hot in)
Ts_out_target = 80.0 + 273.15 # K (deseada)
Tt_in = 20.0 + 273.15  # K (tube cold in)

# Geometría propuesta
d_o = 0.01905         # m (diámetro externo del tubo, 3/4'')
t_wall = 0.00124      # m
d_i = d_o - 2*t_wall  # m
D_shell = 0.35        # m (diámetro interno de la coraza)
L_tube = 3.0          # m (longitud efectiva)
N_t = 50              # número de tubos
k_tube = 16.0         # W/mK (acero inoxidable, aproximado)
Rf_s = 1e-4           # m2K/W (fouling shell)
Rf_t = 1e-4           # m2K/W (fouling tube)

# Factor de corrección LMTD (usuario puede modificar)
F_corr = 0.85         # valor típico inicial para 1 shell pass / 2 tube passes

# --------------------- Funciones auxiliares ----------------------------
def props(fluid, T, p=1e5):
    rho = PropsSI("D", "T", T, "P", p, fluid)
    cp  = PropsSI("C", "T", T, "P", p, fluid)
    k   = PropsSI("L", "T", T, "P", p, fluid)
    mu  = PropsSI("V", "T", T, "P", p, fluid)
    return rho, cp, k, mu

def dittus_boelter(Re, Pr, k, Dh, heating=False):
    n = 0.4 if not heating else 0.3
    Nu = 0.023 * (Re**0.8) * (Pr**n)
    h = Nu * k / Dh
    return Nu, h

# --------------------- Estimación inicial de Q y temperaturas -------------
# Estimamos propiedades en temperaturas de película (iteración sencilla)
Tfilm_s = 0.5*(Ts_in + Ts_out_target)
rho_s, cp_s, k_s, mu_s = props("Water", Tfilm_s)
# primera aproximación de Q por el lado shell (m_dot * cp * dT)
Qdot = m_dot_s * cp_s * (Ts_in - Ts_out_target)
# estimar Tt_out inicial
rho_t, cp_t, k_t, mu_t = props("Water", 0.5*(Tt_in + (Tt_in + Qdot/(m_dot_t*cp_s))))
Tt_out = Tt_in + Qdot/(m_dot_t*cp_t)

# --------------------- Cálculos geométricos -----------------------------
A_o = pi * d_o * L_tube * N_t   # área total de transferencia (exterior tubos)
A_i_tube = pi * d_i**2 / 4.0
A_shell_flow = max(1e-6, pi*(D_shell**2)/4.0 - N_t * pi * d_o**2 / 4.0)  # área aproximada de paso de shell
D_h_shell = max(1e-6, D_shell - N_t * d_o / pi)  # estimación simple del diámetro hidráulico de shell flow

# --------------------- Velocidades y números adimensionales ----------------
u_t = m_dot_t / (rho_t * N_t * A_i_tube)   # velocidad en cada tubo (media) si todos los tubos participan
Re_t = rho_t * u_t * d_i / mu_t
Pr_t = cp_t * mu_t / k_t

u_s = m_dot_s / (rho_s * A_shell_flow)
Re_s = rho_s * u_s * D_h_shell / mu_s
Pr_s = cp_s * mu_s / k_s

# --------------------- Coeficientes convectivos -------------------------
Nu_t, h_t = dittus_boelter(Re_t, Pr_t, k_t, d_i, heating=True)   # el fluido en tubo se calienta -> n=0.3
Nu_s = 0.36 * (Re_s**0.55) * (Pr_s**0.33)
h_s = Nu_s * k_s / D_h_shell

# --------------------- Coeficiente global U ------------------------------
U_inv = 1.0/h_s + (d_o/ d_i) * (log(d_o/d_i) / (2.0 * k_tube)) + (d_o/ d_i) * (1.0/h_t) + Rf_s + Rf_t * (d_o/ d_i)
U = 1.0 / U_inv

# --------------------- LMTD -----------------------------------------------
dT1 = (Ts_in - Tt_out)
dT2 = (Ts_out_target - Tt_in)
if dT1<=0 or dT2<=0:
    raise SystemExit("Diferencias de temperatura no válidas. Revisa entradas.")
dT_lm = (dT1 - dT2) / (log(dT1/dT2))
dT_lm_corr = F_corr * dT_lm

A_req = Qdot / (U * dT_lm_corr)

# Verificación con área disponible
A_available = A_o
L_required = A_req / (pi * d_o)  # longitud equivalente si cambiamos número de tubos
area_ok = A_available >= A_req

# --------------------- NTU -----------------------------------------------
Ch = m_dot_s * cp_s
Ct = m_dot_t * cp_t
Cmin = min(Ch, Ct)
Cmax = max(Ch, Ct)
NTU = U * A_available / Cmin
Cr = Cmin / Cmax
eps = (1.0 - exp(-NTU*(1.0 - Cr))) / (1.0 - Cr * exp(-NTU*(1.0 - Cr)))
Q_NTU = eps * Cmin * (Ts_in - Tt_in)
Ts_out = Ts_in - Q_NTU / Ch
Tt_out_chk = Tt_in + Q_NTU / Ct

# --------------------- Reporte ------------------------------------------
def K2C(T): return T - 273.15
print("=== Intercambiador coraza-tubos (ejemplo didáctico) ===")
print(f"Re_t={Re_t:.0f}, Re_s={Re_s:.0f}, Pr_t={Pr_t:.3f}, Pr_s={Pr_s:.3f}")
print(f"h_t={h_t:.1f} W/m2K, h_s={h_s:.1f} W/m2K, U={U:.1f} W/m2K")
print(f"Q (diseño) = {Qdot/1e3:.3f} kW")
print(f"dT_lm = {dT_lm:.2f} K, F = {F_corr}, dT_lm_corr = {dT_lm_corr:.2f} K")
print(f"A requerida = {A_req:.4f} m2, Área disponible = {A_available:.4f} m2, suficiente? {area_ok}")
print(f"NTU = {NTU:.3f}, Cr = {Cr:.3f}, eps = {eps:.3f}, Q_NTU = {Q_NTU/1e3:.3f} kW")
print(f"Temperaturas (LMTD target): Ts_out_target = {K2C(Ts_out_target):.2f} °C, Tt_out (estim) = {K2C(Tt_out):.2f} °C")
print(f"Temperaturas (NTU result): Ts_out = {K2C(Ts_out):.2f} °C, Tt_out = {K2C(Tt_out_chk):.2f} °C")
print(f"Área total externa tubos (A_o) = {A_o:.4f} m2, Longitud de tubo propuesta = {L_tube:.2f} m, N_t = {N_t}")

print('\\nAdvertencias:')
print('- La correlación shell-side usada es una simplificación (Kern-style).')
print('- Para proyecto real use Bell-Delaware y tablas de factor de corrección F específicas.')
