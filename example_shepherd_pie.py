import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# --------------------------------
# Propiedades del Shepherd's Pie
# --------------------------------
rho = 1052.0          # kg/m3
cp = 3606.0           # J/kgK
k = 0.505             # W/mK
alpha = 1.33e-7       # m2/s

# Geometría térmica: pared plana
L = 0.02              # m (media pared)

# Condiciones
Ti = 20.0             # °C
Tinf = 180.0          # °C
h = 30.0              # W/m2K

# Geometría total del pastel para energía
Rpie = 0.10           # m
Hpie = 0.04           # m
V = np.pi * Rpie**2 * Hpie

# --------------------------------
# Biot
# --------------------------------
Bi = h * L / k

# --------------------------------
# Primera raíz zeta1: zeta*tan(zeta)=Bi
# --------------------------------
def eq_zeta(z):
    return z * np.tan(z) - Bi

zeta1 = fsolve(eq_zeta, 0.9)[0]

# Coeficiente C1
C1 = 4*np.sin(zeta1) / (2*zeta1 + np.sin(2*zeta1))

# --------------------------------
# Funciones
# --------------------------------
def Fo(t):
    return alpha * t / L**2

def theta0_star(t):
    return C1 * np.exp(-(zeta1**2) * Fo(t))

def temperature(x, t):
    theta = theta0_star(t) * np.cos(zeta1 * x / L)
    return Tinf + theta * (Ti - Tinf)

def Q_total(t):
    Q0 = rho * cp * V * (Tinf - Ti)
    ratio = 1.0 - (np.sin(zeta1)/zeta1) * theta0_star(t)
    return Q0 * ratio

# --------------------------------
# 1. T(0,t)
# --------------------------------
t = np.linspace(0, 3600, 300)
Tcenter = np.array([temperature(0.0, tt) for tt in t])

plt.figure()
plt.plot(t/60, Tcenter, linewidth=2)
plt.xlabel("Tiempo (min)")
plt.ylabel("Temperatura en el centro (°C)")
plt.title("T(0,t) en el Shepherd's Pie")
plt.grid(True)
plt.show()

# --------------------------------
# 2. Perfil T(x,t) a 30 min
# --------------------------------
t_profile = 1800.0
x = np.linspace(0, L, 200)
Tx = np.array([temperature(xx, t_profile) for xx in x])

plt.figure()
plt.plot(x*100, Tx, linewidth=2)
plt.xlabel("Posición desde el centro x (cm)")
plt.ylabel("Temperatura (°C)")
plt.title("Perfil T(x,t) a t = 30 min")
plt.grid(True)
plt.show()

# --------------------------------
# 3. Calor total absorbido Q(t)
# --------------------------------
Q = np.array([Q_total(tt) for tt in t])

plt.figure()
plt.plot(t/60, Q/1000, linewidth=2)
plt.xlabel("Tiempo (min)")
plt.ylabel("Calor absorbido Q(t) (kJ)")
plt.title("Calor total absorbido por el pastel")
plt.grid(True)
plt.show()

# --------------------------------
# 4. Valores numéricos a 30 min
# --------------------------------
print(f"Bi = {Bi:.4f}")
print(f"zeta1 = {zeta1:.4f}")
print(f"C1 = {C1:.4f}")
print(f"Fo(30 min) = {Fo(1800):.4f}")
print(f"theta0*(30 min) = {theta0_star(1800):.4f}")
print(f"Tcentro(30 min) = {temperature(0,1800):.2f} °C")
print(f"Q(30 min) = {Q_total(1800)/1000:.2f} kJ")