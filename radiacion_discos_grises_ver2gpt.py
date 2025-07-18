import numpy as np
from scipy.optimize import fsolve

# Constante de Stefan-Boltzmann
sigma = 5.670374419e-8

# Geometría
r1 = 0.1  # m
r2 = 0.5 # m
L = 0.1   # m

# Emisividades
eps1 = 0.8
eps2 = 0.6
eps3 = 1.0  # ambiente

# Temperaturas
T1 = 800.0  # K
T2 = 600.0  # K
T3 = 300.0  # K

# Áreas
A1 = np.pi * r1**2
A2 = np.pi * r2**2
A3 = 1000.0  # ambiente grande

# Cálculo de F12 (discos coaxiales)
Ri = r1 / L
Rj = r2 / L
S = 1 + (1 + Rj**2) / Ri**2
F12 = 0.5 * (S - np.sqrt(S**2 - 4 * (r2/r1)**2))
F21 = (A1 * F12) / A2
F13 = 1 - F12
F23 = 1 - F21
F31 = (A1 * F13) / A3
F32 = (A2 * F23) / A3
F33 = 1 - F31 - F32

# J3 conocido (cuerpo negro)
J3 = sigma * T3**4

# Sistema de radiosidades corregido
def sistema_radiosidades(x):
    J1, J2 = x
    G1 = F12 * J2 + F13 * J3
    G2 = F21 * J1 + F23 * J3
    eq1 = J1 - eps1 * sigma * T1**4 - (1 - eps1) * G1
    eq2 = J2 - eps2 * sigma * T2**4 - (1 - eps2) * G2
    return [eq1, eq2]

# Resolver el sistema
sol = fsolve(sistema_radiosidades, [sigma * T1**4, sigma * T2**4])
J1, J2 = sol

# Irradiaciones
G1 = F12 * J2 + F13 * J3
G2 = F21 * J1 + F23 * J3
G3 = F31 * J1 + F32 * J2 + F33 * J3

# Flujos netos
Q1 = A1 * (J1 - G1)
Q2 = A2 * (J2 - G2)
Q3 = A3 * (J3 - G3)

# Resultados
print(f"Radiosidades:")
print(f"J1 = {J1:.2e} W/m²")
print(f"J2 = {J2:.2e} W/m²")
print(f"J3 = {J3:.2e} W/m² (cuerpo negro)")

print(f"\nFlujos netos:")
print(f"Q1 = {Q1:.2f} W")
print(f"Q2 = {Q2:.2f} W")
print(f"Q3 = {Q3:.2f} W")
print(f"Suma total Q1+Q2+Q3 = {Q1+Q2+Q3:.2e} W (≈ 0 si está bien)")
