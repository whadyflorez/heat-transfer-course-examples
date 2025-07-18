import numpy as np

# Constante de Stefan-Boltzmann
sigma = 5.670374419e-8

# Geometría
r1 = 0.15   # m
r2 = 0.20  # m
L = 0.001    # m

# Áreas
A1 = np.pi * r1**2
A2 = np.pi * r2**2
A3 = 10  # ambiente grande, envolvente

# Emisividades
eps1 = 0.8
eps2 = 0.2
eps3 = 1.0

# Temperaturas (K)
T1 = 300+273.15
T2 = 80+273.15
T3 = 300.0

# Función de forma F_12 usando fórmula para discos coaxiales
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

# Radiosidad del ambiente (negro)
J3 = sigma * T3**4

# Energías emitidas
Eb1 = sigma * T1**4
Eb2 = sigma * T2**4

# Resistencias superficiales
Rs1 = (1 - eps1) / (eps1 * A1)
Rs2 = (1 - eps2) / (eps2 * A2)

# Resistencias espaciales (geométricas)
R12 = 1 / (A1 * F12)
R13 = 1 / (A1 * F13)
R21 = 1 / (A2 * F21)
R23 = 1 / (A2 * F23)

# Sistema lineal en J1 y J2
# Definimos la matriz A y el vector b: A · J = b
A = np.array([
    [1/Rs1 + 1/R12 + 1/R13, -1/R12],
    [-1/R21, 1/Rs2 + 1/R21 + 1/R23]
])

b = np.array([
    Eb1 / Rs1 + J3 / R13,
    Eb2 / Rs2 + J3 / R23
])

# Resolver radiosidades J1 y J2
J1, J2 = np.linalg.solve(A, b)

# Calcular flujos netos
q1 = (J1 - J2)/R12 + (J1 - J3)/R13
q2 = (J2 - J1)/R21 + (J2 - J3)/R23
q3 = (J3 - J1)/R13 + (J3 - J2)/R23  # por conservación

# Mostrar resultados
print("Radiosidades:")
print(f"J1 = {J1:.2e} W/m²")
print(f"J2 = {J2:.2e} W/m²")
print(f"J3 = {J3:.2e} W/m²")

print("\nFlujos netos de calor:")
print(f"q1 = {q1:.2f} W")
print(f"q2 = {q2:.2f} W")
print(f"q3 = {q3:.2f} W")

print(f"\nSuma total de flujos: {q1 + q2 + q3:.2e} W (≈ 0 si todo está bien)")
