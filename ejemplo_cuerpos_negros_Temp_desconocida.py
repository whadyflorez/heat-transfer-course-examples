import numpy as np
from scipy.optimize import fsolve

# Constante de Stefan-Boltzmann
sigma = 5.670374419e-8  # W/m²·K⁴

# Áreas
A = np.array([1.0, 2.0, 3.0])

# Temperaturas conocidas
T1 = 800.0  # K
T3 = 400.0  # K

# Temperatura desconocida: T2
# Supuesto inicial
T2_guess = 500.0

# Calor neto en superficie 2
Q2_deseado = 500.0  # W

# Factores de forma
F = np.array([
    [0.5, 0.2, 0.3],
    [0.1, 0.6, 0.3],
    [0.1, 0.2, 0.7]
])

def ecuacion_T2(T2):
    T = [T1, T2, T3]
    Q2 = 0.0
    for j in range(3):
        if j != 1:
            Q2 += A[1] * F[1, j] * sigma * (T[1]**4 - T[j]**4)
    return Q2 - Q2_deseado

# Resolver la ecuación no lineal
T2_sol = fsolve(ecuacion_T2, T2_guess)[0]

# Calcular los otros Q para ver consistencia
T = [T1, T2_sol, T3]
Q = np.zeros(3)
for i in range(3):
    for j in range(3):
        if i != j:
            Q[i] += A[i] * F[i, j] * sigma * (T[i]**4 - T[j]**4)

# Resultados
print(f"Temperatura T2 que produce Q2 = 500 W: {T2_sol:.2f} K\n")
for i in range(3):
    print(f"Q{i+1} (flujo neto desde superficie {i+1}): {Q[i]:.2f} W")

print(f"\nSuma total de flujos: {np.sum(Q):.2e} W (≈ 0)")
