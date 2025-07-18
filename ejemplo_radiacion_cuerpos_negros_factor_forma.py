import numpy as np

# Constante de Stefan-Boltzmann
sigma = 5.670374419e-8  # W/m²·K⁴

# Áreas (m²) y temperaturas (K)
A = np.array([1.0, 2.0, 3.0])
T = np.array([800.0, 600.0, 400.0])

# Factores de forma F_ij (validados por reciprocidad)
F = np.array([
    [0.5, 0.2, 0.3],
    [0.1, 0.6, 0.3],
    [0.1, 0.2, 0.7]
])

# Cálculo del flujo neto Q_i para cada superficie
Q = np.zeros(3)
for i in range(3):
    for j in range(3):
        if i != j:
            Q[i] += A[i] * F[i, j] * sigma * (T[i]**4 - T[j]**4)

# Mostrar resultados
for i in range(3):
    print(f"Q{i+1} (flujo neto desde superficie {i+1}): {Q[i]:.2f} W")

# Verificación: suma total debe ser cero
print(f"\nSuma total de flujos (Q1+Q2+Q3): {np.sum(Q):.2e} W (≈ 0)")
