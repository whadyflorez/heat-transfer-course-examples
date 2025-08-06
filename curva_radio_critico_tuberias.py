import numpy as np
import matplotlib.pyplot as plt

# --- Parámetros conocidos ---
T_i = 120          # Temperatura de la pared interna (°C)
T_inf = 25         # Temperatura del ambiente (°C)
r1 = 0.05          # Radio externo de la tubería (m)
k = 0.72          # Conductividad térmica del aislamiento (W/m·K)
h = 10            # Coeficiente de convección externo (W/m²·K)
L = 1.0            # Longitud de la tubería (m)

# --- Radios exteriores del aislamiento (r2) ---
r2_array = np.linspace(r1 + 0.001, 0.20, 200)  # desde un poco más que r1 hasta 20 cm

# --- Cálculo de la pérdida de calor q para cada r2 ---
q_array = []
for r2 in r2_array:
    R_cond = np.log(r2 / r1) / (2 * np.pi * k * L)
    R_conv = 1 / (h * 2 * np.pi * r2 * L)
    R_total = R_cond + R_conv
    q = (T_i - T_inf) / R_total
    q_array.append(q)

# --- Gráfica ---
plt.figure(figsize=(8,6))
plt.plot(r2_array * 100, q_array, color='firebrick', linewidth=2)
plt.xlabel('Radio exterior del aislamiento [cm]', fontsize=12)
plt.ylabel('Pérdida de calor q [W]', fontsize=12)
plt.title('Variación de la pérdida de calor con el radio del aislamiento', fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.show()
