import numpy as np
import matplotlib.pyplot as plt

# --- Parámetros ---
T_i = 120      # °C
T_inf = 25     # °C
r1 = 0.05      # m, radio externo del tubo
k1 = 0.72      # W/m·K, conductividad 1er aislante
k2 = 0.9      # W/m·K, conductividad 2do aislante
h = 10         # W/m²·K
L = 1.0        # m, longitud

# --- Rango de radios ---
r2_vals = np.linspace(r1 + 0.001, 0.30, 100)    # radio exterior de la primera capa
r3_vals = np.linspace(0.30, 0.50, 100)          # radio exterior de la segunda capa

# Crear mallas 2D
R2, R3 = np.meshgrid(r2_vals, r3_vals)

# Asegurar que r3 > r2
mask = R3 > R2

# Calcular resistencias y q
R_cond1 = np.log(R2 / r1) / (2 * np.pi * k1 * L)
R_cond2 = np.log(R3 / R2) / (2 * np.pi * k2 * L)
R_conv  = 1 / (h * 2 * np.pi * R3 * L)

R_total = np.where(mask, R_cond1 + R_cond2 + R_conv, np.nan)
q = (T_i - T_inf) / R_total

# --- Gráfica de mapa de calor ---
plt.figure(figsize=(8,6))
cp = plt.contourf(R2*100, R3*100, q, levels=50, cmap='viridis')
# Líneas de contorno superpuestas (blancas o negras, más gruesas)
contour_lines = plt.contour(R2*100, R3*100, q, levels=10, colors='k', linewidths=1.2)
# Etiquetas en las líneas de contorno (opcional)
plt.clabel(contour_lines, inline=True, fontsize=9, fmt="%.0f")
cbar = plt.colorbar(cp)
cbar.set_label('Pérdida de calor q [W]')
plt.xlabel('Radio r₂ del 1er aislante [cm]')
plt.ylabel('Radio r₃ del 2do aislante [cm]')
plt.title('Pérdidas de calor en función de ambos radios de aislamiento')
plt.grid(True)
plt.tight_layout()
plt.show()
