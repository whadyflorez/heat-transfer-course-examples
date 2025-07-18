import numpy as np
import matplotlib.pyplot as plt

# Constante de Stefan-Boltzmann
sigma = 5.670374419e-8  # W/m²·K⁴

# Datos del problema
T1 = 800  # K (pared caliente)
T3 = 300  # K (pared fría)
A = 1.0   # m² (área de cada superficie)

# Cálculo de T2 (temperatura del escudo)
T2_4 = (T1**4 + T3**4) / 2
T2 = T2_4**0.25

# Calor transferido
Q = A * sigma * (T1**4 - T2**4)

# Resultados
print(f"Temperatura del escudo (T2): {T2:.2f} K")
print(f"Flujo de calor por radiación (Q): {Q:.2f} W")

# Gráfica ilustrativa
temps = [T1, T2, T3]
labels = ['Pared caliente', 'Escudo', 'Pared fría']
plt.plot(labels, temps, marker='o')
plt.title('Temperaturas de las superficies (cuerpos negros)')
plt.ylabel('Temperatura (K)')
plt.grid(True)
plt.show()
