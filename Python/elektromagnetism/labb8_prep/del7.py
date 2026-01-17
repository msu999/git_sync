import numpy as np
from tabulate import tabulate
import matplotlib.pyplot as plt

f = 50  # Standardfrekvens i vägguttag
w = 2*np.pi*f

R = 33
R_L = 5
R_t = R + R_L
L = 0.235  # Mer exakt värde från del 6
C = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60])  # Varieras med box

U = np.array([25.8, 25.9, 25.9, 25.8, 25.9, 25.9, 25.9, 26.0, 26.0, 26.0, 26.1, 26.1, 26.1])
I = np.array([0.311, 0.278, 0.248, 0.220, 0.201, 0.185, 0.175, 0.176, 0.186, 0.203, 0.229, 0.259, 0.287])

P = np.array([4.17, 4.17, 4.18, 4.20, 4.21, 4.20, 4.20, 4.22, 4.24, 4.25, 4.26, 4.24, 4.28])

cosphi = P/(U*I)

headers = ["C", "U", "I", "P", "cos(phi)"]
table = np.transpose(np.array([C, U, I, P, cosphi]))

print(tabulate(table, headers))

C_crit = L/(L*w**2 - R_t**2)

print(f"C för bäst faskompensering = {C_crit}")

fig, ax = plt.subplots()

ax.plot(C, cosphi)

plt.show()
