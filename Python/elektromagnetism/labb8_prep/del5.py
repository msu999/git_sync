import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate

L = 83.3e-3  # Mätt
R_L = 94.2  # Mätt
C = 22.2e-9  # Mätt

w_r = np.sqrt(1/(L*C) - (R_L/L)**2)
f_r = w_r/(2*np.pi)
f_c = 1/(2*np.pi*R_L*C)

print(f"w_r = {w_r}")
print(f"f_r = {f_r}")
print(f"f_c = {f_c}")

f = np.array([1e3, 2e3, 3e3, 3.5e3, 3.6e3, 3.7e3, 3.8e3, 3.9e3, 4e3, 5e3, 6e3, 7e3, 8e3, 9e3, 10e3])
w = 2*np.pi*f

U = np.array([1.690, 1.718, 1.724, 1.726, 1.728, 1.728, 1.728, 1.728, 1.728, 1.727, 1.725, 1.725, 1.725, 1.725, 1.726])  # Effektivvärde
I = np.array([0.107, 0.111, 0.115, 0.107, 0.105, 0.104, 0.104, 0.104, 0.105, 0.109, 0.110, 0.111, 0.112, 0.113, 0.114])  # Effektivvärde

headers = ["w", "I"]
table = np.transpose(np.array([f, I]))

print(tabulate(table, headers))

fig, ax = plt.subplots()

ax.set(xlabel = "Frequency (hz)", ylabel = "Current (A)")

ax.plot(f, I, "--", color = "gray")
ax.plot(f, I, "*", color = "red")

plt.show()
