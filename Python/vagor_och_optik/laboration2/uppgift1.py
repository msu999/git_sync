import numpy as np
from tabulate import tabulate
import matplotlib.pyplot as plt

I_0 = 131.3  # Mät

theta = np.arange(-90, 100, 10)
I = np.array([2.3, 7.5, 20.0, 40.0, 63.7, 87.0, 108.8, 125.0, 134.1, 137.7, 136.0, 124.4, 108.8, 87.0, 63.0, 40.0, 26.9, 7.6, 2.4])  # Mät
qout = I/I_0

headers = ["theta", "I", "I/I_0"]
table = np.transpose(np.array([theta, I, qout]))

print(tabulate(table, headers))

thetas = np.linspace(-np.pi, np.pi, 100)

fig_theory, ax_theory = plt.subplots()
fig_exp, ax_exp = plt.subplots()

ax_theory.set(xlabel="theta (deg)", ylabel="I")
ax_exp.set(xlabel="theta (deg)", ylabel="I (V)")

ax_theory.plot(thetas, np.cos(thetas)**2)
ax_exp.plot(theta, qout, "*")
ax_exp.plot(theta, qout, color = "grey")

plt.show()
