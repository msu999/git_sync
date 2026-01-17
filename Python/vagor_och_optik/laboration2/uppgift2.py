import numpy as np
from tabulate import tabulate
import matplotlib.pyplot as plt

theta = np.concatenate((np.arange(10, 50, 5), np.arange(50, 70, 2.5), np.arange(70, 85, 5)))

print(theta)
I_ort = np.array([48.4, 21.0, 8.1, 2.3, 0.7, 1.2, 2.6, 4.3, 6.1, 7.0, 7.6, 8.5, 9.1, 9.7, 10.4, 10.6, 10.9, 11.4, 10.7])
I_par = np.array([116.2, 89.0, 71.7, 56.0, 44.1, 35.9, 29.3, 24.2, 20.4, 18.2, 17.7, 16.8, 16.2, 15.5, 15.2, 14.2, 13.2, 12.0, 10.6])

headers = ["theta", "I_ort", "I_par"]
table = np.transpose(np.array([theta, I_ort, I_par]))

print(tabulate(table, headers))

theta_B = np.deg2rad(60)  # Från graf

n_2 = np.tan(theta_B)

print(f"n_2 från graf = {n_2}")

fig_ort, ax_ort = plt.subplots()
fig_par, ax_par = plt.subplots()

ax_ort.set(xlabel = "theta", ylabel = "I_ort")
ax_par.set(xlabel = "theta", ylabel = "I_par")

ax_ort.plot(theta, I_ort, "*", color = "red")
ax_ort.plot(theta, I_ort, "--", color = "grey")
ax_par.plot(theta, I_par, "*", color = "red")
ax_par.plot(theta, I_par, "--", color = "grey")

plt.show()
