import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate

U = 2.5  # Satt i laborationen
L = 83.3e-3  # Mätt
R_L = 94.2  # Mätt
C = 22.2e-9  # Mätt
R = 1000

f = np.array([1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000])
w = 2*np.pi*f

U_R = np.array([0.736/2, 1.72/2, 3.36/2, 3.76/2, 2.88/2, 2.12/2, 1.4/2, 1.04/2])
T = np.array([992e-6, 512e-6, 332e-6, 252e-6, 202e-6, 166e-6, 126e-6, 100e-6])
dt = np.array([216e-6, 96e-6, 32e-6, 10e-6, 22e-6, 26e-6, 24e-6, 20e-6])

# Ström från utgång 2
I = U_R/R

# Fasförskjutning i grader
phi = 360*dt/T

print(f"phi = {phi}")

headers = ["f", "U_R", "I", "dt", "T", "phi"]
table = np.transpose(np.array([f, U_R, I, dt, T, phi]))

print(tabulate(table, headers))

N = 100

# Numpyvektorer för teoretiska beräkningar
fs = np.linspace(1, 1e4, N)
ws = 2*np.pi*fs
Zs = R + R_L + 1j*(ws*L - 1/(ws*C))
Is = U/abs(Zs)
phis = np.rad2deg(np.arctan((ws*L - 1/(ws*C))/(R + R_L)))

phi_corr = np.array([-78.38709677, -67.5,        -34.69879518, -14.28571429, 39.20792079, 56.38554217,
 68.57142857, 72.0        ])

# Plottning
fig_I_theory, ax_I_theory = plt.subplots()
fig_phi_theory, ax_phi_theory = plt.subplots()
fig_I_exp, ax_I_exp = plt.subplots()
fig_phi_exp, ax_phi_exp = plt.subplots()

ax_I_theory.set(xlabel="frequency (hz)", ylabel="current (A)", title = "Theoretical I")
ax_phi_theory.set(xlabel="frequency (hz)", ylabel="phase shift (degrees)", title = "Theoretical phi")
ax_I_exp.set(xlabel="frequency (hz)", ylabel="current (A)", title = "Experimental I")
ax_phi_exp.set(xlabel="frequency (hz)", ylabel="phase shift (degrees)", title = "Experimental phi")

ax_I_theory.plot(fs, Is)
ax_phi_theory.plot(fs, phis)
ax_I_exp.plot(f, I)
ax_phi_exp.plot(f, phi_corr)

plt.show()
