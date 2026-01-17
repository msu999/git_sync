import numpy as np
from tabulate import tabulate

L_meas = 136.5e-3

f = 50  # Standardspänning i vägguttag
w = 2*np.pi*f

U = np.array([25.1, 26.4, 25.8, 25.8])
I = np.array([0.751, 0.041, 0.311, 0.247])
P = np.array([18.95, 0.002, 4.16, 4.68])

S = U*I
cosphi = P/S
phi = np.arccos(P/S)

headers = ["U", "I", "P", "S", "phi", "cos(phi)"]
part = np.array(["(a)", "(b)", "(c)", "(d)"])
table = np.transpose(np.array([part, U, I, P, S, phi, cosphi]))

print(tabulate(table, headers))

# Fråga 1
C = I[1]/(U[1]*w)

print(f"Beräknat C = {C} F")

# Fråga 2
R = 33
R_L = 5

L = np.sqrt((U[2]/I[2])**2 - (R + R_L)**2)/w

print(f"Beräknat L = {L} H")
