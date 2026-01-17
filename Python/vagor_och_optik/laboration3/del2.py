import numpy as np
from tabulate import tabulate

l_röd = 670e-9  # Står på sidan av lasern
l_grön = 532e-9  # Står på sidan av lasern

L = 197e-2  # Mät

# Uppgift 5:
print("Uppgift 5:")

r_röd = np.array([4e-3, 7.8e-3])  # Radier för första mörka ring
r_grön = np.array([3.2e-3, 6e-3])  # Radier för första mörka ring

D_röd = 1e3*1.22*l_röd*L/r_röd
D_grön = 1e3*1.22*l_grön*L/r_grön

headers = ["(Diametrar) [m]", "Hål 1", "Hål 2"]
columns = np.array(["Röd laser", "Grön laser"])
# table = np.transpose(np.array([columns, D_röd, D_grön]))

print(f"Grön laser:\n Diameter hål 1 = {D_grön[0]}, Diameter hål 2 = {D_grön[1]}\nRöd laser:\n Diameter hål 1 = {D_röd[0]}, Diameter hål 2 = {D_röd[1]}")

# Uppgift 6:
print("\nUppgift 6:")

theta_1 = np.deg2rad(42)  # Mät
d = l_röd/np.sin(theta_1)
N = 1/(d*1e3)

print(f"Avstånd mellan ritsar = {d} m\nAntal ritsar per mm = {N}")


