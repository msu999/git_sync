import numpy as np

# Konstant genom hela laborationen:
l = 532e-9  # Anges på lasern
L = 0.9  # Mät

# Uppgift 1:
print("\nUppgift 1:")

x = (8.65 - 7.4)/2*1e-2  # Mät
a = (l*L/x)*1e3

print(f"Beräknad spaltbredd = {a} mm")

# Uppgift 2:
print("\nUppgift 2:")
I_0 = 110  # Mät med PASCO
I = 7  # Mät med PASCO

print(f"Kvot = {I/I_0}")

# Uppgift 3:
print("\nUppgift 3:")
x_0 = 7.5941e-2
x_1 = 7.6434e-2-x_0  # Avstånd till första minimi i centralbilden
x_2 = 8.2352e-2-x_0  # Avstånd till bildövergång

a = (l*L/x_2)*1e3
d = (l*L/(2*x_1))*1e3

print(f"Spaltbredd = {a} mm\nAvstånd mellan spalter = {d} mm")


