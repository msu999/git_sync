from sympy import *
import numpy as np

x = symbols("x")

exp = np.linspace(1, 6, 6, dtype=int)
p = 0

for k in exp:
    p += x**k

p_sq = p**2

print((1/36)*p_sq.expand())
