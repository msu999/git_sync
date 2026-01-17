from sympy import *
import numpy as np
from scipy.constants import g

pi = np.pi
g = g
d = 560
R = d

omega_g = np.sqrt(g/R)

M = 68000e3 

A = (((d**2)/2 + (R**2))**2)/(R**2)
B = (pi*(d**4)*(omega_g**2))/(2*(R**2))

a_max = 0.02

alpha_max = -B/(2*A) + np.sqrt((B**2)/(4*A**2) + (a_max**2)/A)

T = np.sqrt(2*np.pi/alpha_max)

print(omega_g)
print(omega_g/(2*np.pi))
print(alpha_max)
