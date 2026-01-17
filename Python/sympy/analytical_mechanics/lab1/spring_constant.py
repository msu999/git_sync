import numpy as np
from sympy import *
from propagation_of_error import error

g = 9.82

# For calculating spring constant
m_wh = 19.52e-3  # Weight  mass
m_w = 19.33e-3  # Weight

m_ws = np.array([m_wh, m_wh + m_w, m_wh + 2*m_w])

A = (60.5-44.5)*1e-2  # Spring elongation under its own weight
Bs = np.array([(67-44.5)*1e-2, (73.7-44.5)*1e-2, (80.5-44.5)*1e-2])  # Spring elongation under ()

ks = m_ws*g/(Bs - A)
k = np.mean(ks)
print(f"Spring constant = {k}")

M_W0, M_W1, M_W2 = symbols("M_W0, M_W1, M_W2")
eps_M = symbols("eps_M")
B_0, B_1, B_2 = symbols("B_0, B_1, B_2")
eps_B = symbols("eps_B")
K_1, K_2, K_3 = M_W0*g/(B_0 - A), M_W2*g/(B_2 - A), M_W1*g/(B_1 - A)
K = (1/3)*(K_1 + K_2 + K_3)
k_error = error(K, Matrix([M_W0, M_W1, M_W2, B_0, B_1, B_2]), Matrix([eps_M, eps_M, eps_M, eps_B, eps_B, eps_B])).subs({M_W0 : m_ws[0], M_W1 : m_ws[1], M_W2 : m_ws[2], B_0 : Bs[0], B_1 : Bs[1], B_2 : Bs[2], eps_M : 0.005e-3, eps_B : 0.05e-2})

print(latex(k_error))
