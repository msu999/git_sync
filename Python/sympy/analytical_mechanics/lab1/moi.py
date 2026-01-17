import numpy as np
from sympy import *
from propagation_of_error import error
from spring_constant import k as k_value, k_error as k_error_value

# Symbols & quantities
g = 9.82
m_s, m_v, a, b, d, h, k = symbols("m_s, m_v, a, b, d, h, k") 
eps_mass, eps_dist, eps_k = symbols("eps_mass, eps_dist, eps_k")

# Expressions derived from symbols
A, B, C = a + d - h, h, a

lambda_s = m_s/(A + B)
lambda_v = m_v/B

lambda_A = lambda_s
lambda_B = lambda_s + lambda_v

M = lambda_A*A + lambda_B*B

x_cm = (1/(2*M)) * (lambda_A*A**2 + lambda_B*(B**2 + 2*A*B))

c = x_cm - a

I_0 = ((lambda_A/3) * A**3) + ((lambda_B/3) * (B**3 + 3*A*B**2 + 3*B*A**2))
I_cm = I_0 - M*x_cm**2
I_P = M*(C**2 - 2*C*x_cm) + I_0

omega_plus = sqrt(M*g*c/I_P)
omega_minus = sqrt((M*g*c + 2*k*b**2)/I_P)

K = (omega_minus**2 - omega_plus**2)/(omega_minus**2 + omega_plus**2)

# Defining symbolic vectors for calculation of error
X = Matrix([m_s, m_v, a, b, d, h, k])
error_X = Matrix([eps_mass, eps_mass, eps_dist, 1e-2, eps_dist, eps_dist, eps_k])

# Dictionary with substitutions of symbols with their true values
subs_dict = {m_s : 32.15e-3, m_v : 70.41e-3, a : 2.6e-2, b : 10e-2, d : 53.4e-2, h : 2.0e-2, k : k_value, eps_mass : 0.005e-3, eps_dist : 0.05e-2, eps_k : k_error_value}

c_value = c.subs(subs_dict)
I_P_value = I_P.subs(subs_dict)
omega_plus_value = omega_plus.subs(subs_dict)
omega_minus_value = omega_minus.subs(subs_dict)
K_value = K.subs(subs_dict)

I_P_error = error(I_P, X, error_X)
c_error = error(c, X, error_X)
omega_plus_error = error(omega_plus, X, error_X)
omega_minus_error = error(omega_minus, X, error_X)
K_error = error(K, X, error_X)

#print(latex(c_error))
print(f"c = {c_value}")
print("c error = " + latex(c_error.subs(subs_dict)))

print(latex(I_P_error))
print(f"I_P = {I_P_value}")
print("I_P error = " + latex(I_P_error.subs(subs_dict)))

#print(latex(omega_plus_error))
print(f"omega_plus = {omega_plus_value}")
print("omega_plus error = " + latex(omega_plus_error.subs(subs_dict)))

#print(latex(omega_minus_error))
print(f"omega_minus = {omega_minus_value}")
print("omega_minus error = " + latex(omega_minus_error.subs(subs_dict)))

#print(latex(K_error))
print(f"K = {K_value}")
print("K error = " + latex(K_error.subs(subs_dict)))

# print(f"x_cm = {x_cm}, c = {c}, I_0 = {I_0}, I_cm = {I_cm}, I_P = {I_P}")

# True values
# m_s, m_v = 32.15e-3, 70.41e-3
# a, d, h = 2.6e-2, 53.4e-2, 2.0e-2

