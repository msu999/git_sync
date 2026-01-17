from sympy import *
from propagation_of_error import error

omega_plus, omega_minus = symbols("omega_plus, omega_minus")
plus_error, minus_error = symbols("plus_error, minus_error")
M, g, c, k, b, I_P = symbols("M, g, c, k, b, I_P")

K = (omega_minus**2 - omega_plus**2)/(omega_minus**2 + omega_plus**2)
print(latex(K.subs({omega_plus : sqrt(M*g*c/I_P), omega_minus : sqrt((M*g*c + 2*k*b**2)/I_P)}).expand().simplify()))

X = Matrix([omega_plus, omega_minus])
X_error = Matrix([plus_error, minus_error])

subs_dict = {omega_plus : 4.53, omega_minus : 4.81, plus_error : 0.0002, minus_error : 0.0001}

K_value = K.subs(subs_dict)
K_error = error(K, X, X_error).subs(subs_dict)

print(f"K = {K_value}\n\nK error = {K_error}")
