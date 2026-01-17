import numpy as np
from sympy import *

hbar = symbols("hbar", positive=True)

a_1, a_2, a_3 = symbols("a_1, a_2, a_3", real=True)
a = Matrix([a_1, a_2, a_3])

b = symbols("b", real=True)

c_1, c_2, c_3 = symbols("c_1, c_2, c_3", real=True)
c = Matrix([c_1, c_2, c_3])

N = (2*sqrt(2)*b**3)/pi**(3/4)

x_1, x_2, x_3 = symbols("x_1, x_2, x_3", real=True)
x = Matrix([x_1, x_2, x_3])

Psi = N*exp(I*a.dot(x) - (b**2)*(x - c).dot(x - c))

def p_hat_1(expr):
    return -I*hbar*expr.diff(x_1)

integrand = conjugate(Psi)*p_hat_1(p_hat_1(Psi)).expand().simplify() 

print(f"Integrand before shift = {latex(integrand)}")

shift_subs = {x_1 : x_1 + c_1, x_2 : x_2 + c_2, x_3 : x_3 + c_3}

integrand = integrand.subs(shift_subs).expand().simplify()

print(f"Integrand after shift = {latex(integrand)}")

bra_p_x_sq_ket = Integral(Integral(Integral(integrand, (x_1, -oo, oo)), (x_2, -oo, oo)), (x_3, -oo, oo))

print(f"bra_p_x_sq_ket = {latex(bra_p_x_sq_ket)}")
