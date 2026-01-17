from sympy import *

x = symbols("x", real = True)
k = symbols("k", positive = True)

print(latex(integrate(sqrt(pi*cos(x)**2 + sin(x)**2), (x, 0, 2*pi)).doit()))
