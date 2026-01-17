from sympy import *

#
#  Calculates propagation of error of a expression (function) f(X) where X is
#  a vector (X[1], ..., X[n]) where X[i] has an error of error_X[i].
#

def error(f, X, error_X):
    n = X.shape[0]
    
    result_squared = sympify(0)
    
    for i in range(n):
        result_squared += (f.diff(X[i])*error_X[i])**2

    return sqrt(result_squared)


# Example usage:

#x, y, eps_x, eps_y = symbols(r"x, y, \varepsilon_x, \varepsilon_y")
#X = Matrix([x, y])
#error_X = Matrix([eps_x, eps_y])
#
#f = x**2 + sin(y)
#
#print(latex(error(f, X, error_X)))
