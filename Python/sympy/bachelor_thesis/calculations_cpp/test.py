from sympy import *
from expressions import y
import sys

sys.path.insert(0, "/home/max/Documents/Python/pickler")  # Adds pickler folder to system path

from pickler import storeData, loadData  # Imports pickler functions

expr = loadData("G_1_db")["LHS_hat_integrated"]

ibp_integrand = expr.args[2].args[2].args[0]
ibp_new_integrand = -I*k*ibp_integrand.subs(sin(y), 1 - cos(y)) + ()
ibp_subs = {ibp_integrand : exp(I*k*y)}

print(expr)
