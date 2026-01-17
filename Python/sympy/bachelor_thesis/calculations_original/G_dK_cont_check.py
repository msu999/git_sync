# General imports
from sympy import *

# Import of often used sympy expressions
from expressions import s, beta, eta, t, u, R_0, r_0, Delta, delta

# Import of pickler
import sys
sys.path.insert(0, "/home/max/Documents/Python/pickler")  # Adds pickler folder to system path
from pickler import storeData, loadData  # Imports pickler functions



########################################
# --- Loads calculations from db:s --- #
########################################
# --- Loads relevant db:s --- #
G_db = loadData("G_db")

# --- Loads calculations from said db:s --- #
dK_var = G_db["dK_var"]
#print(latex(dK_var))



################################
# --- Simplification of K' --- #
################################

# --- Args and subargs of dK_var --- #
if False:
    for i, arg in enumerate(dK_var.args):
        print(f"Main arg {i}: {latex(arg)}", end="\n\n")
        
        for j, subarg in enumerate(arg.args):
            print(f"Subarg {i}.{j}: {latex(subarg)}", end="\n\n")

# Notes from above print:
#
#   Term 0 (arg 0) consists of numerator subarg 0.1 * subarg 0.2, and denom subarg 0.0
#
#   Term 1 consists of numerator 1.1*1.2*1.3
#
#   Term 2 consists of numerator -2.2*2.3*2.4

# --- Denominator of terms --- ##
denom = dK_var.args[0].args[0].as_numer_denom()[1]
#print(latex(denom))

denom_cos_coeff = denom.coeff(cos(-beta + s*eta(t) + t*u/r_0))  # Will abbreviate as dcc
dcc_roots = solve(denom_cos_coeff, s)
#print(latex(().expand().simplify()))
dcc_simplified = ((s - dcc_roots[0]) * Delta(beta, t)).expand() * ((s - dcc_roots[1]) * delta(t)).expand() * (-2)
denom_cos_term = dcc_simplified * cos(-beta + s*eta(t) + t*u/r_0)

denom_rest = (denom - denom_cos_term).expand().simplify()
rest_roots = solve(denom_rest, s)
#print(latex(rest_roots), end="\n\n")
rest_simplified = ((s - rest_roots[0]) * (Delta(beta, t) + I*delta(t))).expand().simplify() * ((s - rest_roots[1]) * (Delta(beta, t) - I*delta(t))).expand().simplify()
rest_simplified2 = (R_0 + s*Delta(beta, t))**2 + (r_0 + s*delta(t))**2
#print(latex(rest_simplified2))

denom_simplified = rest_simplified2 + denom_cos_term

control_rest = (denom - denom_simplified).expand().simplify()  # Should evaluate to 0

print(solve(denom, s))

print(f"Denom = {latex(denom)}\n\nDenom cos term = {latex(denom_cos_term)}\n\nDenom rest = {latex(rest_simplified2)}\n\nDenom simplified = {latex(denom_simplified)}\n\nControl rest = {latex(control_rest)}")

#print(latex(dcc_simplified))
# dcc is non-zero iff $R_0 + s\Delta(\beta, t) \neq 0$ and $r_0 + s\delta(t) \neq 0$, if $s \neq -R_0/\Delta$ and $s \neq -r_0/\delta$ but this is true if $s$ is close to $0$ since $R_0/\Delta, r_0/\delta > 0$ for all $\beta, t$. Hence $dK_var$ is continuous 
