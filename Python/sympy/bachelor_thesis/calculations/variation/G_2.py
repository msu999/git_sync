from sympy import *
from custom_functions import ei, perp, vector_integral, rotmat as M
from expressions import x, y, t, f, h, Omega, d_y_ga, f_0, h_0, Omega_0, phi, delta, eta, s, la, d_y_phi, diff_ga_la
import sys

sys.path.insert(0, "/home/max/Documents/Python/pickler")  # Adds pickler folder to system path

from pickler import storeData, loadData  # Imports pickler functions



######################################################
# --- Expressions particular to this calculation --- #
######################################################



#############################################
# --- Construction of $G(R, r, \omega)$ --- #
#############################################

# -- Defining parts of $G_2$ -- #
A = 1/(4*pi)
I = log((1 + f(y)) + h**2 - 2*h*sqrt(1 + f(y))*sin(y))
J = d_y_ga[1]
B = (Omega*M(pi/2)*la)[1]
print(f"I = {latex(I)}\n\nJ = {latex(J)}\n\nB = {latex(B)}", end="\n\n")


# -- Definition of $G_2$ -- #
G_2 = A*Integral(I*J, (y, 0, 2*pi)) + B
print(latex(G_2), end="\n\n")



######################################################
## --- Performing variation of $G_2(f, h, \Omega)$ --- #
######################################################

# -- Substituting original functions for their varied counterparts -- #
# Make substitution: $f \mapsto 0 + s*\phi$, $h \mapsto h_0 + s*\delta$ and $\Omega \mapsto \Omega_0 + s*\eta$
var_subs = {f(x) : s*phi(x), f(y) : s*phi(y), h : h_0 + s*delta, Omega : Omega_0 + s*eta}

I_subs = I.subs(var_subs)
dI_subs = I_subs.diff(s)
dI_0 = dI_subs.subs(s, 0).expand().simplify()
I_0 = I_subs.subs(s, 0).expand().simplify()

J_subs = J.subs(var_subs)
dJ_subs = J_subs.diff(s)
dJ_0 = dJ_subs.subs(s, 0).expand().simplify()
J_0 = J_subs.subs(s, 0).expand().simplify()

B_subs = B.subs(var_subs)
dB_subs = B_subs.diff(s)
dB_0 = dB_subs.subs(s, 0).expand().simplify()

G_2_var = A*Integral(dI_0*J_0 + I_0*dJ_0, (y, 0, 2*pi)) + dB_0

print(f"Complete variation of G_2 = {latex(G_2_var)}", end="\n\n")

print(f"I_0 = {latex(I_0)}\n\ndI_0 = {latex(dI_0)}\n\nJ_0 = {latex(J_0)}\n\ndJ_0 = {latex(dJ_0)}\n\ndB_0 = {latex(dB_0)}")
