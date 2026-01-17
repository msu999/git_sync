from sympy import *
from sympy import I as Imag_unit
from custom_functions import ei, perp, vector_integral, rotmat as M
from expressions import x, y, t, f, h, Omega, d_y_ga, d_x_ga, f_0, h_0, Omega_0, phi, delta, eta, s, la, d_y_phi, diff_ga_la, diff_ga_ga, eps
import sys

sys.path.insert(0, "/home/max/Documents/Python/pickler")  # Adds pickler folder to system path

from pickler import storeData, loadData  # Imports pickler functions



######################################################
# --- Expressions particular to this calculation --- #
######################################################



#############################################
# --- Construction of $G(R, r, \omega)$ --- #
#############################################

# -- Defining parts of $F$ -- #
A = 1/(4*pi)
I = log(2 + f(x) + f(y) - 2*sqrt(1 + f(x))*sqrt(1 + f(y))*cos(x - y))
J = d_y_ga.dot(perp(d_x_ga)).expand().simplify()
B = Omega*f(x).diff(x)
#print(f"I = {latex(I)}\n\nJ = {latex(J)}\n\nB = {latex(B)}")


# -- Definition of $F$ -- #
F_0 = A*Integral(I*J, (y, 0, 2*pi)) + B
F_tilde = -(eps/(2*pi)) * diff_ga_la.dot(d_x_ga)/diff_ga_la.norm()
F = F_0 + eps*F_tilde
#print(latex(F))



#####################################################
# --- Performing variation of $F(f, h, \Omega)$ --- #
#####################################################

# -- Substituting original functions for their varied counterparts -- #
# Make substitution: $f \mapsto 0 + s*\phi$, $h \mapsto h_0 + s*\delta$ and $\Omega \mapsto \Omega_0 + s*\eta$
var_subs = {f(x) : s*phi(x), f(y) : s*phi(y), h : h_0 + s*delta, Omega : Omega_0 + s*eta}

I_subs = I.subs(var_subs)
dI_subs = I_subs.diff(s)
dI_0 = dI_subs.subs(s, 0).expand().simplify()
I_0 = I_subs.subs(s, 0).expand().simplify()

J_subs = J.subs(var_subs)
dJ_subs = J_subs.diff(s)
dJ_0 = dJ_subs.subs(s, 0).expand().simplify().factor()
dJ_0_N, dJ_0_D = dJ_0.as_numer_denom()
sin_coeff = dJ_0_N.coeff(sin(x - y))
cos_coeff = dJ_0_N.coeff(cos(x - y))
dJ_0_simp = (1/dJ_0_D) * (sin_coeff*sin(x - y) + cos_coeff*cos(x - y))
J_0 = J_subs.subs(s, 0).expand().simplify()

B_subs = B.subs(var_subs)
dB_subs = B_subs.diff(s)
dB_0 = dB_subs.subs(s, 0).expand().simplify()

F_tilde_subs = F_tilde.subs(var_subs)
dF_tilde_subs = F_tilde_subs.diff(s)
dF_tilde_0 = dF_tilde_subs.subs(s, 0).expand().simplify()

dIntegrand = dI_0*J_0 + I_0*dJ_0
dIntegrand = dIntegrand.subs(y, x + y).simplify()  # Perform change of variables $y \mapsto x + y$

F_0_var = A*Integral(dIntegrand, (y, 0, 2*pi)) + dB_0

#print(f"Complete variation of F_0 = {latex(F_0_var)}", end="\n\n")

#print(f"I_0 = {latex(I_0)}\n\ndI_0 = {latex(dI_0)}\n\nJ_0 = {latex(J_0)}\n\ndJ_0 = {latex(dJ_0)}\n\ndJ_0_simp = {latex(dJ_0_simp)}\n\ndB_0 = {latex(dB_0)}")



##########################################
# --- Analysis of parts of variation --- #
##########################################

# -- Extracting phi coefficients from integrand -- #
dIntegrand_N, dIntegrand_D = dIntegrand.expand().factor().as_numer_denom()

#print(f"dIntegrand_N = {latex(dIntegrand_N)}", end="\n\n")

phi_x_coeff_int = A*dIntegrand_N.coeff(phi(x))/dIntegrand_D
phi_y_coeff_int = A*dIntegrand_N.coeff(phi(x + y))/dIntegrand_D
d_phi_x_coeff_int = A*dIntegrand_N.coeff(phi(x).diff(x))/dIntegrand_D
d_phi_y_coeff_int = A*dIntegrand_N.coeff(phi(x + y).diff(y))/dIntegrand_D

# -- Extracting phi coefficients from dB_0 -- #
phi_x_coeff_dB = dB_0.coeff(phi(x))
d_phi_x_coeff_dB = dB_0.coeff(phi(x).diff(x))

#print(f"phi(x) integrand coeff = {latex(phi_x_coeff_int)}\n\nphi(x + y) integrand coeff = {latex(phi_y_coeff_int)}\n\nd_phi(x) integrand coeff = {latex(d_phi_x_coeff_int)}\n\nd_phi(x + y) integrand coeff = {latex(d_phi_y_coeff_int)}")

# -- Total phi coefficients -- #
phi_x_coeff = Integral(phi_x_coeff_int, (y, 0, 2*pi)) + phi_x_coeff_dB
d_phi_x_coeff = Integral(d_phi_x_coeff_int, (y, 0, 2*pi)) + d_phi_x_coeff_dB

# -- Expression before application of fourier coefficient operator -- #
A = phi_x_coeff
B = phi_y_coeff_int
C = d_phi_x_coeff
D = d_phi_y_coeff_int
#print(f"A = {latex(A)}\n\nB = {latex(B)}\n\nC = {latex(C)}\n\nD = {latex(D)}\n\n")

# -- Fourier coefficient of LHS -- #
k, phi_hat = symbols(r"k, \hat{\phi}_k")
LHS_hat = (A + Integral(B*exp(Imag_unit*k*y), (y, 0, 2*pi)) + C*Imag_unit*k + Integral(D*Imag_unit*k*exp(Imag_unit*k*y), (y, 0, 2*pi)))*phi_hat

run_integration = False  # Run integration of LHS_hat? (takes more than a minute)

if run_integration:
    # -- Saves results necessary for further calculations in db -- #
    LHS_integrated = LHS_hat.doit().expand().simplify()
    data = {"LHS_hat_integrated" : LHS_integrated}
    storeData(data, "G_1_db")

if not run_integration:
    data = loadData("G_1_db")
    LHS_integrated = data["LHS_hat_integrated"]

print(latex(LHS_integrated.subs(k, -1)))
