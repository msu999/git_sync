from sympy import *
from sympy import I as Imag_unit
from custom_functions import ei, perp, vector_integral, rotmat as M
from expressions import x, y, t, f, ga_y, la, la_1, la_2, omega, d_y_ga, d_x_ga, f_0, omega_0, la_1_0, la_2_0, phi, delta_1, delta_2, eta, s, la, d_y_phi, diff_ga_la, diff_ga_ga, eps
import sys

sys.path.insert(0, "/home/max/Documents/Python/pickler")  # Adds pickler folder to system path

from pickler import storeData, loadData  # Imports pickler functions



######################################################
# --- Expressions particular to this calculation --- #
######################################################



#################################################
# --- Construction of $X_0^1(R, r, \omega)$ --- #
#################################################

# -- Defining parts of $X_0^1$ -- #
I = (-(1/(4*pi)) * log((la - ga_y).dot(la - ga_y))).expand().simplify()
J = d_y_ga[0]
print(f"I = {latex(I)}\n\nJ = {latex(J)}\n")


# -- Definition of $X_2$ -- #
X = Integral(I*J, (y, 0, 2*pi))



###############################################################
# --- Performing variation of $X_0^1(f, \lambda, \omega)$ --- #
###############################################################

# -- Substituting original functions for their varied counterparts -- #
# Make substitution: $f \mapsto 0 + s*\phi$, $\lambda \mapsto \lambda_0 + s*\delta$ and $\omega \mapsto \omega_0 + s*\eta$
var_subs = {f(x, t) : s*phi(x, t), f(y, t) : s*phi(y, t), la_1(t) : la_1_0 + s*delta_1(t), la_2(t) : la_2_0 + s*delta_2(t), omega : omega_0 + s*eta}

I_subs = I.subs(var_subs)
dI_subs = I_subs.diff(s)
dI_0 = dI_subs.subs(s, 0).expand().simplify()
I_0 = I_subs.subs(s, 0).expand().simplify()

J_subs = J.subs(var_subs)
dJ_subs = J_subs.diff(s)
dJ_0 = dJ_subs.subs(s, 0).expand().simplify().factor()
J_0 = J_subs.subs(s, 0).expand().simplify()

print(f"I_0 = {latex(I_0)}\n\nJ_0 = {latex(J_0)}\n\ndI_0 = {latex(dI_0)}\n\ndJ_0 = {latex(dJ_0)}\n")

dIntegrand_0 = dI_0*J_0 + I_0*dJ_0
print(f"dIntegrand_0 before ibp = {latex(dIntegrand_0)}\n")

dJ_0 = Derivative(cos(y)*phi(y, t), y)/2  # By observation

# Integrating second term by parts:
#
#                      dI_0   /
#   I_0*dJ_0  --->  -  ---- * | dJ_0 dy
#                       dy    /
#
dIntegrand_0 = dI_0*J_0 - I_0.diff(y)*integrate(dJ_0, y) 
print(f"dIntegrand_0 after ibp = {latex(dIntegrand_0)}\n")

#F_0_var = Integral(dIntegrand, (y, 0, 2*pi))
#
##print(f"Complete variation of F_0 = {latex(F_0_var)}", end="\n\n")
#
##print(f"I_0 = {latex(I_0)}\n\ndI_0 = {latex(dI_0)}\n\nJ_0 = {latex(J_0)}\n\ndJ_0 = {latex(dJ_0)}\n\ndJ_0_simp = {latex(dJ_0_simp)}\n\ndB_0 = {latex(dB_0)}")
#
#
#
###########################################
## --- Analysis of parts of variation --- #
###########################################
#
## -- Extracting phi coefficients from integrand -- #
#dIntegrand_N, dIntegrand_D = dIntegrand.expand().factor().as_numer_denom()
#
##print(f"dIntegrand_N = {latex(dIntegrand_N)}", end="\n\n")
#
#phi_x_coeff_int = dIntegrand_N.coeff(phi(x))/dIntegrand_D
#phi_y_coeff_int = dIntegrand_N.coeff(phi(x + y))/dIntegrand_D
#d_phi_x_coeff_int = dIntegrand_N.coeff(phi(x).diff(x))/dIntegrand_D
#d_phi_y_coeff_int = dIntegrand_N.coeff(phi(x + y).diff(y))/dIntegrand_D
#
##print(f"phi(x) integrand coeff = {latex(phi_x_coeff_int)}\n\nphi(x + y) integrand coeff = {latex(phi_y_coeff_int)}\n\nd_phi(x) integrand coeff = {latex(d_phi_x_coeff_int)}\n\nd_phi(x + y) integrand coeff = {latex(d_phi_y_coeff_int)}")
#
## -- Total phi coefficients -- #
#phi_x_coeff = Integral(phi_x_coeff_int, (y, 0, 2*pi))
#d_phi_x_coeff = Integral(d_phi_x_coeff_int, (y, 0, 2*pi))
#
## -- Expression before application of fourier coefficient operator -- #
#A = phi_x_coeff
#B = phi_y_coeff_int
#C = d_phi_x_coeff
#D = d_phi_y_coeff_int
##print(f"A = {latex(A)}\n\nB = {latex(B)}\n\nC = {latex(C)}\n\nD = {latex(D)}\n\n")
#
## -- Fourier coefficient of LHS -- #
#k, phi_hat = symbols(r"k, \hat{\phi}_k")
#LHS_hat = (A + Integral(B*exp(Imag_unit*k*y), (y, 0, 2*pi)) + C*Imag_unit*k + Integral(D*Imag_unit*k*exp(Imag_unit*k*y), (y, 0, 2*pi)))*phi_hat
#
#run_integration = False  # Run integration of LHS_hat? (takes more than a minute)
#
#if run_integration:
#    # -- Saves results necessary for further calculations in db -- #
#    LHS_integrated = LHS_hat.doit().expand().simplify()
#    data = {"LHS_hat_integrated" : LHS_integrated}
#    storeData(data, "G_1_db")
#
#if not run_integration:
#    data = loadData("G_1_db")
#    LHS_integrated = data["LHS_hat_integrated"]
#
#print(latex(LHS_integrated.subs(k, -1)))
