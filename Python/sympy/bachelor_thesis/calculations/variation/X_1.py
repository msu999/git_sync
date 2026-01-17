from sympy import *
from sympy import I as Imag_unit
from custom_functions import ei, perp, vector_integral, rotmat as M
from expressions import x, y, t, f, la, la_1, la_2, omega, d_y_ga, d_x_ga, f_0, omega_0, la_1_0, la_2_0, phi, delta_1, delta_2, eta, s, la, d_y_phi, diff_ga_la, diff_ga_ga, eps
import sys

sys.path.insert(0, "/home/max/Documents/Python/pickler")  # Adds pickler folder to system path

from pickler import storeData, loadData  # Imports pickler functions



######################################################
# --- Expressions particular to this calculation --- #
######################################################



#################################################
# --- Construction of $X_{0,1}(R, r, \omega)$ --- #
#################################################

# -- Defining parts of $X_{0,1}$ -- #
I = (1/pi) * log(2 + f(x, t) + f(y, t) - 2*sqrt(1 + f(x, t))*sqrt(1 + f(y, t))*cos(x - y))
J = d_y_ga.dot(perp(d_x_ga)).expand().simplify()
print(f"I = {python(I)}\n\nJ = {latex(J)}\n\n")


# -- Definition of $X_0^1$ -- #
X = Integral(I*J, (y, 0, 2*pi))
X_tilde = -(1/pi) * diff_ga_la.subs(y, x).dot(d_x_ga)/(diff_ga_la.norm().subs(y, x)**2)
X_total = X + eps*X_tilde
print(f"X_tilde = {latex(X_tilde)}\n\n")
#print(f"Esto = {latex(diff_ga_la)}")


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
dJ_0_N, dJ_0_D = dJ_0.as_numer_denom()
sin_coeff = dJ_0_N.coeff(sin(x - y))
cos_coeff = dJ_0_N.coeff(cos(x - y))
dJ_0 = (1/dJ_0_D) * (sin_coeff*sin(x - y) + cos_coeff*cos(x - y))
J_0 = J_subs.subs(s, 0).expand().simplify()

print(f"I(s) = {latex(I_subs)}\n\ndI(s) = {latex(dI_subs)}\n\nJ(s) = {latex(J_subs)}\n\ndJ(s) = {latex(dJ_subs)}\n")
print(f"I_0 = {latex(I_0)}\n\nJ_0 = {latex(J_0)}\n\ndI_0 = {latex(dI_0)}\n\ndJ_0 = {latex(dJ_0)}\n")

##F_tilde_subs = F_tilde.subs(var_subs)
##dF_tilde_subs = F_tilde_subs.diff(s)
##dF_tilde_0 = dF_tilde_subs.subs(s, 0).expand().simplify()

dIntegrand_0 = dI_0*J_0 + I_0*dJ_0
#print(latex(dIntegrand_0))
dIntegrand_0 = dIntegrand_0.subs(y, x + y).simplify()  # Perform change of variables $y \mapsto x + y$
print(f"dIntegrand_0 = {latex(dIntegrand_0)}\n")


Integrand_0 = I_0*J_0
X_tilde_0 = X_tilde.subs(var_subs).subs(s, 0).expand().simplify()
print(f"Integrand_0 = {latex(Integrand_0)} = 0\n\nX_tilde_0 = {latex(X_tilde_0)}\n")


# -- Integrand analysis -- #
N, D = dIntegrand_0.as_numer_denom()
N = N.expand()
A_1 = N.coeff(phi(x, t))
A_1_ = N.coeff(phi(x + y, t))
A_2 = N.coeff(phi(x, t).diff(x))
A_2_ = N.coeff(phi(x + y, t).diff(y))

print(f"A_1 = {latex(A_1)}\n\nA_2 = {latex(A_2)}\n\nA_1_ = {latex(A_1_)}\n\nA_2_ = {latex(A_2_)}\n")

X_1_db = {"A_1" : A_1, "A_2" : A_2, "X_1" : None}

#storeData(X_1_db, "X_1_db")

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
