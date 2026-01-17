from sympy import *
from custom_functions import ei, perp, vector_integral
from expressions import alpha, beta, t, R, r, omega, ga_beta, d_beta_ga, R_0, r_0, omega_0, Delta, delta, eta, s, la



######################################################
# --- Expressions particular to this calculation --- #
######################################################

diff_ga_la = ga_beta - la  # $\gamma(\beta, t) - \lambda(t)$



#############################################
# --- Construction of $G(R, r, \omega)$ --- #
#############################################

# -- Defining parts of $G$ -- #
I = log((diff_ga_la).norm()).expand().simplify()
J = d_beta_ga.dot(ei(omega(t))).expand().simplify()
#print(f"I = {latex(I)}\n\nJ = {latex(J)}")


# -- Performing integration by parts described in document -- #
# $\dot r$ before integration by parts below
r_dot_before_ibp = -1/(2*pi) * Integral(I*J, (beta, 0, 2*pi))#.expand().simplify()

# Simplification of $\dot r$ through integration by parts
I_diff = I.diff(beta)
J_int = R(beta, t)*cos(beta - omega(t))
new_integrand = (I_diff*J_int).expand().simplify()

# $\dot r$ after integration by parts
r_dot = 1/(2*pi) * Integral(new_integrand, (beta, 0, 2*pi))#.expand().simplify()
# print(latex(r_dot))


# -- Definition of $G$ -- #
G = r_dot - r.diff(t)
# print(latex(G))



#####################################################
# --- Performing variation of $G(R, r, \omega)$ --- #
#####################################################

# -- Substituting original functions for their varied counterparts -- #
# G with $R \mapsto R_0 + s*\Delta$, $r \mapsto r_0 + s*\delta$ and $\omega \mapsto \omega_0 + s*\eta$ (G_var)
I_var = I.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t), r(t) : r_0 + s*delta(t), omega(t) : omega_0*t + s*eta(t)}).expand().simplify()
J_var = J.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t), r(t) : r_0 + s*delta(t), omega(t) : omega_0*t + s*eta(t)})
# print(latex(J_var.expand().simplify().factor()))

G_var = G.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t), r(t) : r_0 + s*delta(t), omega(t) : omega_0*t + s*eta(t)})
# print(latex(G_var))


# -- The calculation -- #
I_0 = I_var.subs(s, 0).expand().simplify()
J_0 = J_var.subs(s, 0).expand().simplify()

dI_0 = I_var.diff(s).subs(s, 0).expand().simplify()
dJ_0 = J_var.diff(s).subs(s, 0).expand().simplify()

dIntegrand_0 = (dI_0*J_0 + I_0*dJ_0).simplify()
# print(latex(dIntegrand_0))


# -- Printing of result -- #
print(f"I(s) = {latex(I_var)}\n\nJ(s) = {latex(J_var)}\n\nI_0 = {latex(I_0)}\n\nJ_0 = {latex(J_0)}\n\ndI_0 = {latex(dI_0)}\n\ndJ_0 = {latex(dJ_0)}")

