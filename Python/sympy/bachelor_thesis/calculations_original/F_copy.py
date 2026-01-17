from sympy import *
from custom_functions import ei, perp, vector_integral
from expressions import alpha, beta, t, R, diff_ga_ga, diff_ga_la, d_alpha_ga, d_beta_ga, eps, R_0, Delta, s 

# Construction of $F(R, r, \omega)$
I = log((diff_ga_ga).norm())
J = d_beta_ga.dot(ei(alpha))

epsilon_term = eps*perp(diff_ga_la).dot(perp(ei(alpha)))/(diff_ga_la.norm()**2)

# $F = -frac{1}{2\pi} \cdot \int_0^{2\pi} IJ d\beta - epsilon_term$
F = -1/(2*pi) * (Integral(I*J, (beta, 0, 2*pi)) - epsilon_term)

# Print out I and epsilon_term
print(f"I = {latex(I.expand().simplify())}\nJ = {latex(J.expand().simplify())}\nepsilon term = {latex(epsilon_term.expand().simplify())}")

# F with $R \mapsto R_0 + \Delta$ (F_var)
I_var = I.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t)})
J_var = J.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t)}).expand().simplify()

# F_var = F.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t), r(t) : r_0 + delta(t), omega(t) : omega_0*t + eta(t)})

# Print out F_var
# print(latex(F_var))

# Calculating variation of F
I_0 = I_var.subs(s, 0).expand().simplify()
J_0 = J_var.subs(s, 0).expand().simplify()

dI_0 = I_var.diff(s).subs(s, 0).expand().simplify()
dJ_0 = J_var.diff(s).subs(s, 0).expand().simplify()

dIntegrand_0 = (dI_0*J_0 + I_0*dJ_0).simplify()
#print(f"Integrand I*J = {latex(dIntegrand_0)}\n\n")

#print(f"I(s) = {latex(I_var)}\n\nJ(s) = {latex(J_var)}\n\nI_0 = {latex(I_0)}\n\nJ_0 = {latex(J_0)}\n\ndI_0 = {latex(dI_0)}\n\ndJ_0 = {latex(dJ_0)}")

