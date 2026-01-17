from sympy import *
from custom_functions import ei, perp, vector_integral
from expressions import alpha, beta, t, R, diff_ga_ga, diff_ga_la, d_alpha_ga, d_beta_ga, eps, R_0, Delta, s 

# Construction of $F(R, r, \omega)$
A = 1/(2*pi*R(alpha, t))
I = log((diff_ga_ga).norm())
J = d_beta_ga.dot(perp(d_alpha_ga))

epsilon_term = eps*perp(diff_ga_la).dot(perp(d_alpha_ga))/(diff_ga_la.norm()**2)

# $F = A \cdot \int_0^{2\pi} IJ d\beta - epsilon_term$
F = A * (Integral(I*J, (beta, 0, 2*pi)) - epsilon_term)

# Print out I and epsilon_term
# print(f"I = {latex(I.expand().simplify())}\nepsilon term = {latex(epsilon_term.expand().simplify())}")

# F with $R \mapsto R_0 + \Delta$ (F_var)
A_var = A.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t)}) 
I_var = I.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t)})
J_var = J.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t)})

F_var = F.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t)})

# Print out F_var
# print(latex(F_var))

# Calculating variation of F
A_0 = A_var.subs(s, 0)
I_0 = I_var.subs(s, 0).expand().simplify()
J_0 = J_var.subs(s, 0).expand().simplify()

dA_0 = A_var.diff(s).subs(s, 0)
dI_0 = I_var.diff(s).subs(s, 0).expand().simplify()
dJ_0 = J_var.diff(s).subs(s, 0).expand().simplify()

dIntegrand_0 = (dI_0*J_0 + I_0*dJ_0).simplify()
# print(latex(dIntegrand_0))

#print(f"A(s) = {latex(A_var)}\nI(s) = {latex(I_var)}\nJ(s) = {latex(J_var)}\nA_0 = {latex(A_0)}\nI_0 = {latex(I_0)}\nJ_0 = {latex(J_0)}\ndA_0 = {latex(dA_0)}\ndI_0 = {latex(dI_0)}\ndJ_0 = {latex(dJ_0)}")

