from sympy import *
from custom_functions import ei, perp, vector_integral
from expressions import alpha, beta, t, R, r, omega, ga_beta, d_beta_ga, R_0, r_0, omega_0, Delta, delta, eta, s, la

diff_ga_la = ga_beta - la

# Construction of $G(R, r, \omega)$
A = -1/(2*pi*r(t))
I = log((diff_ga_la).norm())
J = d_beta_ga.dot(perp(ei(omega(t))))

G = A*Integral(I*J, (beta, 0, 2*pi))#.expand().simplify()

# Print out G
# print(latex(G))

# G with $R \mapsto R_0 + s*\Delta$, $r \mapsto r_0 + s*\delta$ and $\omega \mapsto \omega_0 + s*\eta$ (G_var)
A_var = A.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t), r(t) : r_0 + s*delta(t), omega(t) : omega_0*t + s*eta(t)}) 
I_var = I.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t), r(t) : r_0 + s*delta(t), omega(t) : omega_0*t + s*eta(t)})
J_var = J.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t), r(t) : r_0 + s*delta(t), omega(t) : omega_0*t + s*eta(t)})

G_var = G.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t), r(t) : r_0 + s*delta(t), omega(t) : omega_0*t + s*eta(t)})

# Print out G_var
# print(latex(G_var))

# Calculating variation of G
A_0 = A_var.subs(s, 0)
I_0 = I_var.subs(s, 0).expand().simplify()
J_0 = J_var.subs(s, 0).expand().simplify()

dA_0 = A_var.diff(s).subs(s, 0)
dI_0 = I_var.diff(s).subs(s, 0).expand().simplify()
dJ_0 = J_var.diff(s).subs(s, 0).expand().simplify()

#dIntegrand_0 = (dI_0*J_0 + I_0*dJ_0).simplify()
#print(latex(dIntegrand_0))

print(f"A(s) = {latex(A_var)}\nI(s) = {latex(I_var)}\nJ(s) = {latex(J_var)}\nA_0 = {latex(A_0)}\nI_0 = {latex(I_0)}\nJ_0 = {latex(J_0)}\ndA_0 = {latex(dA_0)}\ndI_0 = {latex(dI_0)}\ndJ_0 = {latex(dJ_0)}")

