from sympy import *
from custom_functions import ei, perp, vector_integral
from expressions import alpha, beta, t, R, r, omega, ga_beta, d_beta_ga, R_0, r_0, omega_0, Delta, delta, eta, s, la, u, theta, d_beta_Delta

diff_ga_la = ga_beta - la

# Construction of $H(R, r, \omega)$
A = -1/(2*pi*r(t))
I = log((diff_ga_la).norm()).expand().simplify()
J = d_beta_ga.dot(perp(ei(omega(t)))).expand().simplify()

# print(f"I = {latex(I)}\n\nJ = {latex(J)}")

# $\dot\omega$ before integration by parts
omega_dot_before_ibp = A*Integral(I*J, (beta, 0, 2*pi))#.expand().simplify()

B = -A  # New coefficient in front of integral from integration by parts
I_diff = I.diff(beta)
J_int = R(beta, t)*sin(beta - omega(t))
K = (I_diff*J_int).expand().simplify()  # New integrand obtained from integration by parts

# $\dot\omega$ after integration by parts (vanishing boundary term due to periodicity)
omega_dot = -A*Integral(K, (beta, 0, 2*pi))#.expand().simplify()
#print(f"\\dot\\omega = {latex(omega_dot)}")

# -- Definition of $H$ -- #
H = omega_dot - omega.diff(t)
# print(latex(H))



#####################################################
# --- Performing variation of $H(R, r, \omega)$ --- #
#####################################################

# -- Substituting original functions for their varied counterparts -- #
# Make substitution: $R \mapsto R_0 + s*\Delta$, $r \mapsto r_0 + s*\delta$ and $\omega \mapsto \omega_0 + s*\eta$ (H_var)
B_var = B.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t), r(t) : r_0 + s*delta(t), omega(t) : omega_0*t + s*eta(t)}) 
K_var = K.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t), r(t) : r_0 + s*delta(t), omega(t) : omega_0*t + s*eta(t)})

H_var = H.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t), r(t) : r_0 + s*delta(t), omega(t) : omega_0*t + s*eta(t)})

# Print out H_var
# print(latex(H_var))


# -- Calculating variation of H -- #
B_0 = B_var.subs(s, 0)
K_0 = K_var.subs(s, 0).expand().simplify()

dB_0 = B_var.diff(s).subs(s, 0)
dK_0 = K_var.diff(s).subs(s, 0).expand().simplify()

# Print variation results
#print(f"B(s) = {latex(B_var)}\n\nK(s) = {latex(K_var)}\n\nB_0 = {latex(B_0)}\n\nK_0 = {latex(K_0)}\n\ndB_0 = {latex(dB_0)}\n\ndK_0 = {latex(dK_0)}\n\n")

var_integrand = dB_0*K_0 + B_0*dK_0

# -- Simplifying result -- #
theta_subs = {beta - t*u/r_0 : theta, 2*beta - 2*t*u/r_0 : 2*theta}
N, D = var_integrand.subs(theta_subs).as_numer_denom()
N = N.expand()

#print(f"N = {latex(N)}\n\n")

delta_coeff_N = N.coeff(delta(t)).simplify()
eta_coeff_N = N.coeff(eta(t)).simplify()

rest_N = (N - delta_coeff_N*delta(t) - eta_coeff_N*eta(t)).expand().simplify().expand()

Delta_coeff_N = rest_N.coeff(Delta(beta, t)).simplify()
d_beta_Delta_coeff_N = rest_N.coeff(d_beta_Delta).simplify()

control_rest = (rest_N - Delta_coeff_N*Delta(beta, t) - d_beta_Delta_coeff_N*d_beta_Delta).expand().simplify()

delta_coeff = (delta_coeff_N/D).simplify()
eta_coeff = (eta_coeff_N/D).simplify()
rest = (rest_N/D).simplify()
Delta_coeff = (Delta_coeff_N/D).simplify()
d_beta_Delta_coeff = (d_beta_Delta_coeff_N/D).simplify()

# -- Printing of results -- #
#print(f"K(s) = {latex(K_var)}\n\ndK_0 = {latex(dK_0)}")
print(f"delta numerator coefficient = {latex(delta_coeff_N)}\n\neta numerator coefficient = {latex(eta_coeff_N)}\n\nDelta numerator coefficient = {latex(Delta_coeff_N)}\n\nd_beta_Delta numerator coefficient = {latex(d_beta_Delta_coeff_N)}\n\nTotal numerator rest = {latex(rest_N)}\n\ndelta coefficient = {latex(delta_coeff)}\n\neta coefficient = {latex(eta_coeff)}\n\nDelta coefficient = {latex(Delta_coeff)}\n\nd_beta_Delta coefficient = {latex(d_beta_Delta_coeff)}\n\nTotal rest = {latex(rest)}\n\nControl rest (should be 0) = {control_rest}")
