from sympy import *
from custom_functions import ei, perp, vector_integral


# --- Symbols ---
# Independent variables
t, alpha, beta = symbols(r"t, \alpha, \beta", real=True)

# Misc symbols
dOmega = symbols(r"\partial\Omega", real=True)


# Defining $\gamma(\alpha, t) = R(\alpha, t) \cdot e^{i\alpha}$
R = Function("R", real=True)
ga = R(alpha, t)*ei(alpha)
ga_beta = ga.subs(alpha, beta)


# Defining $\lambda(t) = r(t) \cdot e^{i\omega(t)}$
r = Function("r", real=True)

# Cartesian form for $\lambda(t)$
#la_x = Function(r"\lambda_x", real=True)
#la_y = Function(r"\lambda_y", real=True)
#la = ImmutableMatrix([la_x(t), la_y(t)])

# Polar form for $\lambda$
omega = Function(r"\omega", real=True)
la = r(t)*ei(omega(t))

diff_ga_ga = ga - ga_beta
diff_ga_la = ga - la
d_alpha_ga = ga.diff(alpha)
d_beta_ga = d_alpha_ga.subs(alpha, beta)


# --- Context of circular patch ---
R_0, r_0, u = symbols(r"R_0, r_0, u", real=True)
omega_0 = u/r_0


# --- Context of singular vorticity ---
eps = symbols(r"\varepsilon", real=True)  # Singular vorticity

# Variation
s = symbols(r"s", real=True)  # Variation variable

Delta = Function(r"\Delta")
delta = Function(r"\delta")
eta = Function(r"\eta")

d_beta_Delta = Delta(beta, t).diff(beta)

theta = symbols(r"\theta")  # Symbol for $beta - omega_0*t$

G = -1/(2*pi*r(t)) * vector_integral(log(diff_ga_la.norm())*ga_beta.diff(beta), (beta, 0, 2*pi)).dot(la)

H = -1/(2*pi*r(t)**2) * vector_integral(log(diff_ga_la.norm())*ga_beta.diff(beta), (beta, 0, 2*pi)).dot(perp(la))

# FF = Matrix([F - R(alpha, t).diff(t), G - r(t).diff(t), H - w(t).diff(t)])
