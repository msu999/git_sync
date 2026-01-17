from sympy import *
from custom_functions import ei, perp, vector_integral


# --- Symbols ---
# Independent variables
t, x, y = symbols(r"\theta, x, y", real=True)

# Defining $\gamma(x, t) = \sqrt(1 + f(x, t)) \cdot e^{ix}$
f = Function("f", real=True)
ga = sqrt(1 + f(x, t))*ei(x)
ga_y = ga.subs(x, y)


# Defining $\lambda(t)$
la_1 = Function(r"\lambda_1", real=True)
la_2 = Function(r"\lambda_2", real=True)
la = Matrix([la_1(t), la_2(t)])


# Defining $\omega$
omega = symbols(r"\omega", real = True)


# Common expressions involving $\gamma$ and $\lambda$
diff_ga_ga = ga - ga_y
diff_ga_la = ga_y - la
d_x_ga = ga.diff(x)
d_y_ga = d_x_ga.subs(x, y)


# --- Context of circular patch ---
f_0 = 0
h_0 = symbols(r"h_0", positive=True)
la_1_0 = 0
la_2_0 = h_0
la_0 = ImmutableMatrix([la_1_0, la_2_0])
omega_0 = symbols(r"\omega_0", real=True)



# --- Context of singular vorticity ---
eps = symbols(r"\varepsilon", real=True)  # Singular vorticity

# Variation
s = symbols(r"s", real=True)  # Variation variable

phi = Function(r"\phi")
delta_1 = Function(r"\delta_1")
delta_2 = Function(r"\delta_2")
delta = ImmutableMatrix([delta_1(t), delta_2(t)])
eta = symbols(r"\eta")

d_y_phi = phi(y).diff(y)
