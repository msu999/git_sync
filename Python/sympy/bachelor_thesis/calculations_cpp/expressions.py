from sympy import *
from custom_functions import ei, perp, vector_integral


# --- Symbols ---
# Independent variables
t, x, y = symbols(r"t, x, y", real=True)

# Defining $\gamma_0(x) = R(x) \cdot e^{ix}$
f = Function("f", real=True)
a = symbols(r"a", positive=True)
ga = sqrt(1 + f(x))*Matrix([a*cos(x), sin(x)])
ga_y = ga.subs(x, y)


# Defining $\lambda_0 = (0, h)$
h = symbols("h", real=True)
la = Matrix([0, h])


# Defining $\Omega$
Omega = symbols(r"\Omega", real = True)


# Common expressions involving $\gamma_0$ and $\lambda_0$
diff_ga_ga = ga - ga_y
diff_ga_la = ga_y - la
d_x_ga = ga.diff(x)
d_y_ga = d_x_ga.subs(x, y)


# --- Context of circular patch ---
f_0 = 0
h_0, Omega_0 = symbols(r"h_0, \Omega_0", real=True)


# --- Context of singular vorticity ---
eps = symbols(r"\varepsilon", real=True)  # Singular vorticity

# Variation
s = symbols(r"s", real=True)  # Variation variable

phi = Function(r"\phi")
delta = symbols(r"\delta")
eta = symbols(r"\eta")

d_y_phi = phi(y).diff(y)

theta = symbols(r"\theta")  # Symbol for $y - omega_0*t$
