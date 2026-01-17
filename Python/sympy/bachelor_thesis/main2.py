from sympy import *

# --- Custom functions ---
def ei(alpha):
    return ImmutableMatrix([cos(alpha), sin(alpha)])

def perp(v):
    if isinstance(v, ImmutableMatrix) and v.shape == (2, 1):
        return ImmutableMatrix([-v[1], v[0]])

    return None

def vector_integral(v, t):
    if isinstance(v, ImmutableMatrix) and v.shape == (2, 1):
        return ImmutableMatrix([Integral(v[0], t), Integral(v[1], t)])

    return None

# --- Symbols ---
x, t, dOmega, R_0, alpha, beta, eps, w = symbols("x, t, \\partial\\Omega, R_0, \\alpha, \\beta, \\varepsilon, w", real=True)

# Plugging in $\gamma(\alpha, t) = R(\alpha, t)e^{i\alpha}$ into:
# $\partial_t R(\alpha, t)=\frac{1}{2 \pi R(\alpha, t)}\left[\int_{\partial \Omega} \log |\gamma(\alpha, t)-y| d y-\varepsilon \cdot \frac{(\gamma(\alpha, t)-\lambda(t))^{\perp}}{|\gamma(\alpha, t)-\lambda(t)|^2}\right] \bullet \partial_\alpha \gamma(\alpha, t)^{\perp}$
R = Function("R", real=True)
ga = R(alpha, t)*ei(alpha)
ga_beta = ga.subs(alpha, beta)
r = Function("r", real=True)

# Cartesian form for $\lambda(t)$
#la_x = Function(r"\lambda_x", real=True)
#la_y = Function(r"\lambda_y", real=True)
#la = ImmutableMatrix([la_x(t), la_y(t)])

# Polar form for $\lambda$ 
w = Function(r"\omega", real=True)  # Needed for alternative form for $\lambda$
la = r(t)*ei(w(t))  # Alternative form for $\lambda$

# Log integrand:
diff_ga_ga = ga - ga_beta
log_integrand = log((diff_ga_ga).norm())*ga_beta.diff(beta).expand().simplify()

print(rf"Log integrand in $\partial_t R$ calculation = ${latex(log_integrand)}$" + "\n")

# Epsilon term:
diff_ga_la = ga - la
epsilon_term = (perp(diff_ga_la)/(diff_ga_la.norm()**2)).dot(ga.diff(alpha)).expand().simplify()

print(rf"Epsilon term = ${latex(epsilon_term)}$" + "\n")


# Plugging in $\gamma(\alpha, t) = R(\alpha, t)e^{i\alpha}$ into:
# $-\frac{1}{2 \pi} \int_\gamma \log |\lambda-y| dy$
log_integrand = log((la - ga_beta).norm()).expand().simplify() 

print(rf"Log integrand in $d\lambda/dt$ calculation = ${latex(log_integrand)}$"+ "\n")

a, b = symbols("a, b")

F = 1/(2*pi*R(alpha, t)) * (vector_integral(log((diff_ga_ga).norm())*ga_beta.diff(beta), (beta, 0, 2*pi)) - eps*perp(diff_ga_la)/(diff_ga_la.norm()**2)).dot(ga.diff(alpha))

G = -1/(2*pi*r(t)) * vector_integral(log(diff_ga_la.norm())*ga_beta.diff(beta), (beta, 0, 2*pi)).dot(la)

H = -1/(2*pi*r(t)**2) * vector_integral(log(diff_ga_la.norm())*ga_beta.diff(beta), (beta, 0, 2*pi)).dot(perp(la))

FF = Matrix([F - R(alpha, t).diff(t), G - r(t).diff(t), H - w(t).diff(t)])

R_0, r_0, u, s = symbols("R_0, r_0, u, s")
w_0 = u/r_0
Delta = Function(r"\Delta")
delta = Function(r"\delta")
eta = Function(r"\eta")

g = FF.subs({R(alpha, t) : R_0 + s*Delta(alpha, t), R(beta, t) : R_0 + s*Delta(beta, t), r(t) : r_0 + s*delta(t), w(t) : w_0 + s*eta(t)})

print(latex(g[0].diff(s).simplify()) + "\n")



f_x = Function("f_x")
f_y = Function("f_y")
f = Matrix([f_x(x, t), f_y(x, t)])
g_x = Function("g_x")
g_y = Function("g_y")
g = Matrix([g_x(t), g_y(t)])

test = integrate(f, (x, a, b))
print(latex(test.dot(g)))
