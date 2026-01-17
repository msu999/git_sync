from sympy import *
import numpy as np
import matplotlib.pyplot as plt

k, delta, eps, m, z_0 = symbols(r"k, \delta, \varepsilon, m, z_0")
t = symbols("t")

omega_plus = sqrt((k - eps*sqrt(k/delta))/m)
omega_minus = sqrt((k + eps*sqrt(k/delta))/m)

A = z_0/(2*sqrt(k))

T = Matrix([[-sqrt(k), sqrt(k)], [sqrt(delta), sqrt(delta)]])

Psi = A*Matrix([-cos(omega_plus*t), cos(omega_minus*t)])

Phi = T*Psi
# print(latex(Phi))

constant_subs = {k : 0.001, delta : 100, eps : 0.1, m : 0.2, z_0 : -0.1}

lam_Phi = lambdify(t, Phi.subs(constant_subs))

t_0 = 0
t_1 = 3*60
N = 3000

ts = np.linspace(t_0, t_1, N)
Phis = lam_Phi(ts)

zs = Phis[0][0]
thetas = Phis[1][0]

fig_theta, ax_theta = plt.subplots()
ax_theta.set(xlabel=r"$t$", ylabel=r"$\theta$", title=r"$\theta$-plot")
ax_theta.plot(ts, thetas)

fig_z, ax_z = plt.subplots()
ax_z.set(xlabel=r"$t$", ylabel=r"$z$", title=r"$z$-plot")
ax_z.plot(ts, zs)

plt.show()
