import sympy as sp
import numpy as np

t, g, m_1, m_2, l, L, k, D, x_0 = sp.symbols(r"t, g, m_1, m_2, l, L, k, D, x_0")

phi = sp.Function(r"\phi")(t)
theta = sp.Function(r"\theta")(t)

alpha, beta, gamma = sp.symbols(r"\alpha, \beta, \gamma")

alpha_expr = (k*l**2 - g*L)/L
beta_expr = -k*(l**2)/L
gamma_expr = k*l*(D - x_0)/L

M = sp.Matrix([[alpha/m_1, beta/m_1], [beta/m_2, alpha/m_2]])

w = sp.symbols(r"\omega")

sol = sp.ImmutableDenseMatrix(sp.solve((M - w*sp.eye(2)).det(), w))
w_1, w_2 = sp.symbols(r"\omega_1, \omega_2", nonzero=True)
W_1, W_2 = sp.Matrix([beta, m_1*w_1 - alpha]), sp.Matrix([beta, m_1*w_2 - alpha])

A_inv = sp.Matrix.hstack(W_1, W_2)
A = A_inv.inv()

u = sp.Matrix([gamma/m_1, gamma/m_2])
v = A*u
Omega = sp.diag(w_1, w_2)
Psi_1, Psi_2 = sp.symbols(r"\Psi_1, \Psi_2", cls=sp.Function) 
Psi = sp.Matrix([Psi_1(t), Psi_2(t)])
ics = {Psi_1(0) : 0, Psi_2(0) : 0}
decoupled_sol = sp.Matrix(sp.dsolve(Psi.diff(t, t) - Omega*Psi - v, [Psi_1(t), Psi_2(t)]))

Psi_1 = decoupled_sol[0].rhs
Psi_2 = decoupled_sol[1].rhs

Psi = sp.Matrix([Psi_1, Psi_2])

Phi = sp.ImmutableMatrix(A_inv*Psi).simplify()
print(sp.latex(Phi))

#print(M - w*sp.eye(2))
#print((M - w*sp.eye(2)).det())
#print(sp.latex(sol.subs({alpha : alpha_expr, beta : beta_expr}).simplify()))


