import numpy as np
from sympy import *

U, F, a, x, y, z, w = symbols(r"\theta, \phi, \alpha, x, y, z, w")

r = Matrix([sin(U)*cos(F), sin(U)*sin(F), cos(U)])
u = r.diff(U)
f = r.diff(F)

f_norm = refine(simplify(sqrt(f.dot(f))), Q.positive(sin(U))) 

f /= f_norm

#print(latex(u, mat_delim='('))

T_BE = Matrix.hstack(r, u, f)
T_EB = simplify(T_BE.inv())

#print(latex(T_EB, mat_delim='('))

R = Matrix([[1, 0, 0], [0, cos(a), -sin(a)], [0, sin(a), cos(a)]])

subs = [(cos(F), x/sqrt(x**2 + y**2)), (sin(F), x/sqrt(x**2 + y**2)), (cos(U), z/sqrt(x**2 + y**2 + z**2)), (sin(U), sqrt(x**2 + y**2)/sqrt(x**2 + y**2 + z**2)), (cos(a), w), (sin(a), sqrt(x**2 + y**2 + z**2))]

T_BE = T_BE.subs(subs) 
T_EB = T_EB.subs(subs)
R = R.subs(subs)

res = simplify(expand(T_BE*R*T_EB))

forth = [(y**2, 1 - z**2 - w**2 - x**2), (y, sqrt(1 - z**2 - w**2 - x**2))]
back = [(sqrt(1 - z**2 - w**2 - x**2), y), (1 - z**2 - w**2 - x**2, y**2)]

# ============================================================
# x^2 + y^2 + z^2 + w^2 = 1
# All rearrangements of the terms between LHS and RHS
# ============================================================


# ------------------------------------------------------------
# x^2 = 1 - y^2 - z^2 - w^2
# ------------------------------------------------------------

forth = [
    (x**2, 1 - y**2 - z**2 - w**2)
]
back = [
    (1 - y**2 - z**2 - w**2, x**2)
]

for i in range(3):
    res = simplify(expand(res.subs(forth)))
    res = simplify(expand(res.subs(back)))


# ------------------------------------------------------------
# y^2 = 1 - x^2 - z^2 - w^2
# ------------------------------------------------------------

forth = [
    (y**2, 1 - x**2 - z**2 - w**2)
]
back = [
    (1 - x**2 - z**2 - w**2, y**2)
]

for i in range(3):
    res = simplify(expand(res.subs(forth)))
    res = simplify(expand(res.subs(back)))


# ------------------------------------------------------------
# z^2 = 1 - x^2 - y^2 - w^2
# ------------------------------------------------------------

forth = [
    (z**2, 1 - x**2 - y**2 - w**2)
]
back = [
    (1 - x**2 - y**2 - w**2, z**2)
]

for i in range(3):
    res = simplify(expand(res.subs(forth)))
    res = simplify(expand(res.subs(back)))


# ------------------------------------------------------------
# w^2 = 1 - x^2 - y^2 - z^2
# ------------------------------------------------------------

forth = [
    (w**2, 1 - x**2 - y**2 - z**2)
]
back = [
    (1 - x**2 - y**2 - z**2, w**2)
]

for i in range(3):
    res = simplify(expand(res.subs(forth)))
    res = simplify(expand(res.subs(back)))


# ============================================================
# Two terms on each side
# ============================================================


# ------------------------------------------------------------
# x^2 + y^2 = 1 - z^2 - w^2
# ------------------------------------------------------------

forth = [
    (x**2 + y**2, 1 - z**2 - w**2)
]
back = [
    (1 - z**2 - w**2, x**2 + y**2)
]

for i in range(3):
    res = simplify(expand(res.subs(forth)))
    res = simplify(expand(res.subs(back)))


# ------------------------------------------------------------
# x^2 + z^2 = 1 - y^2 - w^2
# ------------------------------------------------------------

forth = [
    (x**2 + z**2, 1 - y**2 - w**2)
]
back = [
    (1 - y**2 - w**2, x**2 + z**2)
]

for i in range(3):
    res = simplify(expand(res.subs(forth)))
    res = simplify(expand(res.subs(back)))


# ------------------------------------------------------------
# x^2 + w^2 = 1 - y^2 - z^2
# ------------------------------------------------------------

forth = [
    (x**2 + w**2, 1 - y**2 - z**2)
]
back = [
    (1 - y**2 - z**2, x**2 + w**2)
]

for i in range(3):
    res = simplify(expand(res.subs(forth)))
    res = simplify(expand(res.subs(back)))


# ------------------------------------------------------------
# y^2 + z^2 = 1 - x^2 - w^2
# ------------------------------------------------------------

forth = [
    (y**2 + z**2, 1 - x**2 - w**2)
]
back = [
    (1 - x**2 - w**2, y**2 + z**2)
]

for i in range(3):
    res = simplify(expand(res.subs(forth)))
    res = simplify(expand(res.subs(back)))


# ------------------------------------------------------------
# y^2 + w^2 = 1 - x^2 - z^2
# ------------------------------------------------------------

forth = [
    (y**2 + w**2, 1 - x**2 - z**2)
]
back = [
    (1 - x**2 - z**2, y**2 + w**2)
]

for i in range(3):
    res = simplify(expand(res.subs(forth)))
    res = simplify(expand(res.subs(back)))


# ------------------------------------------------------------
# z^2 + w^2 = 1 - x^2 - y^2
# ------------------------------------------------------------

forth = [
    (z**2 + w**2, 1 - x**2 - y**2)
]
back = [
    (1 - x**2 - y**2, z**2 + w**2)
]

for i in range(3):
    res = simplify(expand(res.subs(forth)))
    res = simplify(expand(res.subs(back)))


# ============================================================
# Three terms on one side
# ============================================================


# ------------------------------------------------------------
# x^2 + y^2 + z^2 = 1 - w^2
# ------------------------------------------------------------

forth = [
    (x**2 + y**2 + z**2, 1 - w**2)
]
back = [
    (1 - w**2, x**2 + y**2 + z**2)
]

for i in range(3):
    res = simplify(expand(res.subs(forth)))
    res = simplify(expand(res.subs(back)))


# ------------------------------------------------------------
# x^2 + y^2 + w^2 = 1 - z^2
# ------------------------------------------------------------

forth = [
    (x**2 + y**2 + w**2, 1 - z**2)
]
back = [
    (1 - z**2, x**2 + y**2 + w**2)
]

for i in range(3):
    res = simplify(expand(res.subs(forth)))
    res = simplify(expand(res.subs(back)))


# ------------------------------------------------------------
# x^2 + z^2 + w^2 = 1 - y^2
# ------------------------------------------------------------

forth = [
    (x**2 + z**2 + w**2, 1 - y**2)
]
back = [
    (1 - y**2, x**2 + z**2 + w**2)
]

for i in range(3):
    res = simplify(expand(res.subs(forth)))
    res = simplify(expand(res.subs(back)))


# ------------------------------------------------------------
# y^2 + z^2 + w^2 = 1 - x^2
# ------------------------------------------------------------

forth = [
    (y**2 + z**2 + w**2, 1 - x**2)
]
back = [
    (1 - x**2, y**2 + z**2 + w**2)
]

for i in range(3):
    res = simplify(expand(res.subs(forth)))
    res = simplify(expand(res.subs(back)))

res = factor(res)

print(latex(res))

#print(rf"{latex(T_BE, mat_delim='(')}\\ \cdot {latex(R, mat_delim='(')}\\ \cdot {latex(T_EB, mat_delim='(')}")
