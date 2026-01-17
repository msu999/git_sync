from sympy import *
import numpy as np

u, v = symbols("u v")

s = Matrix([u - (u**3)/3 + u*v**2, v - (v**3)/3 + v*u**2, u**2 - v**2])
U = Matrix([u, v])

s_u = s.diff(u)
s_v = s.diff(v)
s_uu = s_u.diff(u)
s_uv = s_u.diff(v)
s_vv = s_v.diff(v)

E = s_u.dot(s_u).expand().simplify()
F = s_u.dot(s_v).expand().simplify()
G = s_v.dot(s_v).expand().simplify()

N = ImmutableMatrix(s_u.cross(s_v).expand()).simplify()

print(N)

e = s_uu.dot(N).expand().simplify()
f = s_uv.dot(N).expand().simplify()
g = s_vv.dot(N).expand().simplify()

A = Matrix([[E, F], [F, G]])
B = Matrix([[e, f], [f, g]])

W = A.inv()*B

K = W.det()
H = (W.trace()/2).expand().simplify()
k_1 = (e/E).simplify()
k_2 = (g/G).simplify()

print(f"E = {latex(E)}")
print(f"F = {latex(F)}")
print(f"G = {latex(G)}")
print(f"e = {latex(e)}")
print(f"f = {latex(f)}")
print(f"g = {latex(g)}")
print(f"N = {latex(N)}")
print(f"Gaussian curv. = {latex(K)}")
print(f"Mean curv. = {latex(H)}")
print(f"k_1 = {latex(k_1)}")
print(f"k_2 = {latex(k_2)}")
