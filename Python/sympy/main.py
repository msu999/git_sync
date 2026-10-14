from sympy import *
from sympy.abc import i
import numpy as np
import itertools

c_11, c_12, c_21, c_22 = symbols(r"c^1_1, c^2_1, c^1_2, c^2_2")

c_1 = Array([c_11, c_12])
c_2 = Array([c_21, c_22])
c = Array([c_1, c_2]) 
cs = flatten(c)

T_1 = Array([[1, 0], [0, -1]])
T_2 = Array([[0, 1], [1, 0]])
T = Array([T_1, T_2]) 

print(T)

ran = [0, 1]
eqs = zeros(1, len(ran)**3)
for (mu, nu, j) in itertools.product(ran, ran, ran):
    LHS = 0
    RHS = 0
    
    for i in ran:
        LHS += T[i][mu][nu] * c[i][j]

    for (a, b) in itertools.product(ran, ran):
        LHS -= c[mu][a] * T[j][a][b] * c[nu][b]

    I = int(f"{mu}{nu}{j}", 2) 
    eqs[I] = LHS

    print(rf"{latex(LHS)} = {latex(RHS)}\\" + "\n")

sols = solve(eqs, cs, dict=True)

for sol in sols:
    keys = list(sol) 

    c_1_sol = Matrix([sol[keys[0]], sol[keys[1]]])
    c_2_sol = Matrix([sol[keys[2]], sol[keys[3]]])

    print(rf"c_1 = {latex(c_1_sol)} \quad c_2 = {latex(c_2_sol)}\\")
