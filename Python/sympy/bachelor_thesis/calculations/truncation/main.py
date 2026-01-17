import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, '/home/max/Documents/Python/sympy/bachelor_thesis/calculations/variation')
#from X_1 import A_1, A_2

sys.path.insert(0, "/home/max/Documents/Python/pickler")  # Adds pickler folder to system path
from pickler import storeData, loadData  # Imports pickler functions

sys.path.insert(0, "/home/max/Documents/Python/sympy")  # Adds pickler folder to system path
from rk4 import rk4  # Imports RK4

X_1_db = loadData("/home/max/Documents/Python/sympy/bachelor_thesis/calculations/variation/X_1_db")

A_1, A_2 = X_1_db["A_1"], X_1_db["A_2"]

print(f"A_1 = {latex(A_1)}\n\nA_2 = {latex(A_2)}\n")

x, y, k, l = symbols(r"x, y, k, l", real = True)

B_int = A_1*(1 + exp(I*k*y)) + I*k*A_2*(exp(I*k*y) - 1)
B_int_trig = A_1*(1 + cos(k*y) + I*sin(k*y)) + I*k*A_2*(cos(k*y) + I*sin(k*y) - 1)

I_k_B_int = integrate(B_int, (k, 0, k))

for i, t in enumerate(I_k_B_int.args):
    print(f"term {i} = {t}\n")

for i, t in enumerate(I_k_B_int.args[2].args):
    print(f"subterm {i} = {t}\n")

print(f"General B_int = {latex(B_int_trig.expand().simplify())}\n\nB_int integrated wrt k = {latex((I_k_B_int.args[0] + I_k_B_int.args[1] + I_k_B_int.args[2].args[1][0]).expand().simplify())}\n\n")

# 

lam_B_int = lambdify([k, y], B_int)

cur_k = -1
k_subs = {k : cur_k}
print(f"\n--- Current k = {cur_k} ---\n")
print(f"B_int = {latex(B_int.subs(k_subs))}")

eps = 1e-6
N = 100
rk4_int = rk4(lambda x, y : lam_B_int(cur_k, x), np.array(0), eps, 2*np.pi - eps, N, iscomplex=True)
print(f"rk4 result = {rk4_int[1][0][-1]}")

#B = integrate(B_int.subs(k_subs), (y, 0, 2*pi))

#print(f"explicit result = {latex(B)}")

N_phi = 1
m = 3
n = 3

plt.show()
