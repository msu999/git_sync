from rk4 import rk4
import matplotlib.pyplot as plt
import numpy as np

A = 0.1*0.1
p_0 = 1e5
m_s = 0.5

V_st = 0.1*A

rho = 1323
m = rho*V_st
M = 32*1e-3
n = m/M
N_0 = 6.022*1e23
N = n*N_0

k = 1.380*1e-23
T = 300
mu = 100

a = N*k*T/m_s
b = A*p_0/m_s

x_stable = a/b
eps = -100
x_0 = x_stable + eps
v_0 = 0

def friction(v):
    return -mu*np.arctan(10*v)

def f(t, y):
    return np.array([y[1], (a/y[0]) - b + friction(y[1])])

#def f(t, y):
#    return np.array([y[1], -y[0]])

t_0 = 0
t_1 = 100
y_0 = np.array([x_0, v_0])

Iter = int(1e5)

sol = rk4(f, y_0, t_0, t_1, Iter)

fig, ax = plt.subplots()

ax.plot(sol[0], sol[1][0], color="black", label="displacement")
ax.plot(sol[0], sol[0]*0 + x_stable, "--r", label="x_stable = " + str(x_stable))

ax.legend()

plt.show()
