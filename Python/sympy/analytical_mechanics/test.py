from rk4 import rk4
import matplotlib.pyplot as plt
import numpy as np

# -------------------- INSTRUCTIONS -----------------------
#
# We have the following ODE problem:
#
#   y' = f(t, y)
#   y(t_0) = y_0
#
# For example:
#
#   x'' = -x
#   x(0) = 0
#   x'(0) = 1
#
#   (Solution x = sin(t))
#
# Gets rewritten into:
#
#   y_1 = x
#   y_2 = x'
#
# Which yields the ODE problem:
#
#   y' = (y_1', y_2') = (x', x'') = (y_2, -y_1) =: f(t, y)
#
# ---------------------------------------------------------

M = 5.97219e24
G = 6.6743e-11
C = M*G
k = 10
h_0 = 1000*1000
R = 6378*1000
r_0 = R + h_0
w_0 = np.sqrt(C/(r_0**3))

def A(y_0, y_1, y_3):
    return y_0*(y_3**2) - C/(y_0**2) - k*(y_1**2)

def B(y_0, y_1, y_3):
    return (-k*y_0*(y_3**2) - 2*y_1*y_3)/y_0

def f(t, y):
    return np.array([y[1], A(y[0], y[1], y[3]), y[3], B(y[0], y[1], y[3])])

t_0 = 0
t_1 = 1*60
y_0 = np.array([r_0, 0, 0, w_0])

# Number of steps in RK4 iteration
N = 1000

# RK4 function returns sol = (t, y) where y = (y_1, y_2, ...). Hence to plot
# y_i we plot sol[0] against sol[1][i].
sol = rk4(f, y_0, t_0, t_1, N)

r = sol[1][0]
theta = sol[1][2]

x = r*np.cos(theta)
y = r*np.sin(theta)

fig, ax = plt.subplots()
ax.set_aspect("equal")

ax.plot(x, y)

plt.show()
