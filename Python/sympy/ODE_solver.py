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


def f(t, y):
    return np.array([])

t_0 = 0
t_1 = 1
y_0 = np.array([])

# Number of steps in RK4 iteration
N = 100

# RK4 function returns sol = (t, y) where y = (y_1, y_2, ...). Hence to plot
# y_i we plot sol[0] against sol[1][i].
sol = rk4(f, y_0, t_0, t_1, N)

fig, ax = plt.subplots()

ax.plot(sol[0], sol[1][])

plt.show()
