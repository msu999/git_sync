import matplotlib.pyplot as plt
import numpy as np
from euler_methods import euler_forward

C = 0.01

def f(t, y):
    a = y[0]

    return np.array([(np.sqrt(1 + a**2)*C)/(1 + (a**2)/((1 + a**2)**3)), a])

t_0 = -100
t_1 = 10

s_0 = 0
a_0 = 1
y_0 = np.array([a_0, s_0])

h = 0.1

res = euler_forward(f, t_0, t_1, y_0, h)

fig, ax = plt.subplots()

ax.plot(res[0], res[1][0], label="derivative")
ax.plot(res[0], res[1][1], label="function")
ax.legend()

plt.show()
