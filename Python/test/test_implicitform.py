import matplotlib.pyplot as plt
import numpy as np
from newton import newton 
from simpsons import simpsons_total
import itertools

# Constants of equation
# Cs = np.linspace(1, 10, 2)
# Ds = np.linspace(1, 10, 2)

C = 10
D = -2

# CxD = np.array([x for x in itertools.product(Cs, Ds)])

# print(CxD)

def f(t, y, C, D):
    return np.arcsinh(y) + (y/15)*(2/((1 + y**2)**(1/2)) + 1/((1 + y**2)**(3/2)) - 3/((1 + y**2)**(5/2))) - C*t - D

N = 100
ts = np.linspace(0, 10, N + 1)

ys = newton(lambda y : f(ts, y, C, D), ts, 20)
Ys = simpsons_total(ys, ts)

fig, ax = plt.subplots()

ax.plot(ts, ys, label="function")
ax.plot(ts, Ys, label="primitive")
ax.legend()

#for t in CxD:
#    C = t[0]
#    D = t[1]
#    
#    ys = newton(lambda y : f(ts, y, C, D), ts, 20)
#    Ys = np.array([simpsons(ys, h, n, M) for M in range(N + 1)])
#
#    ax.plot(ts, ys, label=f"(C, D) = ({C}, {D})")
#    ax.plot(ts, Ys, label=f"(C, D) = ({C}, {D})")
#    ax.legend()

plt.show()
