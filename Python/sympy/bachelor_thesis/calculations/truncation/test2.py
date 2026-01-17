import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.insert(0, "/home/max/Documents/Python/sympy/")

from rk4 import rk4, simpsons
from numnal import evalf, cross, slicer


def integrator(fs, xs, i):
    Fs = np.empty_like(fs)

#    for (I, x) in np.ndenumerate(xs):

    return None

def f(x):
    return x**2

a = 5
xs = np.linspace(-a, a, 101)
Xs = cross(xs, xs)
print(slicer(Xs, 1, (0, 0)))

fs = evalf(f, xs)

#print(fs)

fs_x_fs = cross(fs, fs)
#print(fs_x_fs)

Fs = np.zeros_like(fs_x_fs)

for (I, f_x_f) in np.ndenumerate(fs_x_fs):
    
    F = 1
    for f in f_x_f:
        F *= f

    Fs[*I] = F

print(f"Fs = {Fs}\n")
gs = simpsons(Fs, Xs, 1)
print(gs)
last_gs = slicer(gs, 0, (0, -1))

#print([:, 0])
#print(Fs[8, 9])

#Fs = simpsons(fs, xs)
#print(Fs)

fig, ax = plt.subplots()

ax.plot(xs, last_gs, color="b", label="numerical")
ax.plot(xs, 2/3 * a**3 * xs**2, color="r", label="true")

ax.legend()

plt.show()
