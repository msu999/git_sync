import numpy as np
import matplotlib.pyplot as plt
from simpsons import simpsons

def f(beta, x, k):
    return np.sin(beta - x)*np.cos(beta - x)*np.sin(k*beta)/(1 - np.cos(beta - x))

fig, ax = plt.subplots()

a = 0
b = 2*np.pi
N = 200

k = 1

resolution = 1000

xs = np.linspace(a, b, N)
print(simpsons(f, a, b, 1, k))
res = np.array([simpsons(f, a, b, resolution, x, k)[1][-1] for x in xs])

#print(res)

Fs = res

ax.plot(xs, Fs)

plt.show()
