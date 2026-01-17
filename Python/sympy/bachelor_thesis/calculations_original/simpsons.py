import numpy as np
import matplotlib.pyplot as plt

def simpsons(f, a, b, N, *args):
    xs = np.linspace(a, b, N)
    Fs = np.zeros(xs.size)

    for i in range(1, xs.size):
        v = xs[i - 1]
        u = xs[i]

        Fs[i] = Fs[i - 1] + ((u - v)/6) * (f(v, *args) + 4*f((v + u)/2, *args) + f(u, *args))

    return (xs, Fs)


## Example integration of f(x, y) = cos(x + y), that is:
##
##    /x                             t = x
##    |  cos(t + y) dt = [sin(t + y)]       = sin(x + y) - sin(y)
##    /0                             t = 0
##
#
#def f(x, y):
#    return np.cos(x + y)
#
#fig, ax = plt.subplots()
#
#a = 0
#b = 2*np.pi
#N = 100
#
#y = np.pi
#
#res = simpsons(f, a, b, N, y)
#
#xs = res[0]
#Fs = res[1]
#
#ax.plot(xs, Fs)
#
#plt.show()
