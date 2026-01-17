import numpy as np

def f(x, y):
    return x + y

def g(x, *args):
    stencil = f(x, *args)
    res = np.zeros(x.size)

    return f(x, *args)

xs = np.linspace(0, 1, 100)
ys = np.linspace(0, -1, 100)

print(g(xs, ys))
