import numpy as np

def diff(f, x, h):
    return (f(x + h) - f(x - h))/(2*h)

def newton(f, x_0, N, diff_h = 0.01):
    x = x_0

    for i in range(N):
        x = x - f(x)/diff(f, x, diff_h)

    return x
