import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from sympy.solvers.solveset import linsolve

import sys

sys.path.insert(0, "/home/max/Documents/Python/") 

from numnal import one


# Single diff quotient functions
def central_diff(_f, f_, dx):
    return (f_ - _f)/(2*dx)

def three_point_forward_diff(f, f_, f__, dx):
    return (-f__ + 4*f_ - 3*f)/(2*dx)

def three_point_backward_diff(__f, _f, f, dx):
    return (__f - 4*_f + 3*f)/(2*dx)


def diff(fs, xs):
    ########################################################
    #             df                                       #
    #  Returns:   --(x)                                    #
    #             dx                                       #
    #                                                      #
    #  at x in xs = (a, a + dx, ... , a + (N - 1)*dx = b)  #
    ########################################################
    

    dfs = np.zeros_like(fs)
    dx = xs[1] - xs[0]

    # Boundary derivatives at a and b evaluated using forward/backward three point diff
    dfs[0] = three_point_forward_diff(fs[0], fs[1], fs[2], dx)
    dfs[-1] = three_point_backward_diff(fs[-3], fs[-2], fs[-1], dx)

    # Derivatives inside (a, b) evaluated using central diff
    N = xs.size
    for n in range(1, N-1):
        dfs[n] = central_diff(fs[n-1], fs[n+1], dx)

    return dfs


def grad(fs, xs):
    grid_shape = xs.shape
    num_vars = len(grid_shape)
    print(f"num_vars = {num_vars}\n\ngrid_shape = {grid_shape}\n\nLast point in xs = {xs[*(grid_shape - np.array([1 for i in grid_shape]))]}", end="\n\n")
    dfs = np.empty(grid_shape, dtype=tuple)
    print(f"One vector test = {one(num_vars, 0)}", end="\n\n")
    dx = (xs[*one(num_vars, 0)] - xs[*np.zeros(num_vars, dtype=int)])[0]
    print(f"dx = {dx}", end="\n\n")

    for (I, x) in np.ndenumerate(xs):
        nabla_f = np.zeros(num_vars)

        for i in range(num_vars):
            dI = one(num_vars, i)

            if I[i] - 1 == -1:
                #print("Forward diff:")
                #print(f"fs[I] = {fs[*I]}, fs[I + dI] = {fs[*(I + dI)]}, fs[I + 2dI] = {fs[*(I + 2*dI)]}", end="\n\n")
                res = three_point_forward_diff(fs[*I], fs[*(I + dI)], fs[*(I + 2*dI)], dx)
                #print(f"res = {res}", end="\n\n")
                nabla_f[i] = res

            elif I[i] + 1 == grid_shape[i]:
                #print("Backward diff:")
                #print(f"fs[I - 2dI] = {fs[*(I - 2*dI)]}, fs[I - dI] = {fs[*(I - dI)]}, fs[I] = {fs[*I]}", end="\n\n")
                res = three_point_backward_diff(fs[*(I - 2*dI)], fs[*(I - dI)], fs[*I], dx)
                #print(f"res = {res}", end="\n\n")
                nabla_f[i] = res 

            else:
                #print("Central diff:")
                #print(f"fs[I - dI] = {fs[*(I - dI)]}, fs[I + dI] = {fs[*(I + dI)]}", end="\n\n")
                res = central_diff(fs[*(I - dI)], fs[*(I + dI)], dx)

                #print(f"res = {res}", end="\n\n")
                nabla_f[i] = res 

        dfs[*I] = nabla_f

    return dfs


def grid(Xs):
    shape = [xs.size for xs in Xs]
    grid = np.empty(shape, dtype=np.ndarray)

    for (I, _) in np.ndenumerate(grid):
        grid[I] = np.array([Xs[var_num][i] for (var_num, i) in enumerate(I)])
    
    return grid


def f(x, y):
    return x * np.exp(-x**2 - y**2)

#fs = np.array([[1, 2, 3], [2, 2, 3], [3, 3, 3]])
xs = np.arange(-2, 2.2, 0.2)
grid = grid((xs, xs))
fs = np.empty_like(grid)

for (I, x) in np.ndenumerate(grid):
    fs[I] = f(*x)

print(grid, end="\n\n")
coord = np.array([2, 1])
print(f"{coord + (1, 1)} = {grid[*coord]}")

dfs = grad(fs, grid)

print(dfs)

#xs = np.linspace(-5, 5, 50)
#fs = f(xs)
#dfs = diff(fs, xs)
#ddfs = diff(dfs, xs)
#dddfs = diff(ddfs, xs)
#
fig = plt.figure()
ax = plt.axes()
#
#ax.plot(xs, fs, label="f")
#ax.plot(xs, dfs, label="df")
#ax.plot(xs, ddfs, label="ddf")
#ax.plot(xs, dddfs, label="dddf")

for (I, x) in np.ndenumerate(grid):
    #print(f"Plotting arrow = {dfs[I]} at x = {x}")
    ax.quiver(*x, *dfs[I], scale=10)

#fig.legend()

plt.show()
