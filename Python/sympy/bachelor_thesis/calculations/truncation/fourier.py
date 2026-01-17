import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.insert(0, "/home/max/Documents/Python/sympy")

from rk4 import rk4  # Imports RK4

#################
#  Integration  #
#################

def integral(f, a, b, step_size=0.01, is_complex=False):
    h = step_size
    xs = np.arange(a, b + h, h)
    vectorized_f = np.vectorize(f)
    fs = vectorized_f(xs)
    #print(f"xs = {xs}, fs = {fs}")

    if is_complex:
        Fs = np.zeros_like(xs, dtype=complex)

    else:
        Fs = np.zeros_like(xs)

    #Fs[1] = h/2 * (fs[0] + fs[1])

    for i in range(2, xs.size):
        dF = h/3 * (fs[i - 2] + 4*fs[i - 1] + fs[i])
        Fs[i] = Fs[i - 2] + dF

    return Fs[-1]

def array_integral(fs, xs, is_complex=False):
    h = xs[1] - xs[0]
    #print(f"xs = {xs}, fs = {fs}")

    if is_complex:
        Fs = np.zeros_like(xs, dtype=complex)

    else:
        Fs = np.zeros_like(xs)

    #Fs[1] = h/2 * (fs[0] + fs[1])

    for i in range(2, xs.size):
        dF = h/3 * (fs[i - 2] + 4*fs[i - 1] + fs[i])
        Fs[i] = Fs[i - 2] + dF

    return Fs[-1]

def array_integral_2d(fs, Xs, Ys, is_complex=False):
    h = Xs[1] - Xs[0]
    k = Ys[1] - Ys[0]
    xs = Xs[:,0]
    ys = Ys[0,:]

    print(f"xs = {xs}\n\nys = {ys}")
    
    Fs_x = np.array([array_integral(fs[i,:], xs, is_complex) for i in range(xs.size)])

    print(Fs_x)

    return array_integral(Fs_x, ys, is_complex)


##############################
#  Fourier series functions  #
##############################

default_step_size = 0.01
default_xs = np.arange(0, 2*np.pi + default_step_size, default_step_size)
def coeff(f, n, step_size=0.01, bound_shave=0):
    if step_size != default_step_size:
        xs = np.arange(0, 2*np.pi + step_size, step_size)
        fs = f(xs)*np.exp(-1j*n*xs)

        integ = array_integral(fs, xs, True)

    else:
        fs = f(default_xs)*np.exp(-1j*n*xs)

        integ = array_integral(fs, default_xs, True)

    #integ = integral(lambda x : f(x)*np.exp(-1j*n*x), bound_shave, 2*np.pi - bound_shave, step_size, True)

    return (1/(2*np.pi)) * integ

def coeff_2d(f, k, l, step_size=0.01, bound_shave=0):
    integ = integral(lambda x : integral(lambda y : f(x, y)*np.exp(-1j*(k*x + l*y)), bound_shave, 2*np.pi - bound_shave, step_size=step_size, is_complex=True), bound_shave, 2*np.pi - bound_shave, step_size=step_size, is_complex=True)

    return 1/(4*np.pi**2) * integ

def array_coeff_2d(f, k, l, step_size=0.01, bound_shave=0):
    if step_size != default_step_size:
        Xs, Ys = np.meshgrid(default_xs, default_xs, indexing="ij")
        fs = f(Xs, Ys)*np.exp(-1j*(k*Xs + l*Ys))

        integ = array_integral_2d(fs, Xs, Ys, True)

    else:
        xs = np.arange(0, 2*np.pi + step_size, step_size)
        Xs, Ys = np.meshgrid(default_xs, default_xs, indexing="ij")
        fs = f(Xs, Ys)*np.exp(-1j*(k*Xs + l*Ys))

        integ = array_integral_2d(fs, Xs, Ys, True)

    #integ = integral(lambda x : integral(lambda y : f(x, y)*np.exp(-1j*(k*x + l*y)), bound_shave, 2*np.pi - bound_shave, True), bound_shave, 2*np.pi - bound_shave, True)

    return 1/(4*np.pi**2) * integ

def series(coeffs, index_shifts, *xs):
    # Corrects type of index_shifts for program below if a single integer has been passed
    if type(index_shifts) == int:
        index_shifts = [index_shifts]

    result = 0

    for (I, c) in np.ndenumerate(coeffs):
        exponent = 0

        for (var_num, x) in enumerate(xs):
            #print(f"I = {I}, var_num = {var_num}, index_shifts = {index_shifts}, x = {x}")
            exponent += (I[var_num] + index_shifts[var_num])*x
        
        result += c*np.exp(1j*exponent)

    return result



####################
#  Testing ground  #
####################

if __name__ == "__main__":
    def f(x, y):
        return np.exp(1j*(x + y)) 

    xs = np.linspace(0, 1, 100)
    Xs, Ys = np.meshgrid(xs, xs, indexing="ij")
    print(f"Xs = {Xs}\n\nYs = {Ys}")
    fs = f(Xs, Ys)
    
    print(f"Coeff = {array_coeff_2d(f, 1, 1, step_size=0.0001)}")
    #print(f"Coeff = {coeff_2d(f, 1, 1, step_size=0.001)}")
