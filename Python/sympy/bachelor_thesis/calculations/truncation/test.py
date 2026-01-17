import time
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider



################
#  PARAMETERS  #
################

## Theoretical parameters ##
eps = 0  # Epsilon in problem

## Numerical parameters ##
N = 100 + 1
h = 2*np.pi/(N - 1)  # Differentiation step size

x_coeff_res = 2  # Resolution of truncation in x (2*x_coeff_res + 1 = nr of x coeffs)
t_coeff_res = 2  # Resolution of truncation in t (2*t_coeff_res + 1 = nr of t coeffs)
num_x_coeffs = 2*x_coeff_res + 1
num_t_coeffs = 2*t_coeff_res + 1

x_shift = -x_coeff_res  # index_(x,py) + x_shift = index_(x,real)
t_shift = -t_coeff_res  # index_(t,py) + t_shift = index_(t,real)



##################################
#  USEFUL QUANTITIES AND ARRAYS  #
##################################

## x and t arrays ##
buffer = 0
x = np.linspace(0, 2*np.pi, N, dtype=np.float128)
t = np.linspace(0, 2*np.pi, N, dtype=np.float128)
X, T = np.meshgrid(x, t, indexing="ij")
Y, _, __ = np.meshgrid(x, x, t, indexing="ij")  # Put Y in front for simple numpy broadcasting
one = 0*Y + 1  # For stretching arrays from X x T domain into Y x X x T domain

ei = np.array([np.cos(X), np.sin(X)])  # e^(i*X) = (cos(X), sin(X))

h_0 = np.pi  # Initial height of rotating point
omega_0 = 1/(2*h_0**2)  # Initial angular velocity

M = -1/(2*h_0**2)
nu = h_0/2  # For use in matrix used in inversion for delta_+-1 and eta



######################
#  USEFUL FUNCTIONS  #
######################

## Index translation functions (from numpy "i, j"-indices to fourier "k, l"-indices using x- and t-shift) ##
def x_index_real_to_py(k):
    # Takes index_(x,real) and returns index_(x,py)
    return int(k-x_shift)

def t_index_real_to_py(k):
    # Takes index_(t,real) and returns index_(t,py)
    return int(k-t_shift)

def x_index_py_to_real(i):
    # Takes index_(x,real) and returns index_(x,py)
    return int(i+x_shift)

def t_index_py_to_real(i):
    # Takes index_(t,real) and returns index_(t,py)
    return int(i+t_shift)


## Vector analysis functions ##
def dot(u, v):
    return u[0]*v[0] + u[1]*v[1]

def norm(u):
    return np.sqrt(dot(u, u))

def perp(u):
    return np.array([-u[1], u[0]])


def riemann(f, axis=0):
    f = np.moveaxis(f, axis, 0)  # Move axis of integration to first position

    res = h * (np.sum(f, axis=0) - f[0])

    print(f"f shape = {f.shape}\nRes shape = {res.shape}\n")

    return res 


## Fourier analysis functions ##
k = np.linspace(-x_coeff_res, x_coeff_res, 2*x_coeff_res + 1, dtype=int)
l = np.linspace(-t_coeff_res, t_coeff_res, 2*t_coeff_res + 1, dtype=int)
KK, LL, XX, TT = np.meshgrid(k, l, x, t, indexing="ij")
series_exponentials_2d = np.exp(1j*(KK*XX + LL*TT))  # Array containing exponentials e^(i(kX + lT))
i_0 = x_index_real_to_py(0)
series_exponentials_1d = series_exponentials_2d[i_0,:,:,:]
coeff_exponentials_2d = np.flip(series_exponentials_2d, (0, 1))
#coeff_exponentials_2d = np.exp(-1j*(KK*XX + LL*TT))
#coeff_exponentials_1d = np.flip(series_exponentials_1d, 0)
one_klxt = 0*TT + 1  # Array in KK x LL x XX x TT space for expanding X x T arrays into said space

def coeffs(f, no_x=False):
    f = f*one_klxt  # Expands f array (in X x T space) into KK x LL x XX x TT space 

    xt_integrand = f * coeff_exponentials_2d

    x_integrand = riemann(xt_integrand, axis=-1)  # Perform integration along TT axis to begin with. Returns the integrand that goes into the XX integral (of shape KK x LL x XX)

    # If no x dependence then we are done and can ignore the KK- and XX-axis (integration along this will just yield a factor 2pi)
    if no_x:
        result = 1/(2*np.pi) * x_integrand[i_0,:,0]

        return result

    # Else we perform integration also along XX-axis
    else:
        result = 1/(4*np.pi**2) * riemann(x_integrand, axis=-1)  # Perform integration along XX axis. Returns coefficients of f (shape KK x LL).

        return result


def coeff(f, k, l, no_x=False):
    #print(f"f = {f}\n\nf shape = {f.shape}\n\nexp = np.exp(-1j*(k*X + l*T))\n\nexp shape = {(np.exp(-1j*(k*X + l*T))).shape}\n\nf * exp shape = {(f * np.exp(-1j*(k*X + l*T))).shape}")
    if no_x:
        k = 0

    xt_integrand = f * np.exp(-1j*(k*X + l*T))

    x_integrand = riemann(xt_integrand, axis=-1)
    #print(f"x integrand shape = {x_integrand.shape}")

    if no_x:
        print(f"x_integrand = {x_integrand}")
        intg = 1/(2*np.pi) * x_integrand[0]

        return intg

    else:
        intg = 1/(4*np.pi**2) * riemann(x_integrand)

        #print(f"intg shape = {intg.shape}")

        return intg


def series(coeffs, no_x=False):
    result = np.zeros_like(X, dtype=np.float128)

    if no_x:
        result = np.sum((series_exponentials_1d.T*coeffs).T, axis=0)

        return result

    else:
        result = np.sum(np.sum((series_exponentials_2d.T*coeffs.T).T, axis=0), axis=0)

        return result


K, L = np.meshgrid(k, l)


f_true_coeffs = np.random.rand(num_t_coeffs)
#print(f"Array_to_sum shape = {array_to_sum.shape}")
f = series(f_true_coeffs, no_x=True)

f_numerical_coeffs = np.zeros((num_x_coeffs, num_t_coeffs))
for (I, _) in np.ndenumerate(f_numerical_coeffs):
    k = x_index_py_to_real(I[0])
    l = t_index_py_to_real(I[1])
    f_numerical_coeffs[*I] = coeff(f, k, l, no_x=True)

#f_numerical_coeffs = coeffs(f, no_x=True)

print(f"True coeffs = {f_true_coeffs}\n\nNumerical coeffs = {f_numerical_coeffs}")

g_coeffs = np.zeros((num_x_coeffs, num_t_coeffs))
i_plus = x_index_real_to_py(1)
i_minus = x_index_real_to_py(-1)
j = t_index_real_to_py(0)
g_coeffs[i_plus, j] = 1/2
g_coeffs[i_minus, j] = 1/2
g = series(g_coeffs)

fig = plt.figure()
ax = plt.axes(projection="3d")

ax.plot(X, T, g, label="Numerical")
ax.plot(X, T, np.cos(X), label="True")

ax.legend()

plt.show()
