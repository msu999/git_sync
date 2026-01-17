import numpy as np
import matplotlib.pyplot as plt
from fourier import coeff, series

import sys

sys.path.insert(0, "/home/max/Documents/Python/") 

from numnal import cross, slicer, evalf
from rk4 import simpsons

def ei(alpha):
    return np.array([np.cos(alpha), np.sin(alpha)])

# Fundamental variables
eps = 0.01

f_0 = 0
h_0 = np.pi
omega_0 = 1/(2*h_0**2)

N_t = 100
t_0 = 0
t_1 = 1
t = np.linspace(t_0, t_1, N_t)

N_x = 100
x_0 = 0
x_1 = 2*np.pi
x = np.linspace(x_0, x_1, N_x)
xt = cross(x, t)
#xyt = cross(xy, t)

la_0 = h_0*ei(omega_0*t)

P_0 = (f_0, la_0, omega_0)


def test_f(x, t):
    return 1


def test_la(x):
    return np.array([np.cos(x), np.sin(x)])


def perp(A):
    result = np.empty_like(A)

    for (I, a) in np.ndenumerate(A):
        result[*I] = np.array([-a[1], a[0]])

    return result

def dot(A, B):
    result = np.empty_like(A)
    shape_A = A.shape

    for I in np.ndindex(shape_A):
        result[*I] = A[*I].dot(B[*I])

    return result


f = evalf(test_f, xt) 
ga = np.sqrt(1 + f) * evalf(lambda x, t : ei(x), xt)
la = evalf(lambda x, t : test_la(x), xt)
diff_ga_la = ga - la
eps_term = -eps/np.pi * dot(perp(ga - la), )
#ga_x_ga = cross(ga, ga)

print(f"x = {x}\n\nf = {f}\n\ngamma = {ga}\n\nla = {la}\n\nga - la = {diff_ga_la}\n\n(ga - la)^P = {perp(ga - la)}")

def Y(P):
    # f : (x, theta) |-> R
    #
    #

    (f, la, omega) = P

    def Y_1(P):
        def X_1(P):
            (f, la, omega) = P
            f_y = lambda y : f(y, )

            ## Creating integrand expression ##
            fxf = cross(f, f)  # Cross f values with themselves to get a function from (x, y) to (f(x), f(y)) ( shape = (f.shape, f.shape), elements = (*, *) ).
            integrand_domain = xyts
            integrand = np.empty_like(fxf)  # Prepare integrand ( shape = fxf.shape, elements = floats ).
            
            for (I, (f_x, f_y)) in np.ndenumerate(fxf):
                integrand[*I] = np.log(-2*np.sqrt(f_x + 1)*np.sqrt(f_y + 1)*np.cos(x - y) + f_x + f_y + 2)/np.pi
            
            integral_part = (1/np.pi) * simpsons(integrand, integrand_domain, 1)
            eps_part

            return 

        return omega*np.gradient(f, axis=1) - X_1(P)

    def Y_4(P):
        delta_1 = P[1][0]

        return delta_1[0]





def p_0(x):
    return (eps/np.pi)*((h_0*np.cos(x))/(h_0**2 - 2*h_0*np.sin(x) + 1))

def B(k):
    if k == 0:
        return 0

    elif abs(k) == 1:
        return 1j*np.sign(k)*(3/2)

    else:
        return 1j*(np.sign(k)*((1 + k**2)/(1 - k**2)) + k)

def a(k):
    return 1/(1j*l*omega_0 - B(k))

fig, ax = plt.subplots()

#ax.plot(gammas[0], gammas[1])
#ax.plot(la_0[0], la_0[1])

plt.show()
