import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from fourier import series, coeff, coeff_2d, array_integral


#A = np.array([[1, 2], [3, 4]])
#A_inv = LA.inv(A)
#v = np.array([1, 1])
#print(f"Test = {A_inv.dot(v)}")


## Parameters
eps = 0.01  # Epsilon in problem
N = 5  # Number of iterations of fixed point scheme
h = 0.1  # Differentiation step size


## x and t arrays
xs = np.arange(0, 2*np.pi + h, h)
ts = np.arange(0, 2*np.pi + h, h)


## Useful functions
def ei(x):
    return np.array([np.cos(x), np.sin(x)])

def perp(v):
    return np.array([-v[1], v[0]])

def integral(f, a, b):
    fs = f(xs)

    return array_integral(fs, xs)


x_coeff_res = 1
t_coeff_res = 1
num_x_coeffs = 2*x_coeff_res + 1
num_t_coeffs = 2*t_coeff_res + 1

x_shift = -x_coeff_res  # index_(x,py) + x_shift = index_(x,real)
t_shift = -t_coeff_res  # index_(t,py) + t_shift = index_(t,real)

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


## Definition of f and functions dependent on f
f_coeffs = np.zeros((num_x_coeffs, num_t_coeffs), dtype=complex)  # Initial condition f(x, t) = 0

def f(x, t):
    return series(f_coeffs, (x_shift, t_shift), x, t)

dt_f_coeffs = np.zeros_like(f_coeffs)
for (k, c) in np.ndenumerate(f_coeffs):
    dt_f_coeffs[*k] = 1j*(k[1] + t_shift)*c

def dt_f(x, t):
    return series(dt_f_coeffs, (x_shift, t_shift), x, t)

dx_f_coeffs = np.zeros_like(f_coeffs)
for (k, c) in np.ndenumerate(f_coeffs):
    dx_f_coeffs[*k] = 1j*(k[0] + x_shift)*c

def dx_f(x, t):
    return series(dx_f_coeffs, (x_shift, t_shift), x, t) 

def ga(x, t):
    return np.sqrt(1 + f(x, t)) * ei(x)

def dt_ga(x, t):
    return (dt_f(x, t)/(2*np.sqrt(1 + f(x, t)))) * ei(x)

def dx_ga(x, t):
    return (dx_f(x, t)/(2*np.sqrt(1 + f(x, t)))) * ei(x) + 1j*ga(x, t)


## Definition of lambda and functions dependent on lambda
h_0 = np.pi  # Initial height of rotating point

la_1_coeffs = np.zeros(num_t_coeffs, dtype=complex)
la_2_coeffs = np.zeros(num_t_coeffs, dtype=complex)

# Initial condition on lambda
la_1_coeffs[t_index_real_to_py(1)] = 1j * h_0/2
la_1_coeffs[t_index_real_to_py(-1)] = -1j * h_0/2
la_2_coeffs[t_index_real_to_py(1)] = h_0/2
la_2_coeffs[t_index_real_to_py(-1)] = h_0/2

def la_1(t):
    return series(la_1_coeffs, t_shift, t)

def la_2(t):
    return series(la_2_coeffs, t_shift, t)

def la(t):
    return np.array([la_1(t), la_2(t)])

dt_la_1_coeffs = np.array([1j*(k + t_shift)*c for (k, c) in enumerate(la_1_coeffs)])
dt_la_2_coeffs = np.array([1j*(k + t_shift)*c for (k, c) in enumerate(la_2_coeffs)])

def dt_la_1(t):
    return series(dt_la_1_coeffs, t_shift, t)

def dt_la_2(t):
    return series(dt_la_2_coeffs, t_shift, t)


## Definition of omega
omega_0 = 1/(2*h_0**2)  # Initial condition on omega 
omega = omega_0


## Definition of phi
phi_coeffs = np.zeros((num_x_coeffs, num_t_coeffs), dtype=complex)  # Initially phi(x, t) = 0

def phi(x, t):
    return series(f_coeffs, (x_shift, t_shift), x, t)

dt_phi_coeffs = np.array([1j*(k[1] + t_shift)*c for (k, c) in np.ndenumerate(f_coeffs)])

def dt_phi(x, t):
    return series(dt_phi_coeffs, (x_shift, t_shift), x, t)

dx_phi_coeffs = np.zeros_like(f_coeffs)
for (k, c) in np.ndenumerate(f_coeffs):
    dx_phi_coeffs[*k] = 1j*(k[0] + x_shift)*c

def dx_phi(x, t):
    return series(dx_phi_coeffs, (x_shift, t_shift), x, t)


## Definition of delta
# Initially delta_1(t) = delta_2(t) = 0
delta_1_coeffs = np.zeros(num_t_coeffs, dtype=complex)
delta_2_coeffs = np.zeros(num_t_coeffs, dtype=complex)

def delta_1(t):
    return series(delta_1_coeffs, t_shift, t)

def delta_2(t):
    return series(delta_2_coeffs, t_shift, t)

def delta(t):
    return np.array([delta_1(t), delta_2(t)])

dt_delta_1_coeffs = np.array([1j*(k + t_shift)*c for (k, c) in enumerate(delta_1_coeffs)])
dt_delta_2_coeffs = np.array([1j*(k + t_shift)*c for (k, c) in enumerate(delta_2_coeffs)])

def dt_delta_1(t):
    return series(dt_delta_1_coeffs, t_shift, t)

def dt_delta_2(t):
    return series(dt_delta_2_coeffs, t_shift, t)


## Y:s
def Y_1(x, t):
    def integrand(y):
            return np.log(LA.norm(ga(x, t) - ga(y, t))) * dx_ga(y, t).dot(perp(dx_ga(x, t)))

    #integrands = integrand(xs)

    X_1 = 1/np.pi * integral(integrand, 0, 2*np.pi) - (eps/np.pi)*(ga(x, t) - la(t)).dot(perp(dx_ga(x, t)))

    return omega*dt_f(x, t) - X_1

def Y_2(t):
    def integrand(y):
        return np.log(LA.norm(la(t) - ga(y, t))) * dx_ga(y, t)[0]

    X_2 = -1/(2*np.pi) * integral(integrand, 0, 2*np.pi)

    return omega*dt_la_1(t) - X_2

def Y_3(t):
    def integrand(y):
        return np.log(LA.norm(la(t) - ga(y, t))) * dx_ga(y, t)[1] 

    X_3 = -1/(2*np.pi) * integral(integrand, 0, 2*np.pi)

    return omega*dt_la_2(t) - X_3 

def Y_4():
    return la_1(0)


#print(f"Test = {coeff_2d(lambda x, y : np.exp(1j*(x + y)), 1, 1)}")

## Functions useful in inversion
def B(k):
    if k == 0:
        return 0

    elif abs(k) == 1:
        return 1j*np.sign(k)*(3/2)

    else:
        return 1j*(np.sign(k)*((1 + k**2)/(1 - k**2)) + k)

def a(k, l):
    print(f"Denominator equal to = {1j*l*omega_0 - B(k)}")
    return 1/(1j*l*omega_0 - B(k))


## Fixed point scheme
# Useful variables
nu = h_0/2

# Setting up H
H_1_coeffs = np.zeros_like(phi_coeffs)
H_2_coeffs = np.zeros_like(delta_1_coeffs)
H_3_coeffs = np.zeros_like(delta_2_coeffs)
H_4 = 0

# Setting up m
m_1_coeffs = np.zeros(num_t_coeffs)
m_2_coeffs = np.zeros(num_t_coeffs)

def m_1(t):
    def integrand(y):
        return (h_0 - np.sin(y))*phi(y, t)/(h_0**2 - 2*h_0*np.sin(y) + 1) 

    return 1/(4*np.pi) * integral(integrand, 0, 2*np.pi)

def m_2(t):
    def integrand(y):
        return np.cos(y)*phi(y, t)/(h_0**2 - 2*h_0*np.sin(y) + 1) 

    return 1/(4*np.pi) * integral(integrand, 0, 2*np.pi)

print(f"Test m_1 = {m_1(0)}")

def H_1(x, t):
    pass

for n in range(N):
    # Useful variables specific to iteration
    delta_1_coeff_sum = 0

    ## Computing H coeffs
    # Computing H_1 coeffs. Then calculating new phi coeffs.
    for I in np.ndindex(phi_coeffs.shape): 
        k = x_index_py_to_real(I[0])
        l = t_index_py_to_real(I[1])

        if (k, l) == (0, 0):
            continue

        H_1_coeffs[*I] = coeff_2d(Y_1, k, l)
        print(f"Loop 1 : Made it through iteration = {I} with a({k}, {l}) = {a(k, l)}")
        phi_coeffs[*I] = a(k, l) * H_1_coeffs[*I]

    # Computing H_2 and H_3 coeffs and H_4. Then calculating new delta_1 and delta_2 coeffs and new value for eta.
    for i in range(num_t_coeffs):
        print(f"Loop 2 : Made it through iteration = {i}")
        k = t_index_py_to_real(i)
        H_2_coeffs[i] = coeff(Y_2, k) 
        H_3_coeffs[i] = coeff(Y_3, k)
        H_4 = delta_1(0)

        m_1_coeffs[i] = coeff(m_1, k)
        m_2_coeffs[i] = coeff(m_2, k)

    # Calculating new delta_1 and delta_2 coeffs.
    for i in range(num_t_coeffs):
        print(f"Loop 3 : Made it through iteration = {i}")
        k = t_index_py_to_real(i)

        if abs(k) == 1:
            continue

        M = np.array([[1j*omega_0*k, -omega_0], [-omega_0, 1j*omega_0*k]])
        M_inv = LA.inv(M)
        
        v = np.array([H_2_coeffs[i] - m_1_coeffs[i], H_3_coeffs[i] - m_2_coeffs[i]])
        
        w = M_inv.dot(v)

        print(f"i = {i}")

        delta_1_coeffs[i] = w[0]
        delta_2_coeffs[i] = w[1]
        delta_1_coeff_sum += delta_1_coeffs[i]

    # Calculating delta_1 and delta_2 coeffs for k = -1, +1 and new value for eta.
    i_plus1 = t_index_real_to_py(1)
    i_minus1 = t_index_real_to_py(-1)
    
    M = np.array([[-nu, -1j*omega_0, -omega_0, 0], [nu, -1j*omega_0, 0, -omega_0], [1j*nu, omega_0, 1j*omega_0, 0], [-1j*nu, -omega_0, 0, -1j*omega_0]])
    M_inv = LA.inv(M)
    
    v = np.array([H_2_coeffs[i_plus1] - m_1_coeffs[i_plus1] + 1j*omega_0*(delta_1_coeff_sum + H_4), H_2_coeffs[i_minus1] - m_1_coeffs[i_minus1], H_3_coeffs[i_plus1] - m_2_coeffs[i_plus1] - omega_0*(delta_1_coeff_sum + H_4), H_3_coeffs[i_minus1] - m_2_coeffs[i_minus1]])

    w = M_inv.dot(v)
    
    eta = w[0]
    delta_1_coeffs[i_minus1] = w[1]
    delta_2_coeffs[i_plus1] = w[2]
    delta_2_coeffs[i_minus1] = w[3]

    f_coeffs = f_coeffs + phi_coeffs
    la_1_coeffs = la_1_coeffs + delta_1_coeffs
    la_2_coeffs = la_2_coeffs + delta_2_coeffs
    omega = omega + eta













print(f"lambda_1_coeffs = {la_1_coeffs}\n\nlambda_2_coeffs = {la_2_coeffs}")

#print(dx_ga(1, 1))
#print(coeff_2d(Y_1, 0, 0))

ts = np.linspace(0, 2*np.pi, 100)
la_1s = la_1(ts)
la_2s = la_2(ts)
#gas = ga(ts, 0)

fig, ax = plt.subplots()

ax.plot(la_1s, la_2s)
#ax.plot(gas[0], gas[1])

plt.show()
