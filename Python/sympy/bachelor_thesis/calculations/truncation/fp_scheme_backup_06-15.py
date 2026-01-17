import time
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, Slider


START = time.time()  # For measuring calculation runtime


################
#  PARAMETERS  #
################

## Theoretical parameters ##
eps = 0  # Epsilon in problem

## Numerical parameters ##
N = 1  # Number of iterations of fixed point scheme
n = 100 + 1  # Number of points per axis in grid
h = 2*np.pi/(n - 1)  # Step size

x_coeff_res = 2**4  # Resolution of truncation in x (2*x_coeff_res + 1 = nr of x coeffs)
t_coeff_res = 2**4  # Resolution of truncation in t (2*t_coeff_res + 1 = nr of t coeffs)
num_x_coeffs = 2*x_coeff_res + 1
num_t_coeffs = 2*t_coeff_res + 1

x_shift = -x_coeff_res  # index_(x,py) + x_shift = index_(x,real)
t_shift = -t_coeff_res  # index_(t,py) + t_shift = index_(t,real)



##################################
#  USEFUL QUANTITIES AND ARRAYS  #
##################################

## x and t arrays ##
x = np.linspace(0, 2*np.pi, n, dtype=np.float128)
t = np.linspace(0, 2*np.pi, n, dtype=np.float128)
X, T = np.meshgrid(x, t, indexing="ij")
Y, _, __ = np.meshgrid(x, x, t, indexing="ij")  # Put Y in front for simple numpy broadcasting
one = 0*Y + 1  # For stretching arrays from X x T domain into Y x X x T domain

ei = np.array([np.cos(X), np.sin(X)])  # e^(i*X) = (cos(X), sin(X))

h_0 = 4*np.sqrt(2)  # Initial height of rotating point

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

    return res

## Fourier analysis functions ##
k = np.linspace(-x_coeff_res, x_coeff_res, 2*x_coeff_res + 1, dtype=int)
l = np.linspace(-t_coeff_res, t_coeff_res, 2*t_coeff_res + 1, dtype=int)
K, L, XX, TT = np.meshgrid(k, l, x, t, indexing="ij")
series_exponentials_2d = np.exp(1j*(K*XX + L*TT))  # Array containing exponentials e^(i(kX + lT))
i_0 = x_index_real_to_py(0)
series_exponentials_1d = series_exponentials_2d[i_0,:,:,:]
coeff_exponentials_2d = np.flip(series_exponentials_2d, (0, 1))
coeff_exponentials_1d = np.flip(series_exponentials_1d, 0)


def series(coeffs, no_x=False):
    print("Started a series calculation") 
    start = time.time()
    result = np.zeros_like(X, dtype=np.float128)

    if no_x:
        result = np.sum((series_exponentials_1d.T*coeffs).T, axis=0)

        return result

    else:
        result = np.sum(np.sum((series_exponentials_2d.T*coeffs.T).T, axis=0), axis=0)

        return result

    end = time.time()
    print(f"Ended a series calculation. Took = {end - start} s")

    return result


def coeff(f, k, l, no_x=False):
    xt_integrand = f * np.exp(-1j*(k*X + l*T))

    x_integrand = riemann(xt_integrand, axis=1)

    if no_x:
        intg = 1/(2*np.pi) * x_integrand[0]

        return intg

    else:
        intg = 1/(4*np.pi**2) * riemann(x_integrand)

        return intg


def old_series(coeffs, no_x=False):
    print("Started a series calculation") 
    start = time.time()
    result = np.zeros_like(X, dtype=np.float128)

    if no_x:
        for (j, c) in enumerate(coeffs):
            l = t_index_py_to_real(j)
            
            result += (c*np.exp(1j*l*T)).real

    else:
        shape = coeffs.shape

        for i in range(shape[0]):
            for j in range(shape[1]):
                k = x_index_py_to_real(i)
                l = t_index_py_to_real(j)

                c = coeffs[i, j]

                result += (c*np.exp(1j*(k*X + l*T))).real

    end = time.time()
    print(f"Ended a series calculation. Took = {end - start} s")

    return result


## Functions useful in inversion ##
def B(k):
    if k == 0:
        return 0

    elif abs(k) == 1:
        return 1j*np.sign(k)*(3/2)

    else:
        return 1j*(np.sign(k)*((1 + k**2)/(1 - k**2)) + k)

def a(k, l):
    return 1/(1j*l*omega_0 - B(k))



###############################################
#  SETUP OF COEFFICIENTS ARRAYS OF FUNCTIONS  #
###############################################

## f setup ##
f_coeffs = np.zeros((num_x_coeffs, num_t_coeffs), dtype=np.complex256)  # Initial condition f(x, t) = 0
dt_f_coeffs = np.zeros_like(f_coeffs, dtype=np.complex256)
dx_f_coeffs = np.zeros_like(f_coeffs, dtype=np.complex256)


## lambda setup ##
la_1_coeffs = np.zeros(num_t_coeffs, dtype=np.complex256)
la_2_coeffs = np.zeros(num_t_coeffs, dtype=np.complex256)

# Initial condition on lambda
la_1_coeffs[t_index_real_to_py(1)] = 1j * h_0/2
la_1_coeffs[t_index_real_to_py(-1)] = -1j * h_0/2
la_2_coeffs[t_index_real_to_py(1)] = h_0/2
la_2_coeffs[t_index_real_to_py(-1)] = h_0/2

dt_la_1_coeffs = np.zeros_like(la_1_coeffs, dtype=np.complex256)
dt_la_2_coeffs = np.zeros_like(la_2_coeffs, dtype=np.complex256)


## omega setup ##
omega_0 = 1/(2*h_0**2) 
omega = omega_0


## phi setup ##
phi_coeffs = np.zeros((num_x_coeffs, num_t_coeffs), dtype=np.complex256)  # Initially phi(x, t) = 0
dt_phi_coeffs = np.zeros_like(f_coeffs)
dx_phi_coeffs = np.zeros_like(f_coeffs)


## delta setup ##
# Initially delta_1(t) = delta_2(t) = 0
delta_1_coeffs = np.zeros(num_t_coeffs, dtype=np.complex256)
delta_2_coeffs = np.zeros(num_t_coeffs, dtype=np.complex256)
dt_delta_1_coeffs = np.zeros_like(delta_1_coeffs, dtype=np.complex256)
dt_delta_2_coeffs = np.zeros_like(delta_2_coeffs, dtype=np.complex256)


## m setup ##
m_1_coeffs = np.zeros(num_t_coeffs, dtype=np.complex256)
m_2_coeffs = np.zeros(num_t_coeffs, dtype=np.complex256)


## H setup ##
H_1_coeffs = np.zeros_like(phi_coeffs, dtype=np.complex256)
H_2_coeffs = np.zeros_like(delta_1_coeffs, dtype=np.complex256)
H_3_coeffs = np.zeros_like(delta_2_coeffs, dtype=np.complex256)
H_4 = 0



######################
#  UPDATE FUNCTIONS  #
######################

## Update f and gamma ##
def update_f():
    print("Started f update") 
    start = time.time()
    global f

    # Update f array from updated f coeffs  
    f = series(f_coeffs)
    #print(f"f = {f}\n")

    # Update functions dependent on f
    update_dt_f()
    update_dx_f()
    update_ga()
    end = time.time()
    print(f"Ended f update. Took = {end - start} s")

def update_dt_f():
    global dt_f

    # Update dt_f coeffs from updated f coeffs
    for (k, c) in np.ndenumerate(f_coeffs):
        dt_f_coeffs[*k] = 1j*(k[1] + t_shift)*c

    # Update dt_f array from updated f array 
    dt_f = np.gradient(f, t, axis=1, edge_order=2)
    #print(f"dt_f = {dt_f}\n")

def update_dx_f():
    global dx_f

    # Update dx_f coeffs from updated f coeffs
    for (k, c) in np.ndenumerate(f_coeffs):
        dx_f_coeffs[*k] = 1j*(k[0] + x_shift)*c

    # Update dx_f array from updated f array 
    dx_f = np.gradient(f, x, axis=0, edge_order=2) 
    #print(f"dx_f = {dx_f}\n")

def update_ga():
    global ga

    # Update gamma array from updated f array 
    ga = np.sqrt(1 + f) * ei
    #print(f"ga = {ga}\n")

    # Update functions dependent on gamma
    update_dt_ga()
    update_dx_ga()

def update_dt_ga():
    global dt_ga

    # Update dt_gamma array from updated f and dt_f arrays 
    dt_ga = (dt_f/(2*np.sqrt(1 + f))) * ei
    #print(f"dt_ga = {dt_ga}\n")

def update_dx_ga():
    global dx_ga

    # Update dx_gamma array from updated f, dx_f and gamma arrays 
    dx_ga = (dx_f/(2*np.sqrt(1 + f))) * ei + perp(ga)
    #print(f"dx_ga = {dx_ga}\n")


## Update lambda ##
def update_la():
    print("Started lambda update") 
    start = time.time()
    global la_1, la_2, la

    # Update lambda arrays from updated lambda coeffs 
    la_1 = series(la_1_coeffs, no_x=True)
    la_2 = series(la_2_coeffs, no_x=True)
    la = np.array([la_1, la_2])
    #print(f"la = {la}\n")

    # Update functions dependent on lambda
    update_dt_la()
    end = time.time()
    print(f"Ended lambda update. Took = {end - start} s")

def update_dt_la():
    global dt_la_1_coeffs, dt_la_2_coeffs, dt_la_1, dt_la_2, dt_la

    # Update dt_lambda coeffs from updated lambda coeffs 
    dt_la_1_coeffs = np.array([1j*(k + t_shift)*c for (k, c) in enumerate(la_1_coeffs)])
    dt_la_2_coeffs = np.array([1j*(k + t_shift)*c for (k, c) in enumerate(la_2_coeffs)])

    # Update dt_lambda arrays from updated dt_lambda coeffs
    dt_la_1 = series(dt_la_1_coeffs, no_x=True)
    dt_la_2 = series(dt_la_2_coeffs, no_x=True)
    dt_la = np.array([dt_la_1, dt_la_2])
    #print(f"dt_la = {dt_la}\n")


## Update phi, m_1 and m_2 ##
def update_phi():
    print("Started phi update") 
    start = time.time()
    global phi, phi_x, phi_y

    phi = series(phi_coeffs)
    #print(f"phi = {phi}\n")

    update_dt_phi()
    update_dx_phi()

    # Setting up integrable phi for m_1 and m_2 
    phi_x = phi*one
    phi_y = np.swapaxes(phi_x, 0, 1)

    # Update m_1 and m_2
    update_m_1()
    update_m_2()
    end = time.time()
    print(f"Ended phi update. Took = {end - start} s")

def update_dt_phi():
    global dt_phi

    # Update dt_phi coeffs from updated phi coeffs
    for (k, c) in np.ndenumerate(phi_coeffs):
        dt_phi_coeffs[*k] = 1j*(k[1] + t_shift)*c

    # Update dt_phi array from updated dt_phi coeffs
    dt_phi = series(dt_phi_coeffs)
    #print(f"dt_phi = {dt_phi}\n")

def update_dx_phi():
    global dx_phi

    # Update dx_phi coeffs from updated phi coeffs
    for (k, c) in np.ndenumerate(f_coeffs):
        dx_phi_coeffs[*k] = 1j*(k[0] + x_shift)*c

    # Update dx_phi array from updated dx_phi coeffs
    dx_phi = series(dx_phi_coeffs)
    #print(f"dx_phi = {dx_phi}\n")

def update_m_1():
    print("Started m_1 update")
    global m_1

    m_1_integrand = (h_0 - np.sin(Y))*phi_y/(h_0**2 - 2*h_0*np.sin(Y) + 1)
    m_1 = 1/(4*np.pi) * riemann(m_1_integrand, axis=0) 
    #print(f"m_1 = {m_1}\n")
    print("Ended m_1 update")

def update_m_2():
    print("Started m_2 update")
    global m_2

    m_2_integrand = np.cos(Y)*phi_y/(h_0**2 - 2*h_0*np.sin(Y) + 1)
    m_2 = 1/(4*np.pi) * riemann(m_2_integrand, axis=0)
    #print(f"m_2 = {m_2}\n")
    print("Ended m_2 update")


## Update delta ##
def update_delta():
    print("Started delta update")
    global delta_1, delta_2, delta

    delta_1 = series(delta_1_coeffs, no_x=True)
    delta_2 = series(delta_2_coeffs, no_x=True)

    delta = np.array([delta_1, delta_2])
    #print(f"delta_1 = {delta_1}\n\ndelta_2 = {delta_2}\n")

    # Update functions dependent on delta
    update_dt_delta()
 
    print("Ended delta update")

def update_dt_delta():
    global dt_delta_1_coeffs, dt_delta_2_coeffs, dt_delta_1, dt_delta_2, dt_delta

    # Update dt_delta coeffs from updated delta coeffs
    dt_delta_1_coeffs = np.array([1j*(k + t_shift)*c for (k, c) in enumerate(delta_1_coeffs)])
    dt_delta_2_coeffs = np.array([1j*(k + t_shift)*c for (k, c) in enumerate(delta_2_coeffs)])

    # Update dt_delta arrays from updated dt_delta coeffs
    dt_delta_1 = series(dt_delta_1_coeffs, no_x=True)
    dt_delta_2 = series(dt_delta_2_coeffs, no_x=True)
    dt_delta = np.array([dt_delta_1, dt_delta_2])
    #print(f"dt_delta = {dt_delta}\n") 


## Update Ys ##
def update_Ys():
    print("Started Ys update")
    global X_1, X_2, X_3, Y_1, Y_2, Y_3, Y_4

    ga_1_x = ga[0]*one
    ga_2_x = ga[1]*one
    ga_x = np.array([ga_1_x, ga_2_x])
    ga_y = np.swapaxes(ga_x, 1, 2)
    dx_ga_1_x = dx_ga[0]*one
    dx_ga_2_x = dx_ga[1]*one
    dx_ga_x = np.array([dx_ga_1_x, dx_ga_2_x])
    dx_ga_y = np.swapaxes(dx_ga_x, 1, 2)

    la_1_integrable = one*la[0]
    la_2_integrable = one*la[1]
    la_integrable = np.array([la_1_integrable, la_2_integrable])
 
    X_1_integrand = np.nan_to_num(np.log(norm(ga_x - ga_y)) * dot(dx_ga_y, perp(dx_ga_x)))  # Collapses the first axis (axis specifying that gamma is a vector). Remaining shape: Y x X x T.
    X_1 = 1/np.pi * riemann(X_1_integrand, axis=0) - (eps/np.pi)*dot(ga - la, perp(dx_ga))/(norm(ga - la)**2)  # Integral performed along first axis (Y-axis) thus producing an array of shape X x T, then epsilon part is added on.
    Y_1 = omega*dt_f - X_1

    X_2_integrand = np.log(norm(la_integrable - ga_y)) * dx_ga_y[0]
    X_2 = -1/(2*np.pi) * riemann(X_2_integrand, axis=0)
    Y_2 = omega*dt_la[0] - X_2

    X_3_integrand = np.log(norm(la_integrable - ga_y)) * dx_ga_y[1]
    X_3 = -1/(2*np.pi) * riemann(X_3_integrand, axis=0)
    Y_3 = omega*dt_la[1] - X_3

    #print(f"X_2 = {X_2}\n\nX_3 = {X_3}\n")

    Y_4 = la[0][0][0]  # Y_4 = la_1(0)

    #print(f"X_1_integrand = {X_1_integrand}\n\nX_1 = {X_1}\n\nX_2 = {X_2}\n\nX_3 = {X_3}\n\nY_1 = {Y_1}\n\nY_2 = {Y_2}\n\nY_3 = {Y_3}\n\nY_4 = {Y_4}\n")

    print("Ended Ys update")


## Main update schemes ##
def update_1():
    update_phi()

def update_2():
    update_delta() 
    update_f()  # Updates gamma as well
    update_la()
    update_Ys()



####################
#  INITIALIZATION  #
####################

update_f()  # Initialize f
update_la()  # Initialize lambda
update_phi()  # Initialize phi
update_delta()  # Initialize delta
update_Ys()  # Initialize Ys



########################
#  FIXED POINT SCHEME  #
########################

for n in range(N):
    print(f"Started {n}:th iteration of main loop.")
    # Useful variables specific to iteration
    delta_1_coeff_sum = 0

    ## Computing H coeffs ##
    # Computing H_1 coeffs. Then calculating new phi coeffs.
    for I in np.ndindex(phi_coeffs.shape): 
        k = x_index_py_to_real(I[0])  # Translates from python index I[0] = i to Fourier index k (index corresponding to x)
        l = t_index_py_to_real(I[1])  # Translates from python index I[1] = j to Fourier index l (index corresponding to l)

        if (k, l) == (0, 0):
            continue

        H_1_coeffs[*I] = coeff(Y_1, k, l)
        #print(f"Loop 1 : Made it through iteration = {I} with a({k}, {l}) = {a(k, l)}")
        phi_coeffs[*I] = a(k, l) * H_1_coeffs[*I]

    #print(f"H_1_coeffs = {H_1_coeffs}\n")
    update_1()  # Update phi array after changing its coefficients (also updates the arrays of functions dependent on phi, i.e. m_1 and m_2)

    ## Computing H_2 and H_3 coeffs and H_4. Then calculating new delta_1 and delta_2 coeffs and new value for eta. ##
    for j in range(num_t_coeffs):
        #print(f"Loop 2 : Made it through iteration = {j}")
        l = t_index_py_to_real(j)  # Translates from python index j to Fourier index l (index corresponding to t)
        H_2_coeffs[j] = coeff(Y_2, 0, l, no_x=True)
        H_3_coeffs[j] = coeff(Y_3, 0, l, no_x=True)
        H_4 = delta[0][0][0]  # H_4 = delta_1(0)

        m_1_coeffs[j] = coeff(m_1, 0, l, no_x=True)
        m_2_coeffs[j] = coeff(m_2, 0, l, no_x=True)

    #print(f"H_2 coeffs = {H_2_coeffs}\n\nH_3 coeffs = {H_3_coeffs}\n")

    ## Calculating new delta_1 and delta_2 coeffs for k != -1, +1. ##
    for j in range(num_t_coeffs):
        #print(f"Loop 3 : Made it through iteration = {j}")
        l = t_index_py_to_real(j)  # Translates from python index j to Fourier index l (index corresponding to t)

        # Skips l = -1, 1 as these cases are dealt with after having calculated all other Fourier coefficients of delta.
        if abs(l) == 1:
            continue

        Mat = np.array([[1j*omega_0*l, M], [M, 1j*omega_0*l]])
        Mat_inv = LA.inv(Mat)
        
        v = np.array([H_2_coeffs[j] - m_1_coeffs[j], H_3_coeffs[j] - m_2_coeffs[j]])
        
        w = Mat_inv.dot(v)

        delta_1_coeffs[j] = w[0]
        delta_2_coeffs[j] = w[1]
        delta_1_coeff_sum += delta_1_coeffs[j] 

    # Calculating delta_1 and delta_2 coeffs for l = -1, +1 and new value for eta.
    j_plus1 = t_index_real_to_py(1)
    j_minus1 = t_index_real_to_py(-1)
    
    Mat = np.array([[-nu, -1j*omega_0, M, 0], [nu, -1j*omega_0, 0, M], [1j*nu, -M, 1j*omega_0, 0], [-1j*nu, M, 0, -1j*omega_0]])
    Mat_inv = LA.inv(Mat)
    
    v = np.array([H_2_coeffs[j_plus1] - m_1_coeffs[j_plus1] + 1j*omega_0*(delta_1_coeff_sum + H_4), H_2_coeffs[j_minus1] - m_1_coeffs[j_minus1], H_3_coeffs[j_plus1] - m_2_coeffs[j_plus1] - omega_0*(delta_1_coeff_sum + H_4), H_3_coeffs[j_minus1] - m_2_coeffs[j_minus1]])

    w = Mat_inv.dot(v)
    
    eta = w[0]
    delta_1_coeffs[j_minus1] = w[1]
    delta_1_coeffs[j_plus1] = -delta_1_coeffs[j_minus1] - delta_1_coeff_sum + H_4
    delta_2_coeffs[j_plus1] = w[2]
    delta_2_coeffs[j_minus1] = w[3]

    #print(f"v = {v}\n\nw = {w}\n")

    j = t_index_real_to_py(0)

    #print(f"delta_2 coeffs = {delta_2_coeffs}\n\ndelta_2 coeff {j} = {delta_2_coeffs[j]}\n") 

    f_coeffs += phi_coeffs
    la_1_coeffs += delta_1_coeffs
    la_2_coeffs += delta_2_coeffs
    omega += eta

    update_2()  # Update delta, f, lambda and Ys after having updated relevant coefficients

    print(f"Ended {n}:th iteration of main loop.")


STOP = time.time()
print(f"Calculation ended. Runtime = {STOP - START} s\n")



#print(f"lambda_1_coeffs = {la_1_coeffs}\n\nlambda_2_coeffs = {la_2_coeffs}")

print(f"la_1 = {la_1}\n\nla_2 = {la_2}\n\nomega = {omega}\n\nT = {2*np.pi/omega}")

fig, ax = plt.subplots()

fig.subplots_adjust(bottom=0.25)

ax.plot(la_1[0,:], la_2[0,:])
ax.set_aspect("equal")

line, = ax.plot(ga[0][:,0], ga[1][:,0])
point, = ax.plot([la_1[0,0]], [la_2[0,0]], "bo")

axtime = fig.add_axes([0.25, 0.1, 0.65, 0.03])
time_slider = Slider(
    ax=axtime,
    label='Time [theta]',
    valmin=0,
    valmax=t.size - 1,
    valinit=0,
    valstep=1
)

def update(val):
    i = int(time_slider.val)
    #print(f"Houston we have a change!\n\ni = {i}\n\nga_1 = {ga[0,:,i]}\n\nga_2 = {ga[1,:,i]}\n") 
    line.set_data(ga[0][:,i], ga[1][:,i])
    point.set_data([la_1[0,i]], [la_2[0,i]])
    ax.relim()
    ax.autoscale_view()
    fig.canvas.draw_idle()

time_slider.on_changed(update)

fig2, ax2 = plt.subplots()

#ax2.plot(t, X_3[0], label="X_3")
#ax2.plot(t, omega*dt_la_2[0], label="la_2")
ax2.plot(t, Y_3[0], label="Y_3")
#ax2.plot(t, series(H_2_coeffs, no_x=True)[0], label="series")

ax2.legend()

fig3 = plt.figure()
ax3 = plt.axes(projection='3d')

ax3.plot(X, T, Y_1, label="Y_1")

ax3.legend()

plt.show()
