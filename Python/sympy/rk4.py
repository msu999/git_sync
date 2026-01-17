import numpy as np
from numnal import one

def rk4(f, y_0, t_0, t_1, N, iscomplex=False): 
    ts = np.linspace(t_0, t_1, N)

    if iscomplex == True:
        ys = np.zeros((N, y_0.size), dtype=complex)

    else:
        ys = np.zeros((N, y_0.size))
    
    h = (t_1 - t_0)/N

    ys[0] = y_0

    for n in range(N-1):
        k_1 = f(ts[n], ys[n])
        k_2 = f(ts[n] + h/2, ys[n] + h*k_1/2)
        k_3 = f(ts[n] + h/2, ys[n] + h*k_2/2)
        k_4 = f(ts[n] + h, ys[n] + h*k_3)

        ys[n + 1] = ys[n] + (h/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)

    return (ts, np.transpose(ys))

def simpsons(fs, xs, axis=0):
    # Integrates along x_"axis" (axis number "axis")

    xs_shape = xs.shape
    print(f"xs shape = {xs_shape}", end="\n\n")
    num_vars = len(xs_shape)

    Fs = np.empty_like(fs)
    print(f"First index = {1, *np.zeros(num_vars - 1, dtype=int)}, Second index = {np.zeros(num_vars, dtype=int)}\n\nFirst x = {xs[1, *np.zeros(num_vars - 1, dtype=int)]}, Second x = {xs[*np.zeros(num_vars, dtype=int)]}", end="\n\n")

    dx = xs[1, *np.zeros(num_vars - 1, dtype=int)] - xs[*np.zeros(num_vars, dtype=int)]   # Assumes constant step size in all directions
    
    if type(dx) != np.float64:
        dx = dx[0]

    print(f"dx = {dx}", end="\n\n")

#    index_ranges = list(xs_shape)
#    index_ranges[axis] = 1  # To guarantee no variation along axis "axis" in outermost loop
#
#    for I in np.ndindex(*index_ranges):
#        for i in range(xs_shape[axis]):
#            mutable_I = list(I)
#            mutable_I[axis] = i
#            J_ = 

    dI = one(num_vars, axis)
    print(f"dI = {dI}", end="\n\n")
 
    for I in np.ndindex(xs_shape):
        #print(f"I = {I}, I - 2*dI = {I - 2*dI}, I - dI = {I - dI}, I + dI = {I + dI}\n\nFs = {Fs}", end="\n\n")
        if I[axis] == 0:
            Fs[*I] = 0*fs[*I]

        elif I[axis] == 1:
            Fs[*I] = dx/2 * (fs[*(I - dI)] + fs[*I])

        elif I[axis] == xs_shape[axis] - 1:
            Fs[*I] = Fs[*(I - dI)] + dx/2 * (fs[*(I - dI)] + fs[*I])

        else: 
            dF = dx/3 * (fs[*(I - 2*dI)] + 4*fs[*(I - dI)] + fs[*I])
            #print(f"Fs[I - 2dI] = {Fs[*(I-2*dI)]}", end="\n\n")
            Fs[*I] = Fs[*(I-2*dI)] + dF

    return Fs
