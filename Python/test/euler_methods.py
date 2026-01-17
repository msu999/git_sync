import numpy as np

def euler_forward(f, t_0, t_1, y_0, h):
    # f should take arguments of the form (t, y)

    ts = np.arange(t_0, t_1 + h, h)

    if ts[-1] != t_1:
        ts = np.append(ts, t_1)
    
    N = ts.size - 1  # Last index of t array
    
    N_funcs = y_0.size  # Number of functions in system of ODE:s

    ys = np.zeros((N + 1, N_funcs))
    ys[0] = y_0

    # Finds y:s for t:s in {t_0, t_0 + h, ... , t_1 - h}
    for i in range(N - 1):
        ys[i + 1] = ys[i] + f(ts[i], ys[i])

    # Finds last y (for t = t_1)
    ys[N] = ys[N - 1] + f(ts[N - 1], ys[N - 1])

    ys = np.transpose(ys)

    return (ts, ys)
