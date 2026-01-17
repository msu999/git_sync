import numpy as np

def simpsons_single(fs, ts, a, b):
    N = ts.size
    
    n_a = max([i for i in range(N) if ts[i] <= a])
    n_b = min([i for i in range(N) if ts[i] >= b]) 

    I = 0
    n = n_a

    while n <= n_b - 2:
        cur_h = (ts[n + 2] - ts[n])/2
        I += (cur_h/3)*(fs[n] + 4*fs[n + 1] + fs[n + 2])

        n += 2

    # If simpsons can't be performed on the last step then perform trapezoidal rule for this one
    if n != n_b:
        I += ((ts[n_b] - ts[n_b - 1])/2)*(fs[n_b - 1] + fs[n_b])

    return I

def simpsons_total(fs, ts):
    Fs = np.array([simpsons_single(fs, ts, ts[0], b) for b in ts.tolist()])

    return Fs
