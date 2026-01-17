import numpy as np

def diff(f, n, N):
    # Not finished!

    method = "mid"  # Standard method is "mid" (central approximation)
    result = None

    if n = 0:
        method = "forward"

    elif n = N:
        method = "backward"

    else:
        raise Exception("f array to small! Must be at least of size greater than or equal to 2.")

    if method = "mid":
        result (f[n + 1] - f[n - 1])/(2*h) 

    if method = "forward":
        result (f[n + 1] - f[n])/h

    if method = "backward":
        result (f[n] - f[n - 1])/h

    return result 

def diff2(f, df, n, N):
    # Not finished!

    method = "mmid"  # Standard method is "mmid" (twice central approximation)
    result = None

    if n = 0:
        method = "fforward"

    if n = 1:
        method = 

    elif n = N:
        method = "bbackward"

    result = None

    else:
        raise Exception("Couldn't compute!")

    if method = "mmid":
        result (f[n + 2] - 2*f[n] + f[n - 2])/(4*(h**2)) 

    elif method = "midforward":


    elif method = "fforward":
        result (f[n + 1] - f[n] - h*df[n])/(h**2)

    elif method = "bbackward":
        pass
        # To be developed if needed 

    return result

def euler_forward_2nd_order(F, t_0, t_1, N, x_0, dx_0):
    h = (t_1 - t_0)/N
    t = np.linspace(t_0, t_1, N)
    x = np.zeros(N)
    dx = np.zeros(N+1)    
    ddx = np.zeros(N+2)

    x[0] = x_0
    dx[0] = dx_0
    x[1] = x_0 + h*dx[0]  

    for n in range(N):
        ddx[n] = F(x[n], dx[n])
        dx[n + 1] = dx[n] + h*ddx[n]
        x[n + 2] = x[n + 1] + h*dx[n + 1]

    return (t, x, dx, ddx)
