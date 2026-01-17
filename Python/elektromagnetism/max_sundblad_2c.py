import numpy as np

def reg(x, y):
    sx = np.sum(x)
    sy = np.sum(y)
    sx2 = np.linalg.norm(x)**2
    sxy = np.dot(x, y)
    n = x.size

    det = n*sx2 - sx**2

    A = (1/det)*np.array([[n, -sx], [-sx, sx2]])
    v = np.array([sxy, sy])

    params = np.matmul(A, v)
    
    S = np.sum((params[0]*x + params[1] - y)**2)

    u_params = np.array([np.sqrt((A[k][k]*S)/(n - 2)) for k in range(2)])

    return (params, u_params)

sample_Us = np.array([100, 90, 80, 70, 60, 50, 40, 30, 20, 10])*(10**(-3))
sample_ts = np.array([0*60+00.0, 0*60+11.7, 0*60+23.2, 0*60+36.8, 0*60+52.5, 1*60+10.9, 1*60+33.9, 2*60+03.1, 2*60+44.5, 3*60+54.9])

sample_xs = sample_ts
sample_ys = np.log(sample_Us)

res = reg(sample_xs, sample_ys)

k = res[0][0]
m = res[0][1]
tau = -1/k
U_0 = np.exp(m)
u_k = res[1][0]
u_m = res[1][1]
u_tau = u_k/(k**2)
u_U0 = u_m*U_0

print(f"k = {k}, m = {m}, u_k = {u_k}, u_m = {u_m}\ntau = {tau}, U_0 = {U_0}, u_tau = {u_tau}, u_U0 = {u_U0}")
