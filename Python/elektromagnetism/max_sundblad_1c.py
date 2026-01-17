import numpy as np
import matplotlib.pyplot as plt

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

# Capacitance in [F]
sample_Cs = np.array([0.432 , 0.232 , 0.165 , 0.129 , 0.107 , 0.094 , 0.083 , 0.078])*(10**(-9))

# Distance between capacitor plates in [m]
sample_ds = np.array([2 , 4 , 6 , 8 , 10 , 12 , 14 , 15])*(10**(-3))

sample_xs = 1/sample_ds
sample_ys = sample_Cs

res = reg(sample_xs, sample_ys)

k = res[0][0]
m = res[0][1]
u_k = res[1][0]
u_m = res[1][1]
u_D = 5e-4
D = 0.20  # Diameter of plate
A = np.pi*((D/2)**2)
eps0 = 8.85418782e-12
epsr = k/(eps0*A)
u_A = u_D*np.pi*D/2
u_epsr = np.sqrt((u_k * (1 / (eps0 * A)))**2 + (u_A * (k / (eps0 * A**2)))**2) 

print(f"k = {k}, m = {m}, u_k = {u_k}, u_m = {u_m}\nD = {D}, A = {A}, eps0 = {eps0}, epsr = {epsr}, u_A = {u_A}, u_D = {u_D}, u_epsr = {u_epsr}")

fig, ax = plt.subplots()

ax.plot(sample_xs, sample_ys, "*")

us = np.linspace(0, 1000, 1000)
vs = k*us+m

ax.plot(us, vs)

ax.set(xlabel='', ylabel='',
       title='')
ax.grid()

fig.savefig("plot.png")
plt.show()
