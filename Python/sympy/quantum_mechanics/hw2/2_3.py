import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

x_0 = 0
v_0 = 1
k_0 = 1

def f_0(x):
    return ((2/np.pi)**(1/4)) * np.exp(-x**2 + 1j*k_0*x)

xs = np.linspace(-2, 2, 1000)
print(np.trapz(f_0(xs)*f_0(xs).conjugate(), x=xs))

def v(x):
    return 0  # v_0*(x**2 - 1)

a = 2j

def b(x):
    return -2*v(x)

def ddx_f(i_t, i_x):
    if i_x == 0:
        return (fs[i_t, i_x + 2*1] - 2*fs[i_t, i_x + 1] + fs[i_t, i_x])/(2*(dx**2))

    elif i_x == N_x - 1:
        return (fs[i_t, i_x - 2*1] - 2*fs[i_t, i_x - 1] + fs[i_t, i_x])/(2*(dx**2))

    else:
        return (fs[i_t, i_x + 1] - 2*fs[i_t, i_x] + fs[i_t, i_x - 1])/(dx**2)

def t_step(i_t, i_x):
    if i_t == 0:
        return -(dt/a)*ddx_f(i_t, i_x) - (b(x)*dt/a)*fs[i_t, i_x] + fs[0, i_x]

    else:
        return -(2*dt/a)*ddx_f(i_t, i_x) - (2*b(x)*dt/a)*fs[i_t, i_x] + fs[i_t - 1, i_x]


N_t = 1000
N_x = 1000

T = 10
X = 10

dt = T/N_t
dx = X/N_x

ts = np.linspace(0, T, N_t)
xs = np.linspace(-X, X, N_x)

#ts, xs = np.meshgrid(t, x)
fs = np.zeros((ts.size, xs.size), dtype=complex)

for i_t, t in enumerate(ts):
    if i_t == ts.size - 1:
        continue
    for i_x, x in enumerate(xs):
        if i_t == 0:
            fs[0, i_x] = f_0(x)

        if i_x == xs.size - 1:
            continue
        fs[i_t + 1, i_x] = t_step(i_t, i_x)

#print(fs)

rho = fs*fs.conjugate()

fig, ax = plt.subplots()


# Create a slider from 0.0 to 1.0 in axes ax_t
# with 0.6 as initial value.
ax_slider = plt.axes([0.2, 0.1, 0.65, 0.05])  # [left, bottom, width, height]
t_slider = Slider(ax_slider, 't index', 0, 100, valinit=0, valstep=1)

line, = ax.plot(xs, rho[0,:])
 
# Create function to be called when slider value is changed
 
def update(val):
    i = t_slider.val
    line.set_ydata(rho[i,:])
    fig.canvas.draw_idle()  # Redraw the figure
 
# Call update function when slider value is changed
t_slider.on_changed(update)

plt.show()
