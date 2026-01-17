import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

omega = 1

A = 1
theta = np.pi/4 
phi_x = np.pi/2
phi_y = 0
diff_phi = phi_x - phi_y

alpha = np.sin(2*theta)*np.sin(diff_phi)

a = A/np.sqrt(2) * np.sqrt(1 + np.sqrt(1 - alpha**2))
b = -A/np.sqrt(2) * alpha/np.sqrt(1 + np.sqrt(1 - alpha**2))

phi = 1/2 * np.arcsin(A**2/(a**2 - b**2) * np.sin(2*theta)*np.cos(diff_phi))

Psi_nc = np.array([A*np.cos(theta)*np.exp(1j*phi_x), A*np.sin(theta)*np.exp(1j*phi_y)])
Psi_c = np.array([a*np.cos(phi) - 1j*b*np.sin(phi), a*np.sin(phi) + 1j*b*np.cos(phi)])

N = 100
t = np.linspace(0, 10, N)
dt = t[1] - t[0]
propagator = np.exp(-1j*omega*t)

x_nc = np.real(Psi_nc[0]*propagator)
y_nc = np.real(Psi_nc[1]*propagator)

x_c = np.real(Psi_c[0]*propagator)
y_c = np.real(Psi_c[1]*propagator)


# --- Plotting ---
def get_arrow(frame):
    x = np.cos(theta)
    y = np.sin(theta)
    z = 0
    u = np.sin(2*theta)
    v = np.sin(3*theta)
    w = np.cos(3*theta)
    return x,y,z,u,v,w


fig, ax = plt.subplots()
ax.set(title="Canonical form")
ax.plot(x_c, y_c, color="blue", label="Canonical form")


fig_, ax_ = plt.subplots()
ax_.set(title="Non-canonical form")
ax_.plot(x_nc, y_nc, color="purple", label="Non-canonical form")


ani_x_vec = ax.quiver(0, 0, x_c[0], y_c[0], scale_units="xy", scale=1, color="black")  # Initial vector direction

time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

# Update function for animation
def update(frame):
    i = frame

    global ani_x_vec
    ani_x_vec.remove()
    ani_x_vec = ax.quiver(0, 0, x_c[i], y_c[i], scale_units="xy", scale=1, color="black")
    
    time_text.set_text(time_template % (i*dt))

    return time_text

# Creating animation:
ani = animation.FuncAnimation(fig, update, frames=t.size, interval=50)
#ani.save('AM_ODE_friction.mp4', fps=fps)

plt.show()

