from rk4 import rk4
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

# -------------------- INSTRUCTIONS -----------------------
#
# We have the following ODE problem:
#
#   y' = f(t, y)
#   y(t_0) = y_0
#
# For example:
#
#   x'' = -x
#   x(0) = 0
#   x'(0) = 1
#
#   (Solution x = sin(t))
#
# Gets rewritten into:
#
#   y_1 = x
#   y_2 = x'
#
# Which yields the ODE problem:
#
#   y' = (y_1', y_2') = (x', x'') = (y_2, -y_1) =: f(t, y)
#
# ---------------------------------------------------------

m = 4000  # Mass of satellite
g = 9.82  # Gravitational acceleration
k = 1e-9  # Proportionality constant for air resistance
h_0 = 1000*1000  # Height above Earth's surface which satellite starts in
R = 6378*1000  # Radius of earth
r_0 = h_0 + R  # Initial height of satellite
omega = np.sqrt(g/r_0)  # Angular velocity of satellite


# --- Animation parameters ---
fps = 30
video_len = 20  # 1 minutes real-time per frame, i.e. 30 minutes real-time per video second with fps = 30
tot_frames = fps*video_len


# --- Problem setup ---
# We reduce problem to an ODE by setting:
#
#   y_0 := r
#   y_1 := r'
#

def A(y_0):
    return (omega**2)*y_0 - g/(1 - k)

# y' = f(t, y) to iterate using RK4
def f(t, y):
    return np.array([y[1], A(y[0])])

# Time interval to solve for y in:
t_0 = 0
t_1 = video_len*5*60

# Initial condition:
y_0 = np.array([r_0, 0])

# Resolution of solution (how many iterations to run per frame)
iterations_per_frame = 10


# --- Solving using RK4 ---
# Number of steps in RK4 iteration
N = iterations_per_frame*tot_frames

#print(f"fps = {fps}\ntotal number of frames = {tot_frames}\niterations per frame = {iterations_per_frame}\nN = {N}")  # For debugging

# RK4 function returns sol = (t, y) where y = (y_1, y_2, ...). Hence to plot
# y_i we plot sol[0] against sol[1][i].
sol = rk4(f, y_0, t_0, t_1, N)


# --- Solution vectors ---
t = sol[0]  # Time vector
dt = t[1] - t[0]  # Time step in iteration

#print(f"dt = {dt}")  # For debugging

r = sol[1][0]  # r-vector
theta = omega*t  # theta-vector

x = r*np.cos(theta)  # x-vector
y = r*np.sin(theta)  # y-vector


# --- Plotting solution ---
# Coordinates of points corresponding to "paperclip" and "cup" respectively
x_1 = x
y_1 = y

fig_rt, ax_rt = plt.subplots()  # Prepares plot of r and theta

ax_rt.plot(t, r, label="r")
ax_rt.plot(t, theta, label="theta")
ax_rt.legend()

fig_xy, ax_xy = plt.subplots()  # Prepares plot of "cup" and "paperclip"

#ax_xy.plot(x_1, y_1, label="xy")
ax_xy.plot(x, y)
ax_xy.legend()

ax_xy.set_xlim([-1.1*r_0, 1.1*r_0])
ax_xy.set_ylim([-1.1*r_0, 1.1*r_0])
ax_xy.set_aspect("equal")

fig_ani, ax_ani = plt.subplots()  # Prepares animation plot
ax_ani.set_xlim([-1.1*r_0, 1.1*r_0])
ax_ani.set_ylim([-1.1*r_0, 1.1*r_0])
ax_ani.set_aspect("equal")

time_template = 'time = %.1fs'
time_text = ax_ani.text(0.05, 0.95, '', transform=ax_ani.transAxes)

r_template = "h = %.1fm above earth"
r_text = ax_ani.text(0.05, 0.01, "", transform=ax_ani.transAxes)

# Circle representing earth
earth = plt.Circle((0, 0), R, color='b')
ax_ani.add_patch(earth)

# Points corresponding to "paperclip" and "cup" respectively:
point, = ax_ani.plot(x[0], y[0], 'o', markersize=3, color="black")  # Initialize first point

# Update function for animation
def update(frame):
    i = frame*iterations_per_frame

    # print(i)  # For debugging

    point.set_data([x[i]], [y[i]])
    
    time_text.set_text(time_template % (i*dt))
    r_text.set_text(r_template % (r[i] - R))

    return point, time_text, r_text

# Creating animation:
brake = 1  # Control animation speed

ani = animation.FuncAnimation(fig_ani, update, frames=tot_frames, interval=50, blit=True)

ani.save('AM_ODE_satellite.mp4', fps=fps)
plt.show()
