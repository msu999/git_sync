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
#   y' = (y_1', y_2') = (x', x'') = (y_2, -y_1) =: f(y)
#
# ---------------------------------------------------------

m = 10  # Mass [kg] of "paperclip"
M = 80  # Mass [kg] of "cup"
L = 50  # Total length of string/rope
g = 9.82  # Gravitational acceleration

# We reduce problem to an ODE by setting:
#
#   y_0 := r
#   y_1 := r'
#   y_2 := theta
#   y_3 := theta'
#

# Equations of motion A and B
def A(y_0, y_2, y_3):
    return (-M*g + m*g*np.sin(y_2) + m*y_0*(y_3**2))/(m + M)

def B(y_0, y_1, y_2, y_3):
    return (g*np.cos(y_2) - 2*y_1*y_3)/y_0

# y' = f(t, y) to iterate using RK4
def f(t, y):
    return np.array([y[1], A(y[0], y[2], y[3]), y[3], B(y[0], y[1], y[2], y[3])])

# Time interval to solve for y in:
t_0 = 0
t_1 = 4

# Initial condition:
y_0 = np.array([30, 0, 0, 0])

# Number of steps in RK4 iteration
N = 1000

# RK4 function returns sol = (t, y) where y = (y_1, y_2, ...). Hence to plot
# y_i we plot sol[0] against sol[1][i].
sol = rk4(f, y_0, t_0, t_1, N)

t = sol[0]  # Time vector
dt = t[1] - t[0]  # Time step in iteration

r = sol[1][0]  # r-vector
theta = sol[1][1]  # theta-vector

print(dt)

x = r*np.cos(theta)  # x-vector
y = r*np.sin(theta)  # y-vector
z = L - r  # z-vector

# Coordinates of points corresponding to "paperclip" and "cup" respectively
x_1 = x
y_1 = y

x_2 = np.zeros(z.size)  # "Cup" only falls straight down, i.e. x_2 = 0 for all t
y_2 = -z

fig_rt, ax_rt = plt.subplots()  # Prepares plot of r and theta

ax_rt.plot(t, r, label="r")
ax_rt.plot(t, theta, label="theta")
ax_rt.legend()

fig_xyz, ax_xyz = plt.subplots()  # Prepares plot of "cup" and "paperclip"

ax_xyz.plot(x_1, y_1, label="xy")
ax_xyz.plot(x_2, y_2, label="z")
ax_xyz.legend()

ax_xyz.set_aspect("equal")

time_template = 'time = %.1fs'
time_text = ax_xyz.text(0.05, 0.9, '', transform=ax_xyz.transAxes)

# Points corresponding to "paperclip" and "cup" respectively:
point_1, = ax_xyz.plot(x_1[0], y_1[0], 'bo', markersize=8)  # Initialize first point
point_2, = ax_xyz.plot(x_2[0], y_2[0], 'bo', markersize=8)  # Initialize second point

# Update function for animation
def update(frame):
    point_1.set_data([x_1[frame]], [y_1[frame]])
    point_2.set_data([x_2[frame]], [y_2[frame]])
    
    time_text.set_text(time_template % (frame*dt))

    return point_1, point_2, time_text

# Creating animation:
brake = 0.5  # Control animation speed

ani = animation.FuncAnimation(fig_xyz, update, frames=len(t), interval=brake*dt*1000, blit=True)

plt.show()
