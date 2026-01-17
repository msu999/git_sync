from importlib.machinery import SourceFileLoader
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

rk4 = SourceFileLoader("rk4", "../rk4.py").load_module().rk4


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


# --- Problem parameters ---
# -- Main parameters --
m = 10  # Mass [kg] of "paperclip"
M = 85  # Mass [kg] of "cup"
L = 15  # Total length of string/rope [m]
g = 9.82  # Gravitational acceleration [m/s^2]
r_0 = 12  # Initial value for r [m]
theta_0 = 0  # Initial angle of "paperclip" [rad]


# -- Experimental parameters
l = 0  # Friction proportionality constant (0 = no friction)
R = 0  # Radius of rod which string/rope coils around (0 = negligible radius) [m]


# We reduce problem to an ODE by setting:
#
#   y_0 := r
#   y_1 := r'
#   y_2 := theta
#   y_3 := theta'
#


# --- Problem setup ---
# Equations of motion A and B
def A(y_0, y_2, y_3):
    return (-M*g + m*g*np.sin(y_2) - m*y_0*(y_3**2))/(m + M)

def B(y_0, y_1, y_2, y_3):
    return (-g*np.cos(y_2) - 2*y_1*y_3 - l*y_0*y_2)/y_0 - (M*g*R + l*(R**2)*y_2)/(m*(y_0**2))

# y' = f(t, y) to iterate using RK4
def f(t, y):
    return np.array([y[1], A(y[0], y[2], y[3]), y[3], B(y[0], y[1], y[2], y[3])])

# Time interval to solve for y in:
t_0 = 0
t_1 = 4

# Initial condition:
y_0 = np.array([r_0, 0, 0, 0])


# --- Animation parameters (+ resolution of solution) ---
fps = 30
tot_frames = fps*(t_1 - t_0)

# Resolution of solution (how many iterations to run per frame, a value of 100 yields nice looking solutions)
iterations_per_frame = 100


# --- Solving using RK4 ---
# Number of steps in RK4 iteration
N = iterations_per_frame*tot_frames

#print(f"fps = {fps}\ntotal number of frames = {tot_frames}\niterations per frame = {iterations_per_frame}\nN = {N}")  # For debugging

# RK4 function returns sol = (t, y) where y = (y_1, y_2, ...). Hence to plot
# y_i we plot sol[0] against sol[1][i].
sol = rk4(f, y_0, t_0, t_1, N)


# --- Solution vectors ---
# -- Main vectors --
t = sol[0]  # Time vector
dt = t[1] - t[0]  # Time step in iteration

#print(f"dt = {dt}")  # For debugging

r = sol[1][0]  # r-vector
theta = sol[1][2]  # theta-vector

x = r*np.cos(theta)
y = r*np.sin(theta)
z = -L + r


# -- Sliced vectors (to remove chaotic behaviour in limit) --
cutoff_r = 0.1  # !! Bound for which values for r which should be plotted. Must be adjusted for different masses m and M (for m = 10 and M = 85 a value r = 0.1 was used for example) !!

# If we model the rod around which the string/rope coils then we only want to plot those value for which r > R.
if R > 0:
    cutoff_r = R

cutoff_index = 0  # Index corresponding to last r to be plotted

# Finds the correct value for cutoff_index
while r[cutoff_index] > cutoff_r:
    cutoff_index += 1

# Vector which are cut off when chaotic behaviour begins
t_sliced = t[:cutoff_index + 1]

r_sliced = r[:cutoff_index + 1]
theta_sliced = theta[:cutoff_index + 1]

x_sliced = x[:cutoff_index + 1] 
y_sliced = y[:cutoff_index + 1]
z_sliced = z[:cutoff_index + 1]


# --- Plotting solution ---
# -- r-theta-plot --
fig_rt, ax_rt = plt.subplots()  # Prepares plot of r and theta
ax_rt.set(title=rf"$r\theta$-plot till $t \approx {t[cutoff_index]} \: s$", xlabel=r"$\theta$ [$\text{rad}$]", ylabel=r"$r$ [$m$]")

ax_rt.plot(theta_sliced, r_sliced)


# -- xyz-plot --
fig_xyz, ax_xyz = plt.subplots()  # Prepares plot of "cup" and "paperclip"
ax_xyz.set(title=rf"Rörelseplot till $t \approx {t[cutoff_index]} \: s$", xlabel=r"[$m$]", ylabel=r"[$m$]")
ax_xyz.set_xlim(-1.1*L, 1.1*L)
ax_xyz.set_ylim(-1.1*L, 1.1*L)
ax_xyz.set_aspect("equal")

ax_xyz.plot(x_sliced, y_sliced, label="\"Gemets\" rörelse")
ax_xyz.plot(np.zeros(len(z_sliced)), z_sliced, label="\"Koppens\" rörelse")


# -- Animation in xyz-plot --
# Coordinates of points corresponding to "paperclip" and "cup" respectively
x_1 = x
y_1 = y

x_2 = np.zeros(z.size)  # "Cup" only falls straight down, i.e. x_2 = 0 for all t
y_2 = z

time_template = 'time = %.1fs'
time_text = ax_xyz.text(0.05, 0.9, '', transform=ax_xyz.transAxes)

# Points corresponding to "paperclip" and "cup" respectively:
#point_1, = ax_xyz.plot(x_1[0], y_1[0], 'bo', markersize=8)  # Initialize first point
#point_2, = ax_xyz.plot(x_2[0], y_2[0], 'bo', markersize=8)  # Initialize second point

# Update function for animation
def update(frame):
    i = frame*iterations_per_frame

    # print(i)  # For debugging

    point_1.set_data([x_1[i]], [y_1[i]])
    point_2.set_data([x_2[i]], [y_2[i]])
    
    time_text.set_text(time_template % (i*dt))

    return point_1, point_2, time_text

# Creating animation:
brake = 1  # Control animation speed (1 = real-time in saved video)

#ani = animation.FuncAnimation(fig_xyz, update, frames=tot_frames, interval=brake*dt*1000, blit=True)
#ani.save('AM_ODE_friction.mp4', fps=fps)


plt.show()
