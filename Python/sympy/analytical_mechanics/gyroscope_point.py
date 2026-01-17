import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.animation as animation
import numpy as np

phi_0 = 0
theta_0 = 0
psi_0 = 0

R = 2
r = 1
origin = [0, 0, 0]

def R_point(phi, theta, psi):
    return [R*np.sin(theta)*np.cos(phi), R*np.sin(theta)*np.sin(phi), R*np.cos(theta)]

def r_point(phi, theta, psi):
    return [r*(-np.cos(psi)*np.sin(phi) - np.sin(psi)*np.cos(theta)*np.cos(phi)), r*(np.cos(psi)*np.cos(phi) + np.sin(psi)*np.cos(theta)*np.sin(phi)), r*np.sin(psi)*np.sin(theta)]

fig = plt.figure()
ax = fig.add_subplot(projection="3d")

R_vec = ax.quiver(*origin, *R_point(phi_0, theta_0, psi_0), color='r')
r_vec = ax.quiver(*origin, *r_point(phi_0, theta_0, psi_0), color='r')

plt.subplots_adjust(bottom=0.25)  # Leave space for the slider

# Define slider positions
ax_phi = plt.axes([0.2, 0.1, 0.65, 0.03])
ax_theta = plt.axes([0.2, 0.15, 0.65, 0.03])
ax_psi = plt.axes([0.2, 0.2, 0.65, 0.03])

# Create sliders
slider_phi = Slider(ax_phi, 'Phi (°)', 0, 360, valinit=phi_0)
slider_theta = Slider(ax_theta, 'Theta (°)', 0, 360, valinit=theta_0)
slider_psi = Slider(ax_psi, 'Psi (°)', 0, 360, valinit=psi_0)

# Update function for animation
def update(val):
    phi = slider_phi.val
    theta = slider_theta.val
    psi = slider_psi.val
    
    global R_vec, r_vec
    R_vec.remove()
    r_vec.remove()

    R_vec = ax.quiver(*origin, *R_point(phi, theta, psi), color='r')
    r_vec = ax.quiver(*origin, *r_point(phi, theta, psi), color='r')

    fig.canvas.draw_idle()

# Connect sliders to update function
slider_phi.on_changed(update)
slider_theta.on_changed(update)
slider_psi.on_changed(update)

plt.show()
