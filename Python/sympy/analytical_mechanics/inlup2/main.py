#from sympy import *
import numpy as np
from scipy.constants import c

def gamma_factor(v):
    return 1/np.sqrt(1 - (v/c)**2)

def beta_factor(v):
    return v/c

def v_solve(l_1, l_2):
    return c*(1 - l_1/l_2)

def eigenlength(l, v, theta):
    # Consider a system moving with a velocity v relative to a stationary system. Let the moving system emit light. The wavelength l' of the light in the moving reference frame can be calculated as:
    #
    #   l' = l/(gamma * (1 + beta*cos(theta)))
    #
    # where gamma is the gamma factor, beta = v/c and theta is the angle that line connecting the origin of the stationary system with the moving object makes with the axis of movement when the light was emitted.

    gamma = gamma_factor(v)
    beta = beta_factor(v)

    return l/(gamma * (1 + beta*np.cos(theta)))

def wavelength(l_prime, v, theta):
    gamma = gamma_factor(v)
    beta = beta_factor(v)

    return l_prime * (gamma * (1 + beta*np.cos(theta)))

l_1 = 550e-9
l_2 = 720e-9

theta_1 = -np.pi
theta_2 = np.pi/2
theta_3 = 0

v = v_solve(l_1, l_2)

l_1_prime = eigenlength(l_1, v, theta_1)
l_2_prime = eigenlength(l_2, v, theta_2)

l_3 = wavelength(l_1_prime, v, theta_3)

print(f"v = {v/c} * c\nFirst eigenlength: {l_1_prime}\nSecond eigenlength: {l_2_prime}\nThird wavelength: {l_3}")
