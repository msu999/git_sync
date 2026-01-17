from sympy import *
from expressions import alpha, beta, t, R, r, omega, ga_beta, d_beta_ga, R_0, r_0, omega_0, Delta, delta, eta, s, la, u, theta, d_beta_Delta
import sys

sys.path.insert(0, "/home/max/Documents/Python/pickler")  # Adds pickler folder to system path

from pickler import storeData, loadData  # Imports pickler functions

xi = Function(r"xi")
data = {"dK_var" : xi}

storeData(data, "G_db")
loadData("G_db", showKeys=True)
