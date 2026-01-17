import sympy as sp
import numpy as np

L = 1
g = 9.82
k = 1
m = 1

t, A, B, C, D = sp.symbols("t, A, B, C, D", real=True)

w_0 = sp.sqrt(g/L)
w_plus = w_0
w_minus = w_0 + 2*k/m

psi_plus = A*sp.cos(w_plus*t) + B*sp.sin(w_plus*t)
psi_minus = C*sp.cos(w_minus*t) + D*sp.sin(w_minus*t)
Psi = sp.Matrix([psi_plus, psi_minus])

A = 1/2 * np.array([[1, 1], [1, -1]])

Phi = A*Psi
Phi_t = Phi.diff(t)

# Initial conditions
phi_0, theta_0 = sp.symbols(r"\phi_0, \theta_0")
omega_plus, omega_minus = sp.symbols(r"\omega_+, \omega_-")
alpha_A = sp.symbols(r"\alpha_A")  #2*sp.pi/10
alpha_B = sp.symbols(r"\alpha_B")  #2*sp.pi/10
alpha_C = sp.symbols(r"\alpha_C")  #2*sp.pi/10

v = sp.sqrt(2) * sp.ImmutableMatrix([phi_0, theta_0, 0, 0])
v_A = sp.sqrt(2)*sp.Matrix([alpha_A, alpha_A, 0, 0])
v_B = sp.sqrt(2)*sp.Matrix([-alpha_B, alpha_B, 0, 0])
v_C = sp.sqrt(2)*sp.Matrix([-alpha_C, 0, 0, 0])
M = sp.ImmutableMatrix([[1, 0, 1, 0], [1, 0, -1, 0], [0, omega_plus, 0, omega_minus], [0, omega_plus, 0, -omega_minus]])

print(M.det())

res = M.solve(v)
res_A = M.subs({phi_0 : alpha_A, theta_0 : alpha_A})
res_B = M.subs({phi_0 : -alpha_B, theta_0 : alpha_B})
res_C = M.subs({phi_0 : -alpha_C, theta_0 : 0})

print(f"Allmänt : [A, B, C, D] = {sp.latex(res)}")
print(f"Fall A : [A, B, C, D] = {sp.latex(res_A)}")
print(f"Fall B : [A, B, C, D] = {sp.latex(res_B)}")
print(f"Fall C : [A, B, C, D] = {sp.latex(res_C)}")

# Measured values 
w_plus_exp = np.array([4.53, 4.53, 4.53])  # Measure in A and plug in +- 0.00024, 0.00026, 0.00032 rad/s
w_minus_exp = np.array([4.81, 4.81, 4.81])  # Measure in B and plug in +-0.00014, 0.00014, 0.00021 rad/s

m_pendulum_weight = 70.41e-3 # Mass of weight on pendulum
m_pendulum_rod = 32.15e-3 # Mass of pendulum rod
a = 0.026 # Distance from top of rod to rotary sensor [m]
L = 0.524 # Distance from rotary sensor to cm of weight [m]
b = 0.1  # Distance from rotary sensor to spring [m]

# For calculating spring constant
m_spring = 15.67e-3  # Spring mass
m_weight_hanger = 19.52e-3  # Weight  mass
m_weight = 19.33e-3  # Weight

m_weights = np.array([m_weight_hanger, m_weight_hanger + m_weight, m_weight_hanger + 2*m_weight])

A = (60.5-44.5)*1e-2  # Spring elongation under its own weight
Bs = np.array([(67-44.5)*1e-2, (73.7-44.5)*1e-2, (80.5-44.5)*1e-2])  # Spring elongation under ()

k = m_weights*g/(Bs - A)
print(f"Spring constant = {np.mean(k)}")

K = (w_minus_exp**2 - w_plus_exp**2)/(w_minus_exp**2 + w_plus_exp**2)
