from sympy import *

rho_0, eps, R_E, r, theta, phi = symbols(r"\rho_0, \varepsilon, R_E, r, \theta, \phi")

R = R_E*(1 - (eps/2)*(3*cos(theta)**2 - 1))
rho = rho_0

x = r*sin(theta)*cos(phi)
y = r*sin(theta)*sin(phi)
z = r*cos(theta)

I_x = Integral(Integral(Integral(rho*(y**2 + z**2)*(r**2)*sin(theta), (r, 0, R)), (theta, 0, pi)), (phi, 0, 2*pi))
I_z = Integral(Integral(Integral(rho*(y**2 + x**2)*(r**2)*sin(theta), (r, 0, R)), (theta, 0, pi)), (phi, 0, 2*pi))

eps_small = {eps**2 : 0, eps**3 : 0, eps**4 : 0, eps**5 : 0}

print(latex(I_x.doit().expand().simplify().subs(eps_small).simplify()), end="\n\n")
print(latex(I_z))
