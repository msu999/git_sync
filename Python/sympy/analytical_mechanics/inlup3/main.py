from sympy import *

t, M, R, d, w_g, T = symbols(r"t, M, R, d, w_g, T", real=True)
alpha = symbols(r"\alpha", positive=True)
T_exp = sqrt(2*pi/alpha)

w_r = -abs(alpha*(t - T/2)) + alpha*T/2

w = Matrix([w_r*sin(w_g*t), w_g, w_r*cos(w_g*t)])
w_dot = w.diff(t)

w_x, w_y, w_z = w[0], w[1], w[2]
w_x_dot, w_y_dot, w_z_dot = w_dot[0], w_dot[1], w_dot[2]

I_x = M*(d**2/2 + R**2)
I_y = 2*M*R**2
I_z = I_x

tau_x = I_x*w_x_dot + (I_z - I_y)*w_y*w_z
tau_y = I_y*w_y_dot + (I_x - I_z)*w_x*w_z
tau_z = I_z*w_z_dot + (I_y - I_x)*w_x*w_y

print(latex(w_z_dot))
print(latex((w_x*w_y).expand().simplify()))

#n, d = tau_x.simplify().subs(abs(T - 2*t), -T + 2*t).args[1][0].as_numer_denom()

#print(latex(div(n, d)))
