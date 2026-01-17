from sympy import *

def ei(alpha):
    return ImmutableMatrix([cos(alpha), sin(alpha)])

def perp(v):
    if isinstance(v, ImmutableMatrix) and v.shape == (2, 1):
        return ImmutableMatrix([-v[1], v[0]])

    return None

def vector_integral(v, t):
    if isinstance(v, ImmutableMatrix) and v.shape == (2, 1):
        return ImmutableMatrix([Integral(v[0], t), Integral(v[1], t)])

    return None

def rotmat(theta):
    return ImmutableMatrix([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])

