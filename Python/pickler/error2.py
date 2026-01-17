import pickle
from sympy.physics.mechanics import dynamicsymbols

q1 = dynamicsymbols('q1')

data = {"q1" : q1}

b = pickle.dumps(data)
