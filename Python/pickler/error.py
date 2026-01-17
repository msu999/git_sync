from sympy import Function
import pickle

A = Function("B")
data = {"A" : A}

b = pickle.dumps(data)
