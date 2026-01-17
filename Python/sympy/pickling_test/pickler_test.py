from sympy import *
from pickler import *

x = symbols("x")
f = Function("f")

expr1 = f(x) + x**2
expr2 = expr1.diff(x)

data = {"expr1" : expr1, "expr2" : expr2}

storeData(data, "db")
loadData("db", showKeys = True)
