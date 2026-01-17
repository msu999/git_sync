from sympy.core.function import UndefinedFunction 
from sympy import *


# Checks whether sympy expression contains a vector (object with attribute vector=True)
def contains_vector(expr):
    ans = False

    expr = sympify(expr)  # Sympifies expression such that it can be analyzed using the "atoms" method

    # Checks whether any atom of Function or Symbol class is a vector
    for s in expr.atoms(Function).union(expr.atoms(Symbol)):
        if hasattr(s, "vector") and s.vector:
            ans = True

    return ans


#class VectorFunction(Function):
#    def __new__(cls, name, **kwargs):
#        obj = super().__new__(cls, name, **kwargs)
#        return obj
#   
#    def _eval_Mul(self, other):
#        return Mul()


class Perp(Expr):
    def __new__(cls, expr):
        return super().__new__(cls, expr)

    def __init__(self, expr):
        expr = sympify(expr)
        self.expr = expr

    def _eval_doit(self, **hints):
        return Perp(self.expr.doit(**hints))

    def _eval_subs(self, old, new):
        new_expr = self.expr.subs(old, new)

        if isinstance(new_expr, ImmutableMatrix) and new_expr.shape == (2, 1):
            return Matrix([-new_expr[1], new_expr[0]])  # Compute perpendicular vector
        
    def _eval_simplify(self, **kwargs):
        expr = self.expr

        if isinstance(expr, ImmutableMatrix) and expr.shape == (2, 1):
            return Matrix([-expr[1], expr[0]])  # Compute perpendicular vector

        # Return unevaluated function for symbolic expr
        return Perp(expr)

    def _latex(self, printer):
        return r"(%s)^{\perp}" % printer._print(self.args[0])

class Norm(Expr):
    def __new__(cls, expr):
        return super().__new__(cls, expr)

    def __init__(self, expr):
        expr = sympify(expr)
        self.expr = expr

    def _eval_doit(self, **hints):
        return Norm(self.expr.doit(**hints))

    def _eval_subs(self, old, new):
        expr = self.expr
        new_expr = expr.subs(old, new)

        if isinstance(new_expr, ImmutableMatrix) and new_expr.shape == (2, 1):
            return sqrt(new_expr[0]**2 + new_expr[1]**2)

        # Return symbolic expression if expr is not a 2D vector
        return Norm(new_expr)

    def _latex(self, printer):
        return r"\left|%s\right|" % printer._print(self.args[0])

class Int(Expr):
    _empty = symbols(r"\text{}")

    def __new__(cls, expr, x, a=_empty, b=_empty):
        expr = sympify(expr)
        x = sympify(x)
        a = sympify(a)
        b = sympify(b) 
        obj = Expr.__new__(cls, expr, x, a, b)
        obj.expr = sympify(expr)
        obj.x = x
        obj.a = a
        obj.b = b

        return obj

    def _doit(self, **hints):
        return Int(self.expr.doit(**hints), self.x, self.a, self.b)

    def _eval_simplify(self, **kwargs):
        expr, x, a, b = self.expr, self.x, self.a, self.b
        new_expr, new_x, new_a, new_b = expr.simplify(**kwargs), x.simplify(**kwargs), a.simplify(**kwargs), b.simplify(**kwargs)

        return Int(new_expr, new_x, new_a, new_b)

    def _eval_subs(self, old, new):
        return Int(self.expr.subs(old, new), self.x, self.a, self.b)
     
    def diff(self, x):
        return Int(self.expr.diff(x), self.x, self.a, self.b)

    def _latex(self, printer=None):
        return rf"\int_{{{printer._print(self.a)}}}^{{{printer._print(self.b)}}} {printer._print(self.expr)} \, d{printer._print(self.x)}"

class Bullet(Expr):
    def __new__(cls, a, b):
        a = sympify(a)
        b = sympify(b)
        obj = Expr.__new__(cls, a, b)
        obj.a = a
        obj.b = b

        return obj
   
    def _eval_doit(self, **hints):
        return Bullet(self.a.doit(**hints), self.b.doit(**hints))

    def _eval_Mul(self, other): 
        print("Multiplication with Bullet!")
        # Case: a * Bullet(b, c) -> Bullet(a*b, c)
        if not contains_vector(other):
            return Bullet((self.a * other), self.b).simplify()

        return None  # Default multiplication

    def _eval_subs(self, old, new):
        a, b = self.a, self.b
        new_a = a.subs(old, new)
        new_b = b.subs(old, new)
        
        if isinstance(new_a, MatrixExpr) and isinstance(new_b, MatrixExpr):
            if new_a.shape == (2, 1) and new_b.shape == (2, 1):
                return (new_a.T * new_b)[0].simplify(**kwargs)
        
        return Bullet(new_a, new_b)

    def _eval_simplify(self, **kwargs):
        a, b = self.a, self.b

        print(f"Running simplification for {self}")

        # Distribute over addition: Bullet(a + b, c) → Bullet(a, c) + Bullet(b, c)
        if isinstance(a, Add):  
            print("Dealing with add")
            return Add(*[Bullet(arg, b).simplify(**kwargs) for arg in a.args])
        
        if isinstance(b, Add):  
            print("Dealing with add")
            return Add(*[Bullet(a, arg).simplify(**kwargs) for arg in b.args])

        if isinstance(a, MatrixExpr) and isinstance(b, MatrixExpr):
            if a.shape == (2, 1) and b.shape == (2, 1):
                return (a.T * b)[0].simplify(**kwargs)
       
        # Allow Bullet object (scalar product) to commute with Int object (integral)
        if isinstance(a, Int) and not isinstance(b, Int):
            print(f"a = {a} is an integral and as such we will simplify!")
            return Int(Bullet(a.expr, b).simplify(**kwargs), a.x, a.a, a.b)

        if not isinstance(a, Int) and isinstance(b, Int):
            return Int(Bullet(a, b.expr).simplify(**kwargs), b.x, b.a, b.b)

        return Bullet(a, b)


    def _latex(self, printer=None):
        return rf"{{{printer._print(self.a)}}} \bullet {{{printer._print(self.b)}}}"
        #return r"%s \bullet %s" % (printer._print(self.a), printer._print(self.b))

#    def __mul__(self, other):
#        """Override * operator to use _eval_Mul."""
#        result = self._eval_Mul(other)
#        return result if result is not None else Mul(self, other)


# --- Custom functions ---
def ei(alpha):
    return Matrix([cos(alpha), sin(alpha)])

# --- Symbols ---
x, t, dOmega, R_0, alpha, beta, eps, w = symbols("x, t, \\partial\\Omega, R_0, \\alpha, \\beta, \\varepsilon, w", real=True)

DR = Function(r"\Delta R", real=True)
R = Function("R", real=True)
r = Function("r", real=True)

ga = Function("\\gamma", real=True, vector=True)
la = Function("\\lambda", real=True, vector=True)
#ga_exp = R*ei(alpha)
#la_exp = r*ei(u*t)
u = -1/(2*pi) * Int(log(Norm(x - ga(beta, t))) * ga(beta, t).diff(beta), beta, 0, 2*pi) 

#print(latex(diff(u_integrand, t)))
# 
d_tR = 1/(2*pi*R(alpha, t)) * Bullet(Int(log(Norm(ga(alpha, t) - ga(beta, t))) * ga(beta, t).diff(beta), beta, 0, 2*pi) - eps*Perp(ga(alpha, t) - la(t))/(Norm(ga(alpha, t) - la(t))**2), Perp(ga(alpha, t).diff(alpha)))

test = ga(alpha)*Bullet(2, 2)

#print(test.simplify())

#print(latex(d_tR))
print(latex(d_tR.doit().subs({ga(alpha, t) : (R_0 + DR(alpha, t))*ei(alpha), ga(beta, t) : (R_0 + DR(beta, t))*ei(beta), la(t) : r(t)*ei(w*t)}).simplify()))
