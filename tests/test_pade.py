import sympy
from sympy.functions import exp
from sympy import S, poly

def exp_taylor_coeff(order):
    x = sympy.symbols('x')
    _series = exp(x).series(n=order+1)
    return tuple(_series.coeff(x, p) for p in range(order+1))

def create_pade_system(coeffs, p, q):
    """returns a tuple of co-efficients (p, q) in terms of q0"""
    sz = len(coeffs)
    _p = sympy.Matrix(p + (S(0),) * (sz - len(p)))
    _q = sympy.Matrix(q + (S(0),) * (sz - len(q)))
    A = sympy.Matrix(sz, sz, lambda i, j: coeffs[:i+1][-(j+1)] if j <= i else S(0))
    eqns = list(A * _q - _p)
    (sol,) = sympy.linsolve(eqns, *p, *q[1:])
    return (sol[:len(p)], q[0:1] + sol[len(p):])

def main():
    L, M = 12, 12
    coeffs = exp_taylor_coeff(L + M)
    p = sympy.symbols(' '.join('p%d' % n for n in range(0, L+1)))
    q = sympy.symbols(' '.join('q%d' % n for n in range(0, M+1)))
    _coeffP, _coeffQ = create_pade_system(coeffs, p, q)

    # convention is to scale results by q[0], multiply out fractions
    xf = lambda x: x.subs({q[0]: S(1295295050649600)})
    coeffP, coeffQ = tuple(map(xf, _coeffP)), tuple(map(xf, _coeffQ))
        
    print("%r\n%r" % (coeffP, coeffQ))
    
if __name__ == '__main__':
    main()
    
