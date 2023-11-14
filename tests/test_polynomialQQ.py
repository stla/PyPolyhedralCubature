from pypolyhedralcubature.polyhedralcubature import(
        integratePolynomialOnPolytope, getAb
    )
from sympy import Poly
from sympy.abc import x, y, z
from fractions import Fraction

def test_polynomialQQ():
    # integral bounds
    i1 = (x >= -5) & (x <= 4)
    i2 = (y >= -5) & (y <= 3 - x)
    i3 = (z >= -10) & (z <= 6 - x - y)
    # get matrix-vector representation of these inequalities
    A, b = getAb([i1, i2, i3], [x, y, z])
    # polynomial to integrate
    P = Poly(x*(x+1) - y*z**2, domain = "QQ")
    # integral of P over the polytope
    I_P = integratePolynomialOnPolytope(P, A, b)
    # compare results
    assert I_P == Fraction(2315691, 40)