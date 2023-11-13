from pypolyhedralcubature.polyhedralcubature import integrateOnPolytope, getAb
from sympy.abc import x, y, z

def test_function():
    # integral bounds
    i1 = (x >= -5) & (x <= 4)
    i2 = (y >= -5) & (y <= 3 - x)
    i3 = (z >= -10) & (z <= 6 - x - y)
    # get matrix-vector representation of these inequalities
    A, b = getAb([i1, i2, i3], [x, y, z])
    # function to integrate
    f = lambda x, y, z : x*(x+1) - y*z**2
    # integral of f over the polytope defined by the bounds
    g = lambda v : f(v[0], v[1], v[2])
    I_f = integrateOnPolytope(g, A, b)
    # compare results
    assert abs(57892.275 - I_f["integral"]) < 1e-5