import numpy as np
from pypoman import compute_polytope_vertices
from scipy.spatial import Delaunay
from pysimplicialcubature.simplicialcubature import (
    integrateOnSimplex,
    integratePolynomialOnSimplex,
)
from sympy import linear_eq_to_matrix, And, LessThan


def integrateOnPolytope(
    f, A, b, dim=1, maxEvals=10000, absError=0.0, tol=1.0e-5, rule=3
):
    """
    Integration a function over a convex polytope.
    
    Parameters
    ----------
    f : function
        The function to be integrated.
    A : array-like
        matrix of the coefficients of the linear inequalities (see README for an example)
    b: vector-like
        vector of the upper bounds of the linear inequalities
    dim : integer
        The dimension of the values of `f`.
    maxEvals : integer
        Maximum number of calls to `f`.
    absError : number
        Desired absolute error.
    tol : number
        Desired relative error.
    rule : integer 
        Integer between 1 and 4 corresponding to the integration rule; a 2*rule+1 degree rule will be applied.
        
    Returns
    -------
    dictionary
        The value of the integral is in the field `"integral"` of the returned value.
        
    Examples
    --------
    >>> from pypolyhedralcubature.polyhedralcubature import *
    >>> from sympy.abc import x, y, z
    >>> # integral bounds
    >>> i1 = (x >= -5) & (x <= 4)
    >>> i2 = (y >= -5) & (y <= 3 - x)
    >>> i3 = (z >= -10) & (z <= 6 - x - y)
    >>> # get matrix-vector representation of these inequalities
    >>> A, b = getAb([i1, i2, i3], [x, y, z])
    >>> # function to integrate
    >>> f = lambda x, y, z : x*(x+1) - y*z**2
    >>> # integral of f on the polytope defined by the bounds
    >>> g = lambda v : f(v[0], v[1], v[2])
    >>> I_f = integrateOnPolytope(g, A, b)
    >>> I_f["integral"]

    """
    vertices = compute_polytope_vertices(A, b)
    dlnay = Delaunay(vertices)
    tetrahedra = np.asarray(vertices)[dlnay.simplices]
    return integrateOnSimplex(
        f,
        tetrahedra,
        dim=dim,
        maxEvals=maxEvals,
        absError=absError,
        tol=tol,
        rule=rule,
        info=False,
    )


def integratePolynomialOnPolytope(P, A, b):
    """
    Integration a polynomial over a convex polytope.
    
    Parameters
    ----------
    P : function
        The function to be integrated.
    A : array-like
        matrix of the coefficients of the linear inequalities (see README for an example)
    b: vector-like
        vector of the upper bounds of the linear inequalities
        
    Returns
    -------
    number
        The exact value of the integral of P over the polytope defined by A and b.
        
    Examples
    --------
    >>> from pypolyhedralcubature.polyhedralcubature import *
    >>> from sympy import Poly
    >>> from sympy.abc import x, y, z
    >>> # integral bounds
    >>> i1 = (x >= -5) & (x <= 4)
    >>> i2 = (y >= -5) & (y <= 3 - x)
    >>> i3 = (z >= -10) & (z <= 6 - x - y)
    >>> # get matrix-vector representation of these inequalities
    >>> A, b = getAb([i1, i2, i3], [x, y, z])
    >>> # polynomial to integrate
    >>> P = Poly(x*(x+1) - y*z**2, domain = "RR")
    >>> # integral of P over the polytope defined by the bounds
    >>> integratePolynomialOnPolytope(P, A, b)

    """
    vertices = compute_polytope_vertices(A, b)
    dlnay = Delaunay(vertices)
    tetrahedra = np.asarray(vertices)[dlnay.simplices]
    integral= 0.0
    for tetrahedron in tetrahedra:
        integral = integral + integratePolynomialOnSimplex(P, tetrahedron)
    return integral


def __getAb0(inequalities, symbols, required_type):
    # assumptions:
    # 1. all inequalities are written with the same relational: < or <= or > or >=
    # 2. For each inequality, LHS and RHS are linear in `symbols`
    def get_ineq_with_correct_type(i, required_type):
        if type(i) != required_type:
            i = i.reversed
        return i

    # extract all inequalities, process them so they are all of the same type and
    # terms containing `symbols` are on the LHS.
    ineq = []
    for i in inequalities:
        if isinstance(i, And):
            ineq.extend([get_ineq_with_correct_type(a, required_type) for a in i.args])
        else:
            ineq.append(get_ineq_with_correct_type(i, required_type))
    # at this point, all inequalities should be of the same type.
    # rewrite them as expressions: LHS - RHS
    equations = [i.lhs - i.rhs for i in ineq]
    return linear_eq_to_matrix(equations, symbols)


def getAb(inequalities, symbols):
    A, b = __getAb0(inequalities, symbols, LessThan)
    return np.array(A, dtype="float"), np.array(b, dtype="float")[:, 0]
