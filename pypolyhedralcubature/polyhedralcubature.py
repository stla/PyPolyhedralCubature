import numpy as np
from pypoman import compute_polytope_vertices
from scipy.spatial import Delaunay
from pysimplicialcubature.simplicialcubature import integrateOnSimplex, integratePolynomialOnSimplex
from sympy import linear_eq_to_matrix
from sympy.core.relational import And, LessThan

def integrateOnPolytope(f, A, b, dim=1, maxEvals=10000, absError=0.0, tol=1.0e-5, rule=3):
    vertices = compute_polytope_vertices(A, b)
    dlnay = Delaunay(vertices)
    tetrahedra = np.asarray(vertices)[dlnay.simplices]
    return integrateOnSimplex(f, tetrahedra, dim=dim, maxEvals=maxEvals, absError=absError, tol=tol, rule=rule, info=False)

def integratePolynomialOnPolytope(P, A, b):
    vertices = compute_polytope_vertices(A, b)
    dlnay = Delaunay(vertices)
    tetrahedra = np.asarray(vertices)[dlnay.simplices]
    return integratePolynomialOnSimplex(P, tetrahedra)

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