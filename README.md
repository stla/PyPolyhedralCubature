# pypolyhedralcubature

<!-- badges: start -->
[![Documentation status](https://readthedocs.org/projects/pypolyhedralcubature/badge/)](http://pypolyhedralcubature.readthedocs.io)
<!-- badges: end -->

*Multiple integration over a convex polytope.*

___

This package allows to evaluate a multiple integral whose integration 
bounds are some linear combinations of the variables, e.g.

$$\int\_{-5}^4\int\_{-5}^{3-x}\int\_{-10}^{6-x-y} f(x, y, z) \text{d}z \text{d}y \text{d}x.$$

In other words, the domain of integration is given by a set of linear 
inequalities:

$$\left\{\begin{matrix} -5  & \leq & x & \leq & 4 \\\ -5  & \leq & y & \leq & 3-x \\\ -10 & \leq & z & \leq & 6-x-y \end{matrix}\right..$$

These linear inequalities define a convex polytope (in dimension 3, a 
polyhedron). 
In order to use the package, one has to get the *matrix-vector* representation 
of these inequalities, of the form

$$A {(x,y,z)}' \leqslant b.$$

The matrix $A$ and the vector $b$ appear when one rewrites the linear 
inequalities above as:

$$\left\{\begin{matrix} -x & \leq & 5 \\\ x & \leq & 4 \\\ -y & \leq & 5 \\\ x+y & \leq & 3 \\\ -z & \leq & 10 \\\ x+y+z & \leq & 6 \end{matrix}\right..$$

The matrix $A$ is given by the coefficients of $x$, $y$, $z$ at the 
left-hand sides, and the vector $b$ is made of the upper bounds at the 
right-hand sides:

```python
import numpy as np
A = np.array([
  [-1, 0, 0], # -x
  [ 1, 0, 0], # x
  [ 0,-1, 0], # -y
  [ 1, 1, 0], # x+y
  [ 0, 0,-1], # -z
  [ 1, 1, 1]  # x+y+z
])
b = np.array([5, 4, 5, 3, 10, 6])
```

The function `getAb` provided by the package allows to get $A$ and $b$ in a 
user-friendly way:

```python
from pypolyhedralcubature.polyhedralcubature import getAb
from sympy.abc import x, y, z
# linear inequalities defining the integral bounds
i1 = (x >= -5) & (x <= 4)
i2 = (y >= -5) & (y <= 3 - x)
i3 = (z >= -10) & (z <= 6 - x - y)
# get the matrix-vector representation of these inequalities
A, b = getAb([i1, i2, i3], [x, y, z])
```

Now assume for example that $f(x,y,z) = x(x+1) - yz^2$. Once we have $A$ and 
$b$, here is how to evaluate the integral of $f$ over the convex polytope:

```python
from pypolyhedralcubature.polyhedralcubature import integrateOnPolytope
# function to integrate
f = lambda x, y, z : x*(x+1) - y*z**2
# integral of f over the polytope defined by the linear inequalities
g = lambda v : f(v[0], v[1], v[2])
I_f = integrateOnPolytope(g, A, b)
I_f["integral"]
# 57892.275000000016
```

In the case when the function to be integrated is, as in our current example, 
a polynomial function, it is better to use the `integratePolynomialOnPolytope` 
function provided by the package. This function implements a procedure 
calculating the exact value of the integral. Here is how to use it:

```python
from pypolyhedralcubature.polyhedralcubature import integratePolynomialOnPolytope
from sympy import Poly
from sympy.abc import x, y, z
# polynomial to integrate
P = Poly(x*(x+1) - y*z**2, domain = "RR")
# integral of P over the polytope 
integratePolynomialOnPolytope(P, A, b)
# 57892.2750000001
```

## Acknowledgments

I am grateful to the StackOverflow user @Davide_sd for the help he provided 
regarding the `getAb` function.
  