Recursive Projection
********************

Description of the algorithm
============================
Recursive projection is an algorithm designed to compute the projection of a convex set. In general, it computes the approximation of a smooth set :math:`P`.
To do so, it generates converging polyhedral approximations of :math:`P`, :math:`P_{inner}` and :math:`P_{outer}`.
The algorithm only needs an "oracle", i.e. an algorithm or optimization problem that yields the extremal point of :math:`P` in a given direction :math:`d`.
As a generic optimization problem:

.. math::
  max. \qquad &d^T x\\
  s.t. \qquad &x \in P

The solution :math:`x^*` of this problem is an extremal point in the direction :math:`d`.
By solving repeatedly the above problem, we obtain:

* The convex hull of all :math:`x^*` is :math:`P_{inner}`.
* The intersection of all halfspaces :math:`\{x \in \mathbb{R}^n | d^T x \leq d^T x^*\}` forms :math:`P_{outer}`

Now, the important point is how to choose the approriate sequence of directions :math:`d`.
To do so, we compute the *uncertainty volumes*, i.e. the cuts of :math:`P_{outer}` by the supporting hyperplanes of :math:`P_{inner}`.
The direction :math:`d` is chosen to be perpendicular to the supporting hyperplane that forms the largest uncertainty volume.

This is very useful to compute projections.
Consider a convex body :math:`P` in :math:`\mathbb{R}^{n+m}`.
Computing the projection of this body onto :math:`\mathbb{R}^n` can be done by specifying the following optimization problem:

.. math::
  max \qquad &d^T x\\
  s.t. \qquad &(x, y) \in P

This is particularly interesting when :math:`m >> n`. In this case, computing the direct projection (see for example `this page <https://scaron.info/teaching/projecting-polytopes.html>`_) is prohibitively expensive as the complexity is exponential in :math:`m+n`. In our case, the complexity depends on the class of optimization problem being solved, but is typically polynomial in the dimension.

For more information, please refer to `this paper <https://hal-lirmm.ccsd.cnrs.fr/lirmm-01477362>`_.

How to use this class
=====================

This class is intended for developpers and researchers who wish to implement a new class of problems.
If you are looking to compute stability or robust stability polygons and polyhedrons, please use :class:`stability.StabilityPolygon`.
If you wish to compute the projection of a set of linear inequalities, please use :class:`stability.LinearProjection`.
In general, one only needs to override the :py:func:`stability.RecursiveProjectionProblem.solve` method.

Let us have a look at an example (available in `sphere.py`)::

   import stability as stab


   class SphereProjection(stab.RecursiveProjectionProblem):
     """Try to approximate a sphere of radius r"""

     def __init__(self, radius):
       """:param radius: Radius of the sphere
          :type radius: double"""
       stab.RecursiveProjectionProblem.__init__(self, dimension=3)
       self.radius = radius

     def solve(self, d):
       """We are computing a sphere so the extremal point in direction d is just r*d"""
       return self.radius*d

   if __name__ == '__main__':
     sphere = SphereProjection(1.0)
     sphere.compute(solver='cdd', mode=stab.Mode.iteration, maxIter=50)

In this example we:

* Implement a class that extends RecursiveProjectionProblem
* Override solve to return the extremal point in the provided direction
* Instanciate that object and compute the approximation using `cdd` as our double-description package
* By default, this call will print the precision at each iteration and plot the result

See below for details of the API.


Class API
=========

.. autoclass:: stability.RecursiveProjectionProblem
  :members:

.. class:: stability.Mode

  All polygon computations should select a mode
  of operation.

  * precision: will try to reach desired precision,
    however many iterations are required.
  * iteration: will iterate for a number of iterations
    and stop, regardless of accuracy.
  * best: will try to reach desired precision under the
    given number of iterations.
