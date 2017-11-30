Linear Projection
*****************

Principle
=========

Computing the projection of a convex set bounded by linear equalities and inequalities is a particular case of :ref:`recursive_projection`.
Indeed, in this case :math:`x \in P` can be directly written as:

.. math::

  A x &\leq b \\
  C x &= d

And thus, finding extremal points amounts to solving Linear Programs (LP).
Denoting the affine projection onto a smaller space by :math:`y = E x + f` (same convention as `Stéphane Caron <https://scaron.info/teaching/projecting-polytopes.html>`_), finding extremal points corresponding to a direction :math:`\delta` is done by solving:

.. math::

  max & \qquad \delta^T (E x + f)\\
  s.t. & \qquad A x \leq b\\
  & \qquad C x = d

A specific class is shown below.

Example of usage
================

The following example (contained in `hypercube.py`) shows how to project a 6D hypercube in 3D, resulting in a cube::

    import stability as stab
    import numpy as np

    if __name__ == '__main__':

       A = np.vstack((np.eye(6), -np.eye(6)))
       b = np.ones(12,)

       linear_proj = stab.LinearProjection(3, A, b, None, None)
       linear_proj.compute(stab.Mode.precision, solver='cdd', epsilon=1e-3)

Important notes:

* You need to specify the dimension you are projecting onto
* If you do not specify the projection operator :math:`E, f`, it will default to projecting on the first `dimension` dimensions.

Class API
=========

.. autoclass:: stability.LinearProjection
  :members:

