Computation of stability polygons
=================================

This package provides an easy interface to compute stability polygons.

One should create a StabilityPolygon by setting the robotMass, like so::

   import stabilipy
   poly = stabilipy.StabilityPolygon(57.5)

By default, a 3D robust static polyhedron is defined, see
:class:`stability.StabilityPolygon` for more information.

It is then necessary to create some contacts, and insert them in the polygon::

   pos = [[[0., 0., 1.]], [[1., 0., 0.]]]
   normals = [[[0., 0., 1.]], [[0.1, 0.1, 1.]]]
   mu = 0.7

   contacts = [stabilipy.Contact(mu, np.array(p).T,
                                 stabilipy.utils.normalize(np.array(n).T))
               for p, n in zip(pos, normals)]
   poly.contacts = contacts

Note that the normals *must* be of norm 1. You can now launch the computation::

   poly.compute(stability.Mode.best, epsilon=1e-3, maxIter=50)

Compute takes many arguments, see :func:`stability.StabilityPolygon.compute` but
the most importants are:

* Mode that determines if you want to reach a desired precision, iterate a
  number of times or best of both.
* epsilon sets the precision
* maxIter the number of iterations
* A number of plot_something keyword arguments are available.

Contents:

.. toctree::
   :maxdepth: 2

   recursive_projection
   stabpoly
   linear_projection
   backends
   constraints

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

