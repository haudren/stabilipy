Backends
********

These are the suitable backends for static stability polygon 
computation. Some backends are restricted to 2D/3D cases.

They all take as argument, a geometry engine. For now, only scipy
is supported. Others are defined at least partially in geomengines.py:

* scipy: default and the only one supported as of now.
  We use its bindings of qhull.
* CGAL: Was supported but the available python bindings are too slow.
* Shapely: Does not support 3D properly
* Qhull-Sch: Custom bindings to qhullcpp for sch, that are not really
  usable as of now.

.. automodule:: stability.backends
   :members:
