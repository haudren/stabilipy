Stability Polygon and contacts
******************************

.. autoclass:: stability.StabilityPolygon
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
