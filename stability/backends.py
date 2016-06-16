#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy.spatial import ConvexHull
from scipy.spatial.qhull import QhullError
import numpy as np

import hashlib
import random

try:
  import cdd
  CDD_AVAILABLE = True
except ImportError:
  CDD_AVAILABLE = False

try:
  import pyparma
  from pyparma.utils import floatize, fractionize
  PPL_AVAILABLE = True
except ImportError:
  PPL_AVAILABLE = False

from plain_polygon import DoublePolygon

class SteppingException(Exception):
  def __init__(self, m):
    self.message = m

  def __str__(self):
    return self.message

class CDDBackend(object):

  """Using the CDD backend for polygon computation. This is the most polyvalent
  backend. Works on floating-point numbers. Requires pycddlib"""

  def __init__(self, geomengine='scipy'):
    """Default constructor.

    :param geomengine: underlying geometry engine. Only scipy is supported
    """

    if not CDD_AVAILABLE:
      raise ValueError("Cannot find cddlib. Make sure you install pycddlib.")
    self.name = 'cdd'

    if geomengine == 'scipy':
      self.volume_convex = self.scipy_volume_convex

  def scipy_volume_convex(self, hrep):
    try:
      points = np.array(cdd.Polyhedron(hrep).get_generators())[:, 1:]
    except RuntimeError:
      return 0

    if points.shape[0] < points.shape[1]+1:
      return 0
    try:
      ch = ConvexHull(points)
    except QhullError:
      return 0
    return ch.volume

  def build_polys(self, poly):
    if poly.outer is None:
      A = np.vstack(poly.directions)
      b = np.vstack(poly.offsets)
      poly.outer = cdd.Matrix(np.hstack([b, -A]))
      poly.outer.rep_type = cdd.RepType.INEQUALITY
    else:
      poly.outer.extend(np.hstack((poly.offsets[-1], -poly.directions[-1])))
      poly.outer.canonicalize()

    if poly.inner is None:
      A = np.vstack([p.T for p in poly.points])
      b = np.ones((len(poly.points), 1))
      poly.inner = cdd.Matrix(np.hstack((b, A)))
      poly.inner.rep_type = cdd.RepType.GENERATOR
    else:
      poly.inner.extend(np.hstack(([[1]], poly.points[-1].T)))
      poly.inner.canonicalize()

  def find_direction(self, poly, plot=False):
    self.build_polys(poly)

    volumes = []

    try:
      ineq = np.array(cdd.Polyhedron(poly.inner).get_inequalities())
    except RuntimeError:
      raise SteppingException('Numerical inconsistency found')

    for line in ineq:
      key = hashlib.sha1(line).hexdigest()
      if key in poly.volume_dic:
        volumes.append(poly.volume_dic[key])
      else:
        if key in poly.hrep_dic:
          A_e = poly.hrep_dic[key]
        else:
          A_e = poly.outer.copy()
          A_e.extend(cdd.Matrix(-line.reshape(1, line.size)))
          A_e.canonicalize()
          poly.hrep_dic[key] = A_e

        if plot:
          poly.reset_fig()
          poly.plot_polyhedrons()
          poly.plot_polyhedron(A_e, 'm', 0.5)
          poly.show()

        vol = self.volume_convex(A_e)
        poly.volume_dic[key] = vol
        volumes.append(vol)
        poly.vrep_dic[key] = np.array(cdd.Polyhedron(A_e).get_generators())

    maxv = max(volumes)
    alli = [i for i, v in enumerate(volumes) if v == maxv]
    i = random.choice(alli)
    key = hashlib.sha1(ineq[i, :]).hexdigest()
    return -ineq[i, 1:]

  def scipy_triangulate_polyhedron(self, hrep):
    points = np.array(cdd.Polyhedron(hrep).get_generators())[:, 1:]
    ch = ConvexHull(points)
    return [c for c in ch.points.T], ch.simplices

  def scipy_convexify_polyhedron(self, hrep):
    points = np.array(cdd.Polyhedron(hrep).get_generators())[:, 1:]
    ch = ConvexHull(points)
    return ch.points

class ParmaBackend(object):

  """Backend using the Parma Polyhedra Library
  This is the most precise, and thus slow backend. Works on integer (unlimited
  precision through the use of GMP). Requires pyparma."""

  def __init__(self, geomengine='scipy'):
    """Default constructor.

    :param geomengine: underlying geometry engine. Only scipy is supported
    """

    if not PPL_AVAILABLE:
      raise ValueError("Cannot find the ppl. Make sure you install pyparma.")
    self.name = 'parma'

    if geomengine == 'scipy':
      self.volume_convex = self.scipy_volume_convex

  def scipy_volume_convex(self, poly):
    vrep = poly.vrep()
    if vrep.size != 0:
      points = floatize(poly.vrep()[:, 1:])
    else:
      return 0

    try:
      ch = ConvexHull(points)
    except QhullError:
      return 0
    return ch.volume

  def build_polys(self, poly):
    if poly.outer is None:
      A = np.vstack(poly.directions)
      b = np.vstack(poly.offsets)
      poly.outer = pyparma.Polyhedron(hrep=fractionize(np.hstack([b, -A])))
    else:
      poly.outer.add_ineq(fractionize(np.hstack((poly.offsets[-1],
                                                 -poly.directions[-1]))))
    if poly.inner is None:
      A = np.vstack([p.T for p in poly.points])
      b = np.ones((len(poly.points), 1))
      poly.inner = pyparma.Polyhedron(vrep=np.hstack([b, fractionize(A)]))
    else:
      poly.inner.add_generator(np.hstack(([[1]],
                                          fractionize(poly.points[-1].T))))

  def find_direction(self, poly, plot=False):
    self.build_polys(poly)
    volumes = []

    ineq = poly.inner.hrep()

    for line in ineq:
      key = hashlib.sha1(line).hexdigest()
      if key in poly.volume_dic:
        volumes.append(poly.volume_dic[key])
      else:
        A_e = poly.outer.copy()
        A_e.add_ineq(-line)
        vol = self.volume_convex(A_e)
        poly.volume_dic[key] = vol
        volumes.append(vol)
        poly.vrep_dic[key] = A_e.vrep()

      if plot:
        poly.reset_fig()
        poly.plot_polyhedrons()
        poly.plot_polyhedron(A_e, 'm', 0.5)
        poly.show()

    maxi, maxv = max(enumerate(volumes), key=lambda x: x[1])
    alli = [i for i, v in enumerate(volumes) if v == maxv]
    i = random.choice(alli)
    return floatize(-ineq[i, 1:])

  def scipy_triangulate_polyhedron(self, poly):
    points = floatize(poly.vrep()[:, 1:])
    ch = ConvexHull(points)
    return [c for c in ch.points.T], ch.simplices

  def scipy_convexify_polyhedron(self, poly):
    points = floatize(poly.vrep()[:, 1:])
    ch = ConvexHull(points)
    return ch.points

class PlainBackend(object):

  """Plain Backend using the cdd backend for initialization.
  This is the simplest, fastest backend. However, only works on
  2D polygons."""

  def __init__(self, geomengine='scipy'):
    """Default constructor.

    :param geomengine: underlying geometry engine. Only scipy is supported.
    """

    if not CDD_AVAILABLE:
      raise ValueError("Cannot find cddlib. Make sure you install pycddlib.")
    self.name = 'parma'
    self.doublepoly = None
    self.edge_i = None

    if geomengine == 'scipy':
      self.volume_convex = self.scipy_volume_convex

  def scipy_volume_convex(self, poly):
    points = poly.vertices
    try:
      ch = ConvexHull(points)
    except QhullError:
      print("Unable to find convex hull")
      return 0
    return ch.volume

  def build_polys(self, poly):
    if poly.outer is None or poly.inner is None:
      poly.outer = None
      poly.inner = None
      self.doublepoly = None
      self.edge_i = None

    if poly.outer is None:
      A = np.vstack(poly.directions)
      b = np.vstack(poly.offsets)
      outer = cdd.Matrix(np.hstack([b, -A]))
      outer.rep_type = cdd.RepType.INEQUALITY
    else:
      i1, i2 = self.doublepoly.inner.edges[self.edge_i][1:]
      ei1, ei2 = self.doublepoly.correspondences[i1], self.doublepoly.correspondences[i2]
      self.doublepoly.addInequalityToOuter(poly.offsets[-1],
                                           -poly.directions[-1],
                                           ei1, ei2)
    if poly.inner is None:
      A = np.vstack([p.T for p in poly.points])
      b = np.ones((len(poly.points), 1))
      inner = cdd.Matrix(np.hstack((b, A)))
      inner.rep_type = cdd.RepType.GENERATOR
    else:
      self.doublepoly.addPointToInner(poly.points[-1], self.edge_i)

    if self.doublepoly is None:
      self.doublepoly = DoublePolygon(inner, outer)
      poly.outer = self.doublepoly.outer
      poly.inner = self.doublepoly.inner
    else:
      self.doublepoly.addCorrespondence()

  def find_direction(self, poly, plot=False):
    if self.edge_i is not None:
      self.build_polys(poly)

    volumes = []

    for edge in self.doublepoly.inner.edges:
      i1, i2 = edge[1:]
      i_e1 = self.doublepoly.correspondences[i1]
      i_e2 = self.doublepoly.correspondences[i2]

      i3 = self.doublepoly.outer.commonPoint(i_e1, i_e2)

      v1 = self.doublepoly.inner.vertices[i1, :]
      v2 = self.doublepoly.inner.vertices[i2, :]
      v3 = self.doublepoly.outer.vertices[i3, :]
      A_e = np.vstack([v1, v2, v3])

      x, y = zip(*A_e)
      if plot:
        poly.reset_fig()
        poly.plot_polyhedrons()
        poly.plot_polyhedron(x, y, 'm', 0.5)
        poly.show()

      vol = 0.5*abs(x[0]*y[1] + x[1]*y[2] + x[2]*y[0]
                    - x[1]*y[0] - x[2]*y[1] - x[0]*y[2])
      volumes.append(vol)

    i, a = max(enumerate(volumes), key=lambda x: x[1])
    self.edge_i = i
    return -self.doublepoly.inner.inequalities[self.doublepoly.inner.edges[i][0], 1:]

  def scipy_triangulate_polyhedron(self, poly):
    points = poly.vertices
    ch = ConvexHull(points)
    return [c for c in ch.points.T], ch.simplices

  def scipy_convexify_polyhedron(self, poly):
    points = poly.vertices
    ch = ConvexHull(points)
    return ch.points
