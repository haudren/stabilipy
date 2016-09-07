#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy.spatial import ConvexHull, HalfspaceIntersection
from scipy.spatial.qhull import QhullError
from scipy.optimize import linprog
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

def cdd_invalidate_vreps(backend, poly):
  offset = poly.offsets[-1]
  direction = poly.directions[-1]
  keys = []

  for key, vrep in poly.vrep_dic.iteritems():
    try:
      valid = all(offset+vrep[:, 1:].dot(direction.T) < 0)
    except IndexError as e:
      print "Error :( "
      valid = True
    if not valid:
      keys.append(key)
      if key in poly.hrep_dic:
        poly.hrep_dic[key].extend(np.hstack((offset, -direction)))

  for key in keys:
    if key in poly.hrep_dic:
      A_e = poly.hrep_dic[key]
    else:
      A_e = poly.outer.copy()
      A_e.extend(cdd.Matrix(-line.reshape(1, line.size)))
      A_e.canonicalize()
      poly.hrep_dic[key] = A_e
      vol = backend.volume_convex(A_e)
      poly.volume_dic[key] = vol
      volumes.append(vol)
      poly.vrep_dic[key] = np.array(cdd.Polyhedron(A_e).get_generators())

def parma_invalidate_vreps(backend, poly):
  offset = poly.offsets[-1]
  direction = poly.directions[-1]
  keys = []

  for key, vrep in poly.vrep_dic.iteritems():
    try:
      valid = all(offset+vrep[:, 1:].dot(direction.T) < 0)
    except IndexError as e:
      print "Error :( "
      valid = True
    if not valid:
      keys.append(key)
      if key in poly.hrep_dic:
        poly.hrep_dic[key].add_ineq(fractionize(np.hstack((offset, -direction))))

  for key in keys:
    if not key in poly.hrep_dic:
      A_e = poly.outer.copy()
      A_e.add_ineq(-line)
      vol = self.volume_convex(A_e)
      poly.volume_dic[key] = vol
      volumes.append(vol)
      poly.vrep_dic[key] = A_e.vrep()
      poly.hrep_dic[key] = A_e

def capping_invalidate_vreps(poly):
  offset = poly.offsets[-1]
  direction = poly.directions[-1]

  keys = []

  for key, vrep in poly.hull_dic.iteritems():
    valid = all(offset+vrep.points.dot(direction.T) < 0)
    if not valid:
      keys.append(key)

  for key in keys:
    hs = np.zeros((1, poly.outer.halfspaces.shape[1]))
    hs[:, :-1] = poly.directions[-1]
    hs[:, -1] = -poly.offsets[-1]
    try:
      poly.hrep_dic[key].add_halfspaces(hs)
      poly.hull_dic[key] = ConvexHull(poly.hrep_dic[key].intersections)
    except QhullError:
      del poly.hrep_dic[key]
      del poly.hull_dic[key]

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
    self.last_hrep = None

    if geomengine == 'scipy':
      self.volume_convex = self.scipy_volume_convex

  def scipy_volume_convex(self, hrep):
    try:
      points = np.array(cdd.Polyhedron(hrep).get_generators())[:, 1:]
    except RuntimeError:
      return 0
    except IndexError:
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
      if self.last_hrep is not None:
        self.last_hrep.extend(np.hstack((poly.offsets[-1], -poly.directions[-1])))
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
    self.last_hrep = poly.hrep_dic[key]
    return -ineq[i, 1:]

  def invalidate_vreps(self, poly):
    cdd_invalidate_vreps(self, poly)

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
        poly.hrep_dic[key] = A_e
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

  def invalidate_vreps(self, poly):
    parma_invalidate_vreps(self, poly)

  def scipy_triangulate_polyhedron(self, poly):
    points = floatize(poly.vrep()[:, 1:])
    ch = ConvexHull(points)
    return [c for c in ch.points.T], ch.simplices

  def scipy_convexify_polyhedron(self, poly):
    points = floatize(poly.vrep()[:, 1:])
    ch = ConvexHull(points)
    return ch.points

class QhullBackend(object):

  """Using the Qhull backend for polygon computation. This is an experimental backend
  that sould yield better performance. Works on floating-point numbers. Requires scipy."""

  def __init__(self, geomengine='scipy'):
    """Default constructor.

    :param geomengine: underlying geometry engine. Only scipy is supported
    """

    self.name = 'qhull'
    self.last_hrep = None
    self.feasible_point = None
    self.feasible_point_str = None

    if geomengine == 'scipy':
      self.volume_convex = self.scipy_volume_convex

  def scipy_volume_convex(self, hrep):
    if isinstance(hrep, HalfspaceIntersection):
      return ConvexHull(hrep.intersections).volume
    return hrep.volume

  def build_polys(self, poly):
    if poly.inner is None:
      A = np.vstack([p.T for p in poly.points])
      poly.inner = ConvexHull(A, incremental=True)
      self.feasible_point = np.sum(A, axis=0)/A.shape[0]
    else:
      try:
        poly.inner.add_points(poly.points[-1].T, restart=False)
      except QhullError:
        print "Incremental processing failed"
        A = np.vstack([p.T for p in poly.points])
        poly.inner = ConvexHull(A, incremental=True)

    if poly.outer is None:
      if poly.inner is None:
        raise AssertionError("Inner polygon should have been initialized")
      A = np.vstack(poly.directions)
      b = np.vstack(poly.offsets)
      poly.outer = HalfspaceIntersection(np.hstack((A, -b)), self.feasible_point, incremental=True)
    else:
      hs = np.zeros((1, poly.outer.halfspaces.shape[1]))
      hs[:, :-1] = poly.directions[-1]
      hs[:, -1] = -poly.offsets[-1]
      poly.outer.add_halfspaces(hs)

  def find_direction(self, poly, plot=False):
    self.build_polys(poly)

    volumes = []

    ineq = poly.inner.equations

    for line in ineq:
      key = hashlib.sha1(line).hexdigest()
      if key in poly.hull_dic:
        volumes.append(poly.hull_dic[key].volume)
      else:
        #TODO: How to compute outer initial vrep ?
        # We use CDD for now, but is there any way to get
        # initial feasible point ?

        #A_e = cdd.Matrix(np.hstack((-poly.outer.equations[:, -1:],
        #                            -poly.outer.equations[:, :-1])))
        #A_e.rep_type = cdd.RepType.INEQUALITY
        #m = np.hstack((line[-1:], line[:-1]))
        #A_e.extend(cdd.Matrix(m.reshape((1, m.size))))
        #points = np.array(cdd.Polyhedron(A_e).get_generators())[:, 1:]

        A_e = np.vstack((poly.outer.halfspaces, -line))
        c = np.zeros((A_e.shape[1],))
        c[-1] = -1
        res = linprog(c, A_ub=np.hstack((A_e[:, :-1], np.ones((A_e.shape[0], 1)))),
            b_ub=-A_e[:, -1:], bounds=(None, None))
        if res.success:
          feasible_point = res.x[:-1]
          try:
            poly.hrep_dic[key] = HalfspaceIntersection(A_e, feasible_point, incremental=True)
            poly.hull_dic[key] = ConvexHull(poly.hrep_dic[key].intersections)
            volumes.append(poly.hull_dic[key].volume)
          except QhullError:
            volumes.append(0.)
        else:
          print "Error : {}".format(res.message)
          volumes.append(0.)

        if plot and res.success:
          poly.reset_fig()
          poly.plot_polyhedrons()
          poly.plot_polyhedron(poly.hull_dic[key], 'm', 0.5)
          poly.show()

    maxv = max(volumes)
    alli = [i for i, v in enumerate(volumes) if v == maxv]
    i = random.choice(alli)
    return ineq[i, :-1]

  def invalidate_vreps(self, poly):
    capping_invalidate_vreps(poly)

  def scipy_triangulate_polyhedron(self, hrep):
    vrep = hrep
    if isinstance(hrep, HalfspaceIntersection):
      vrep = ConvexHull(hrep.intersections)
    return [c for c in vrep.points.T], vrep.simplices

  def scipy_convexify_polyhedron(self, hrep):
    vrep = hrep
    if isinstance(hrep, HalfspaceIntersection):
      vrep = ConvexHull(hrep.intersections)
    return vrep.points


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
