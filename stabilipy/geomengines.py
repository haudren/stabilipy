#!/usr/bin/env python2
# -*- coding: utf-8 -*-

# Copyright 2015-2017 CNRS-AIST JRL

# This file is part of StabiliPy.

# StabiliPy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# StabiliPy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with StabiliPy.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division
from builtins import zip
from builtins import range
import numpy as np
import cdd
from scipy.spatial import ConvexHull
from functools import reduce

try:
  import shapely.geometry as geom
  SHAPELY_AVAILABLE = True
except ImportError:
  SHAPELY_AVAILABLE = False

try:
  from CGAL.CGAL_Polyhedron_3 import Polyhedron_3
  from CGAL.CGAL_Kernel import Point_3, Tetrahedron_3
  from CGAL.CGAL_Convex_hull_3 import convex_hull_3
  CGAL_AVAILABLE = True
except ImportError:
  CGAL_AVAILABLE = False

try:
  from qhull_sch import convexVolume, convexHull
  QHULL_SCH_AVAILABLE = True
except ImportError:
  QHULL_SCH_AVAILABLE = False

def area_convex(hrep):
  try:
    gen = np.array(cdd.Polyhedron(hrep).get_generators())
  except RuntimeError:  # Whenever numerical inconsistency is found
    return 0
  #If the polygon is empty or degenerate, return 0
  if gen.shape[0] < 3:
    return 0
  p = np.vstack([gen, gen[0, :]])
  x, y = p[:, 1], p[:, 2]
  poly = geom.Polygon(list(zip(x, y)))
  return poly.convex_hull.area

def convex_hull(x_p, y_p):
  poly = geom.Polygon(list(zip(x_p, y_p)))
  x_c, y_c = poly.convex_hull.exterior.xy
  return x_c, y_c

def add_cgal(p1, p2):
  return Point_3(p1.x() + p2.x(),
                 p1.y() + p2.y(),
                 p1.z() + p2.z())

def mult_cgal(p, a):
  return Point_3(a*p.x(), a*p.y(), a*p.z())

def as_list(p):
  return [p.x(), p.y(), p.z()]

def convexify_polyhedron(hrep):
  return scipy_convexify_polyhedron(hrep)

def qhull_convexify_polyhedron(hrep):
  gen = np.array(cdd.Polyhedron(hrep).get_generators())
  #If the polygon is empty or degenerate, return 0
  if gen.shape[0] < 3:
    return 0

  points = gen[:, 1:].tolist()
  pout = convexHull(points)
  return np.array(pout)

def scipy_convexify_polyhedron(hrep):
  points = np.array(cdd.Polyhedron(hrep).get_generators())[:, 1:]
  ch = ConvexHull(points, qhull_options='QbB')
  return ch.points

def cgal_convexify_polyhedron(hrep):
  gen = np.array(cdd.Polyhedron(hrep).get_generators())
  #If the polygon is empty or degenerate, return 0
  if gen.shape[0] < 3:
    return 0

  points = [Point_3(x, y, z) for x, y, z in gen[:, 1:]]
  poly = Polyhedron_3()
  convex_hull_3(points, poly)
  lines = [np.array([as_list(p)]) for p in poly.points()]
  return np.vstack(lines)

def tetrahedron_from_facet(facet, center):
  points = []
  he = facet.halfedge()
  points.append(he.vertex().point())
  h = next(he)
  for i in range(facet.facet_degree()-1):
    points.append(h.vertex().point())
    h = next(h)
  points.append(center)
  return Tetrahedron_3(*points)

def qhull_volume_convex(hrep):
  gen = np.array(cdd.Polyhedron(hrep).get_generators())
  #If the polygon is empty or degenerate, return 0
  if gen.shape[0] < 3:
    return 0
  points = gen[:, 1:].tolist()
  return convexVolume(points)

def cgal_volume_convex(hrep):
  gen = np.array(cdd.Polyhedron(hrep).get_generators())
  #If the polygon is empty or degenerate, return 0
  if gen.shape[0] < 3:
    return 0
  #p = np.vstack([gen, gen[0, :]])
  p = gen

  points = [Point_3(x, y, z) for x, y, z in p[:, 1:]]
  poly = Polyhedron_3()
  convex_hull_3(points, poly)
  points = list(poly.points())

  center = mult_cgal(reduce(add_cgal, points),
                     1/float(len(points)))

  tetrahedrons = [tetrahedron_from_facet(f, center) for f in poly.facets()]

  return sum([abs(t.volume()) for t in tetrahedrons])
