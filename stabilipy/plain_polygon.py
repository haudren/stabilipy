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
from __future__ import print_function
from builtins import object
import cdd
import numpy as np
import scipy.spatial as spatial

class DoublePolygon(object):

  """DoublePolygon holds the necessary information to compute stability polygons"""

  def __init__(self, cddinner, cddouter):
    """Default constructor, uses cdd representation of
    inner/outer for initialization

    cddinner: cdd matrix for the initial inner approximation
    cddouter: cdd matrix for the initial outer approximation"""

    polyinner = cdd.Polyhedron(cddinner)
    polyouter = cdd.Polyhedron(cddouter)

    self.inner = Polygon(np.array(polyinner.get_generators())[:, 1:],
                         np.array(polyinner.get_inequalities()))
    self.outer = Polygon(np.array(polyouter.get_generators())[:, 1:],
                         np.array(polyouter.get_inequalities()))

    self.correspondences = [None]*self.inner.vertices.shape[0]

    for i, v in enumerate(self.inner.vertices):
      dists = []
      for j, e in enumerate(self.outer.edges):
        ineq = self.outer.inequalities[e[0], :]
        dist = abs(ineq[0] + v.dot(ineq[1:]))
        dists.append((j, dist))

      j_min, dist = min(dists, key=lambda x: x[1])
      assert dist < 1e-7, "No inequality matches this vertex : {}".format(v)
      self.correspondences[i] = j_min

    print(self.correspondences)
    assert(all([c is not None for c in self.correspondences]))

  def addPointToInner(self, point, edge_index):
    """ Add a point to the inner approximation. This point was found
        in the direction perpendicular to the edge argument"""
    edge = self.inner.edges[edge_index]
    new_i = self.inner.vertices.shape[0]
    self.inner.vertices = np.vstack([self.inner.vertices, point.T])

    self.inner.addEdgeFromPoints(edge[1], new_i, edge[0])
    self.inner.addEdgeFromPoints(new_i, edge[2], None)

  def addInequalityToOuter(self, offset, direction, i1, i2):
    e1 = self.outer.edges[i1]
    e2 = self.outer.edges[i2]

    if e1[1] == e2[2]:
      v_i = e1[1]
      next_edge = i1
      prev_edge = i2
    elif e1[2] == e2[1]:
      v_i = e1[2]
      next_edge = i2
      prev_edge = i1
    else:
      raise ValueError("Edges {}, {} do not share a common opposite vertex".format(e1, e2))

    new_i = self.outer.inequalities.shape[0]

    #print offset, direction

    self.outer.inequalities = np.vstack([self.outer.inequalities,
                                         np.hstack([offset, direction])])

    #Add an edge for the new inequality: both indexes will be updated
    new_ei = len(self.outer.edges)
    self.outer.edges.append([new_i, 0, 0])

    self.outer.addVertexFromEdges(prev_edge, new_ei, v_i)
    self.outer.addVertexFromEdges(new_ei, next_edge, None)

  def addCorrespondence(self):
    """Update correspondence table by adding a last vertex of inner
    to last edge of outer correspondence"""
    self.correspondences.append(len(self.outer.edges) - 1)

class Polygon(object):

  """Polygon class for use by DoublePolygon"""

  def __init__(self, points, inequalities):
    """Initialize the polygon from points and inequalities."""

    qhull = spatial.ConvexHull(points)
    ordered_points = [points[i, :] for i in qhull.vertices]
    self.vertices = np.vstack(ordered_points)
    self.inequalities = inequalities
    self.edges = []
    for j, ineq in enumerate(inequalities):
      extremals = []
      for i, vertex in enumerate(self.vertices):
        if abs(ineq[0] + vertex.dot(ineq[1:])) < 1e-7:
          extremals.append(i)

      assert(len(extremals) == 2)

      if (extremals[1] - extremals[0]) == 1:
        self.edges.append([j, extremals[0], extremals[1]])
      else:
        self.edges.append([j, extremals[1], extremals[0]])

  def addEdgeFromPoints(self, i1, i2, index):
    p1 = self.vertices[i1, :].reshape((2, 1))
    p2 = self.vertices[i2, :].reshape((2, 1))

    if index is None:
      self.edges.append([len(self.edges), i1, i2])
    else:
      self.edges[index] = [index, i1, i2]

    vec = (p2 - p1)/np.linalg.norm(p2-p1)

    normal = np.array([[vec.item(1)],
                       [-vec.item(0)]])

    if index is None:
      self.inequalities = np.vstack([self.inequalities,
                                     np.hstack([normal.T.dot(p1), -normal.T])])
    else:
      self.inequalities[index, :] = np.hstack([normal.T.dot(p1), -normal.T])

  def addVertexFromEdges(self, i1, i2, index):
    v = self.intersection(self.edges[i1][0],
                          self.edges[i2][0])
    if index is None:
      self.vertices = np.vstack([self.vertices, v.T])
      #print i1, i2
      #print self.edges[i1], self.edges[i2]
      self.edges[i1][2] = self.vertices.shape[0]-1
      self.edges[i2][1] = self.vertices.shape[0]-1
    else:
      self.vertices[index, :] = v.T
      self.edges[i1][2] = index
      self.edges[i2][1] = index

  def intersection(self, i1, i2):
    ineq1 = self.inequalities[i1, :]
    ineq2 = self.inequalities[i2, :]

    b1, a1 = ineq1[0], ineq1[1:]
    b2, a2 = ineq2[0], ineq2[1:]

    A = np.vstack([a1, a2])
    b = np.vstack([-b1, -b2])
    x = np.linalg.lstsq(A, b)
    #print "ineqs"
    #print ineq1
    #print ineq2
    #print "A, b + x"
    #print np.hstack([A, b])
    #print x
    #print A.dot(x[0]) - b
    return x[0]

  def commonPoint(self, i1, i2):
    e1 = self.edges[i1]
    e2 = self.edges[i2]

    if e1[1] == e1[2]:
      return e1[1]
    elif e2[1] == e1[2]:
      return e2[1]
    else:
      raise ValueError("Edges do not share a common vertex")
