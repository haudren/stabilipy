#!/usr/bin/env python
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

import stabilipy as stab
import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection

from scipy.optimize import linprog

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def _make_polygon(mass, gravity, contacts, radius):
  polygon = stab.StabilityPolygon(mass, dimension=2, radius=radius)
  polygon.contacts = contacts
  polygon.gravity_envelope = [gravity]
  return polygon

def feasible_point(A):
  c = np.zeros((A.shape[1],))
  c[-1] = -1

  res = linprog(c, A_ub=np.hstack((A[:, :-1], np.ones((A.shape[0], 1)))),
      b_ub=-A[:, -1:], bounds=(None, None))
  if res.success:
    return res.x[:-1]
  else:
    raise RuntimeError("No interior point found")

def _intersect(hulls, feas_point=None):
  all_eq = np.vstack([h.equations for h in hulls])

  if feas_point is None:
    #points = np.vstack([h.points[h.vertices, :] for h in hulls])
    #feas_point = np.sum(points, axis=0)/points.shape[0]
    feas_point = feasible_point(all_eq)

  intersected = ConvexHull(HalfspaceIntersection(all_eq, feas_point).intersections)
  return intersected

class PrismIntersection():

  def __init__(self, mass, polytope, contacts, radius=1.5):
    self.mass = mass
    self.polytope = polytope
    self.contacts = contacts
    self.radius = radius

    self.initialized = False

    self.polygons = [_make_polygon(mass, s, contacts, self.radius)
                      for s in self.polytope]
    self.prisms = None
    self.figure = None
    self.axes = []
    self.threedax = None

  def initialize(self):
    for polygon in self.polygons:
      polygon.select_solver('plain')
      polygon.make_problem()
      polygon.init_algo()
      polygon.build_polys()

    self.initialized = True

  def addTorqueConstraint(self, begin, end, point, ub, lb=None):
    """Add a limit on torque at a point over a set of contacts defined by
    begin:end """
    for polygon in self.polygons:
      polygon.addTorqueConstraint(polygon.contacts[begin:end], point, ub, lb)

  def addForceConstraint(self, begin, end, flim):
    """Add a limit on sum of forces over a set of contacts defined by
    begin:end """
    for polygon in self.polygons:
      polygon.addForceConstraint(polygon.contacts[begin:end], flim)

  def compute(self, mode, epsilon, maxIter):
    for polygon in self.polygons:
      polygon.compute(mode, epsilon=epsilon, maxIter=maxIter, solver='plain',
          record_anim=False, plot_init=False, plot_step=False, plot_final=False)
    self.initialized = True

  def sample(self, point):
    if not self.initialized:
      self.initialize()

    trans_p = [point[:2] + point.item(2)*np.array([s[0,0], s[1,0]])/(s[2]+9.81) for s in self.polytope]

    iters = 0
    for res in self._sample(trans_p):
      iters += res[1]
      if not res[0]:
        return False, iters

    return True, iters

  def _sample(self, trans_p):
    for pi, polygon in zip(trans_p, self.polygons):
      yield polygon.sample(pi, plot_final=False)

  def polyhedron(self):
    self.prisms = []
    for s, polygon in zip(self.polytope, self.polygons):
      points = ConvexHull(polygon.polyhedron()).points
      origin = np.array([[-s[0]/(s[2]+9.81), -s[1]/(s[2]+9.81), 1]])
      top = np.hstack((points, np.zeros((points.shape[0], 1)))) + self.radius*origin
      bot = np.hstack((points, np.zeros((points.shape[0], 1)))) - self.radius*origin
      self.prisms.append(ConvexHull(np.vstack((top, bot))))

    return _intersect(self.prisms)

  def outer_polyhedron(self):
    self.outer_prisms = []
    for s, polygon in zip(self.polytope, self.polygons):
      points = ConvexHull(polygon.outer_polyhedron()).points
      origin = np.array([[-s[0]/(s[2]+9.81), -s[1]/(s[2]+9.81), 1]])
      top = np.hstack((points, np.zeros((points.shape[0], 1)))) + self.radius*origin
      bot = np.hstack((points, np.zeros((points.shape[0], 1)))) - self.radius*origin
      self.outer_prisms.append(ConvexHull(np.vstack((top, bot))))

    return _intersect(self.outer_prisms)

  def plot(self):
    self.figure = plt.figure()
    nraxes = len(self.polytope)
    nrows = 2
    ncols = nraxes // nrows + 1
    if(nraxes % nrows != 0):
      ncols += 1

    self.axes = [plt.subplot2grid((nrows, ncols), (i % nrows, i // nrows), projection="3d")
        for i in range(nraxes)]

    self.threedax = plt.subplot2grid((nrows, ncols), (0, ncols-1), rowspan=3, projection="3d")

    for polygon, ax in zip(self.polygons, self.axes):
      polygon.ax = ax
      polygon.plot_polyhedrons()

    polyhedron = self.polyhedron()

    coords = [c for c in polyhedron.points.T]
    surf = self.threedax.plot_trisurf(*coords, triangles=polyhedron.simplices, color="red", alpha=0.3, shade=True)
    surf.set_edgecolor('red')

    self.set_axes_properties()

  def set_axes_properties(self):
    for s, ax in zip(self.polytope, self.axes):
      ax.set_aspect("equal")
      ax.set_xlabel("x(m)")
      ax.set_ylabel("y(m)")
      ax.set_title("s : {}".format(s.T))

  def show(self):
    plt.show()
