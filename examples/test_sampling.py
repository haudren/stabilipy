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
from builtins import map
from builtins import zip
import numpy as np
import stabilipy as stab

import matplotlib.pyplot as plt
import cdd

#normals = map(stab.normalize,
#              [np.array([[0, 0, 1]]).T,
#               np.array([[0, 0, 1]]).T,
#               np.array([[1, 0, 0]]).T])
#

def main():
  #normals = map(stab.normalize,
  #              [np.array([[-0.7, 0, 1]]).T,
  #               np.array([[0.1, 0.5, 1]]).T,
  #               np.array([[0.5, -0.5, 1]]).T])

  #pos = [np.array([[0, 0, 0]]).T,
  #       np.array([[0, 1, 0]]).T,
  #       np.array([[1, 0, 0.2]]).T]

  #normals = map(stab.normalize,
  #              [np.array([[0., 0., 1.]]).T,
  #               np.array([[0., 0., 1.]]).T,
  #               np.array([[0., 0., 1.]]).T,
  #               np.array([[0., 0., 1.]]).T])

  #pos = [np.array([[0.1, 0.155, -0.105]]).T,
  #       np.array([[-0.07, 0.155, -0.105]]).T,
  #       np.array([[-0.07, 0.055, -0.105]]).T,
  #       np.array([[0.1, 0.055, -0.105]]).T]

  #normals = [np.array([[0, 0, 1]]).T]*3

  normals = list(map(stab.normalize,
                [np.array([[0.5, 0.0, 1]]).T,
                 np.array([[0.5, 0.0, 1]]).T,
                 np.array([[0.0, 0.0, 1]]).T]))

  pos = [np.array([[0.1, -0.2, 0]]).T,
         np.array([[-0.1, -0.2, 0]]).T,
         np.array([[0.2, 0, 1]]).T]

  bar = sum(pos)/float(len(pos))
  pos = [p-bar for p in pos]

  mu = 0.5
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  poly = stab.StabilityPolygon(60)
  poly.contacts = contacts

  poly.compute(stab.Mode.precision, epsilon=2e-2,
               record_anim=False, plot_final=False)

  polyhedron = poly.polyhedron()

  m = np.min(polyhedron, axis=0)
  v = np.max(polyhedron, axis=0) - np.min(polyhedron, axis=0)

  points = np.random.rand(10000, 3)
  points = points*v+m
  stability = np.array([poly.sample(p.T[:, np.newaxis]) for p in points])

  poly.reset_fig()
  poly.plot()

  #fig = plt.figure()
  #ax = fig.add_subplot('121', aspect='equal', projection='3d')
  stable = points[stability]
  unstable = points[~stability]

  #x, y, z = zip(*stable)
  #poly.ax.scatter(x, y, z, c='g', marker='o')
  #ax = fig.add_subplot('122', aspect='equal', projection='3d')
  if unstable.size > 0:
    pass
    #x, y, z = zip(*unstable)
    #poly.ax.scatter(x, y, z, c='r', marker='^')

  vol_in = stab.volume_convex(poly.inner)
  vol_out = stab.volume_convex(poly.outer)
  vol_sampling = v[0]*v[1]*v[2]

  ineq_in = np.array(cdd.Polyhedron(poly.inner).get_inequalities())
  ineq_out = np.array(cdd.Polyhedron(poly.outer).get_inequalities())

  #H format : b - A x > 0 stored as b | -A
  test_in_points = ineq_in[:, [0]] + ineq_in[:, 1:].dot(points.T)
  in_points = np.all(test_in_points > 0, axis=0)

  test_out_points = ineq_out[:, [0]] + ineq_out[:, 1:].dot(points.T)
  out_points = np.all(test_out_points > 0, axis=0)

  und_points = np.logical_and(out_points, np.logical_not(in_points))
  criterion = np.logical_or(und_points, np.equal(in_points, stability))

  c = np.all(criterion)
  print("Criterion : {}".format(c))
  if not c:
    print("Criterion is false...")
    print("Number of failures : {}".format(
          criterion.size - np.count_nonzero(criterion)))
    print("Failure points : {}".format(points.T[~criterion]))
  else:
    print("Hurray all decidable points are confirmed!")

  print("Volume inner : {} | Volume outer : {} | Volume sampling : {}".format(
        vol_in, vol_out, vol_sampling))
  print("Ratio vol. in./sampling : {} | Ratio points in/total : {}".format(
        vol_in/vol_sampling, stable.shape[0]/points.shape[0]))
  print("Ratio vol. undecidable/sampling : {} | Ratio undec. points/total: {}".format(
        (vol_out - vol_in)/vol_sampling, np.count_nonzero(und_points)/points.shape[0]))

  plt.show()

main()
