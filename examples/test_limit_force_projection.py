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

from __future__ import print_function
from builtins import zip
from builtins import range
import stabilipy as stab
import numpy as np
import cdd
import copy
from stabilipy.linear_cone import build_cone
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

from stairs_contact import pos, normals

np.set_printoptions(linewidth=1000000, precision=2, threshold=1e12)


def main():
  global pos, normals

  #bar = sum(pos)/float(len(pos))
  #pos = [p-bar for p in pos]

  mu = 0.2
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  poly = stab.StabilityPolygon(60, dimension=2,
                               force_lim=1000.)

  poly.contacts = contacts[0:8]

  for c in poly.contacts[0:4]:
    c.r[2] = 0.
    c.n = np.array([[0.4], [0.4], [np.sqrt(1-2*(0.4**2))]])
    #c.n = np.array([[0., 1., 0.]]).T

  for c in poly.contacts[4:]:
    c.r[2] = 0.
    #c.n = np.array([[0., 0., 1.]]).T
    c.n = np.array([[-0.4], [-0.4], [np.sqrt(1-2*(0.4**2))]])

  poly.reset_fig()
  poly.plot_contacts()
  poly.show()

  sol = 'plain'

  #Compute the unconstrained and save ineqs
  poly.compute(stab.Mode.iteration, maxIter=20, epsilon=2e-3,
               solver=sol, plot_error=False, plot_step=False,
               plot_init=False, plot_final=False)
  poly_ineq = poly.backend.doublepoly.inner.inequalities

  radius = 0.11
  fc = 4
  nc = 8

  for c in range(fc, nc):
    poly.addForceConstraint([poly.contacts[c]], radius)

  poly.compute(stab.Mode.iteration, maxIter=20, epsilon=2e-3,
               solver=sol, plot_error=False, plot_step=False,
               plot_init=False, plot_final=False)

  poly.plot()
  poly.show()

  assert(poly.dimension == 2)
  assert(len(poly.gravity_envelope) == 1)

  A1, A2, t = poly.A1, poly.A2, poly.t

  sphere_ineq = np.array([[1., -1., 1.],
                          [1., -1., -1.],
                          [1., 1., -1.],
                          [1., 1., 1.],
                          [-1., 1., -1.],
                          [-1., 1., 1.],
                          [-1., -1., -1.],
                          [-1., -1., 1.]])

  sphere = np.zeros((8*(nc-fc), 1+poly.nrVars()))
  for contact_id in range(fc, nc):
    line = 8*(contact_id-fc)
    col = 1+3*contact_id
    sphere[line:line+8, col:col+3] = sphere_ineq
    sphere[line:line+8, 0] = radius*poly.mass*9.81

  nr_lines = poly_ineq.shape[0]
  exp_poly_ineq = np.hstack([poly_ineq[:, 0:1], np.zeros((nr_lines, poly.nrVars()-2)), poly_ineq[:, 1  :]])

  eq = np.hstack((t, -A1, -A2))

  mat = cdd.Matrix(sphere, number_type='fraction')
  mat.rep_type = cdd.RepType.INEQUALITY
  mat.extend(exp_poly_ineq)
  mat.extend(eq, linear=True)
  print("Let's goooooo")
  cdd_poly = cdd.Polyhedron(mat)
  vertices = np.array(cdd_poly.get_generators())

  print(vertices.shape)
  if len(vertices.shape) > 1:
    point_mask = vertices[:, 0] == 1
    points = vertices[point_mask, -2:]
    rays = vertices[~point_mask, -2:]
    hull = ConvexHull(points)

    poly.reset_fig()
    poly.plot_contacts()
    poly.plot_polyhedron(poly.inner, 'blue', 0.5)
    poly.ax.plot(hull.points[hull.vertices, 0],
                 hull.points[hull.vertices, 1],
                 'red', label='cdd', marker='^', markersize=10)
    for ray in rays:
      if np.linalg.norm(ray) > 1e-10:
        print(ray)
        pp = np.vstack([i*ray for i in np.linspace(0.01, 1)])
        poly.ax.plot(pp[:, 0], pp[:, 1], 'red')
      else:
        print("This is a zero ray")
    poly.show()
  else:
    print("No vertices")

  #if sol == 'plain':
  #  ineq = [l/abs(l[0]) for l in poly.inner.inequalities]
  #  print np.vstack(ineq)
  #elif sol == 'cdd':
  #  poly = cdd.Polyhedron(poly.inner)
  #  print poly.get_inequalities()
  #return poly

if __name__ == '__main__':
  main()
