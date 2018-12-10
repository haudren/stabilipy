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
from builtins import zip
from builtins import range
import stabilipy as stab
import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import cdd
import os

from pyparma import Polyhedron as PPL_Poly
from fractions import Fraction

#from crapahut_contacts import pos, normals
from crapahut_bricks_contacts import points, normals
import shapes

np.set_printoptions(linewidth=1000)

def main():
  global points, normals
  pos = points

  #bar = sum(pos)/float(len(pos))
  #pos = [p-bar for p in pos]

  mu = 0.7
  frame_dir = 'interp_vs_cstr'
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  poly = stab.StabilityPolygon(100, dimension=2, radius=100.)
  #poly.contacts = contacts

  #poly.reset_fig()
  #poly.plot_contacts()
  #poly.show()

  #poly.make_problem()
  #poly.check_sizes()


  #poly.contacts = contacts[0:4] + contacts[8:]
  poly.contacts = contacts[4:]

  poly.compute(stab.Mode.best, maxIter=50, epsilon=1e-2, solver='plain',
               plot_init=False, plot_final=False, plot_step=False,
               plot_direction=False, plot_error=False)

  p1 = poly.polygon()

  poly.reset()
  poly.contacts = contacts
  #poly.contacts += [stab.Contact(mu, p, n) for p, n in zip(*lfc())]

  #point = np.array([[0.36510907, 0.31419711, 0.73607441]]).T

  #poly.addTorqueConstraint(contacts[-4:-2],
  #                         point,
  #                         10*np.ones((3, 1)))

  #gravities = np.linspace(0, 2, 10)
  #shape = [
  #            np.array([[-1., 0, 0]]).T,
  #            np.array([[1., 0, 0]]).T,
  #            np.array([[0, 1., 0]]).T,
  #            np.array([[0, -1., 0]]).T
  #        ]

  #for gravity in gravities:
  #  poly.gravity_envelope = [gravity*s for s in shape]

  for k in range(0, 4):
    poly.addForceConstraint(contacts[k:k+1], 1.0)

  poly.compute(stab.Mode.best, maxIter=50, epsilon=1e-2, solver='plain',
               plot_init=False, plot_final=False, plot_step=False,
               plot_direction=False, plot_error=False)

  p0 = poly.polygon()
  interp = shapes.PolygonInterpolator(p0, p1)

  for i in range(3):
    x, y = list(zip(*[[c.r[0], c.r[1]] for c in poly.contacts[4*i:4*(i+1)]]))
    print(x)
    x, y = list(x), list(y)
    x.append(x[0])
    y.append(y[0])
    plt.plot(x, y)
    with open(os.path.join(frame_dir, 'contact_{}'.format(i)), 'w') as f:
      for xi, yi in zip(x, y):
        f.write('{} {}{}'.format(xi[0], yi[0], os.linesep))

  plt.plot(*p1.exterior.coords.xy)
  plt.plot(*p0.exterior.coords.xy)
  plt.show()

  #poly.addTorqueConstraint(contacts[-2:],
  #                         point,
  #                         10*np.ones((3, 1)))
  #poly.reset_fig()
  #poly.plot_contacts()

  A = np.array(p0.exterior.coords)
  np.savetxt(os.path.join(frame_dir, 'p0'), A)
  A = np.array(p1.exterior.coords)
  np.savetxt(os.path.join(frame_dir, 'p1'), A)

  f = open('top_lel.txt', 'w')
  nr_tests = 400
  plot_stuff = False
  for i in range(nr_tests+1):
    poly.clearConstraints()
    cur_percent = (nr_tests-i)/nr_tests
    for k in range(0, 4):
      poly.addForceConstraint(contacts[k:k+1], cur_percent)

    poly.compute(stab.Mode.best, maxIter=50, epsilon=1e-2, solver='plain',
                 plot_init=False, plot_final=False, plot_step=False,
                 plot_direction=False, plot_error=False)

    cur_p = poly.polygon()

    fname = 'cstr_{}'.format(nr_tests-i)
    A = np.array(cur_p.exterior.coords)
    np.savetxt(os.path.join(frame_dir, fname), A)

    cur_area = cur_p.area
    max_j = 1.
    min_j = 0.
    j = 0.5

    pi = interp.fast_interpolate(j)
    interp_area = pi.area

    nr_iter = 0
    while abs(cur_area - interp_area) > 1e-12 and nr_iter < 1000:
      if cur_area > interp_area:
        max_j = j
        j = (min_j+j)/2
      else:
        min_j = j
        j = (max_j+j)/2
      pi = interp.fast_interpolate(j)
      interp_area = pi.area
      #print cur_area - interp_area, max_j, j, min_j
      nr_iter += 1
    f.write("{} {} {} {}\n".format(cur_percent, j, cur_area, interp_area))

    fname = 'interp_{}'.format(nr_tests-i)
    A = np.array(pi.exterior.coords)
    np.savetxt(os.path.join(frame_dir, fname), A)

    #To watch the log file grow
    f.flush()
    print(cur_area - interp_area)

    if plot_stuff:
      poly.plot()
      poly.ax.plot(*pi.exterior.coords.xy)
      poly.show()

    #plt.plot(*pi.exterior.coords.xy)
    #plt.show()
  f.close()

  ##A = np.vstack([p.T for p in poly.points])
  #print zip(*A)
  #qhull = ConvexHull(A)
  #coords, tri = [c for c in qhull.points.T], qhull.simplices

  #dirs = np.vstack(poly.directions[:-1])
  #offsets = np.vstack(poly.offsets[:-1])
  #mat = np.hstack([offsets, -dirs])
  #out = cdd.Matrix(mat, number_type='fraction')
  #print np.array(out)
  #out.rep_type = cdd.RepType.INEQUALITY
  #polyhedron = cdd.Polyhedron(out)
  #gen = np.array(polyhedron.get_generators())

  #f = np.vectorize(lambda x: Fraction(str(x)).limit_denominator(1000))
  #ppl_poly = PPL_Poly(hrep=f(mat))
  #p3 = ppl_poly.vrep()[:, 1:]
  #qhull = ConvexHull(p3)
  #coords3, tri3 = [c for c in qhull.points.T], qhull.simplices

  #print np.array(out).shape
  #print np.array(poly.outer).shape
  ##print np.array(out) - np.vstack([np.array(poly.outer), np.zeros((1, 4))])

  #p2 = gen[:, 1:]
  #qhull = ConvexHull(p2)
  #coords2, tri2 = [c for c in qhull.points.T], qhull.simplices

  #plt.close("all")
  #fig = plt.figure()
  #ax = fig.add_subplot('111', aspect='equal', projection='3d')
  #x, y, z = zip(*A)
  #ax.scatter(x, y, z)
  #x, y, z = zip(*p3)
  #ax.scatter(x, y, z, color='b', marker='^', s=120.)

  #ax.plot_trisurf(*coords, triangles=tri, color='r', alpha=0.5)

  #ax.plot_trisurf(*coords2, triangles=tri2, color='b', alpha=0.1)

  #ax.plot_trisurf(*coords3, triangles=tri3, color='g', alpha=0.1)

  #plt.show()

  return poly

if __name__ == '__main__':
  main()
