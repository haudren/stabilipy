from __future__ import print_function
from builtins import zip
from builtins import range
import stabilipy as stab
import numpy as np
import sys
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


from unilateral_contacts import pos, normals

import matplotlib.pyplot as plt

azim = 48.9035087719
elev = 31.5350877193
xlim = [-0.95389899, 0.95389899]
ylim = [-0.95389899, 0.95389899]
zlim = [-0.95389899, 0.95389899]

def main(margin):
  mu = 0.5
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  contacts[2].mu = 0.5

  polyhedron = stab.StabilityPolygon(200, dimension=3, radius=1.5)
  polyhedron.contacts = contacts

  shape = [
              np.array([[-1., 0, 0]]).T,
              np.array([[1., 0, 0]]).T,
              np.array([[0, 1., 0]]).T,
              np.array([[0, -1., 0]]).T
          ]

  polytope = [margin*s for s in shape]

  polyhedron.gravity_envelope = polytope

  polyhedron.select_solver('qhull')
  polyhedron.make_problem()
  polyhedron.init_algo()
  polyhedron.build_polys()

  #polyhedron2 = stab.StabilityPolygon(200, dimension=3, radius=1.5)
  #polyhedron2.contacts = contacts
  #polyhedron2.gravity_envelope = polytope
  #polyhedron2.compute(stab.Mode.best, epsilon=1e-3, maxIter=100, solver='cdd', plot_final=False)

  points = np.array([1.,1.,3.])*(np.random.random((10**3,3)) - 0.5)
  #np.savetxt('random_points.txt', points)
  #points = np.loadtxt('random_points.txt')

  truths, nrIter = list(zip(*[polyhedron.sample(point, plot_final=False, plot_step=False)
                      for point in points]))
  truths = np.array(truths)

  print(sum(nrIter))

  polyhedron.plot()
  #polyhedron2.ax = polyhedron.ax
  #polyhedron2.plot_polyhedrons()

  polyhedron.ax.plot(*list(zip(*points[truths, :])), linestyle="none", marker="*", markerfacecolor="green")
  polyhedron.ax.plot(*list(zip(*points[~truths, :])), linestyle="none", marker="x", markerfacecolor="red")

  polyhedron.show()

  A = np.vstack((
      np.array(list(range(points.shape[0]))),
      np.cumsum(np.array(nrIter))))

  np.savetxt('nr_iter.txt', A.T)

print("Margin : {}".format(sys.argv[1]))

main(float(sys.argv[1]))
