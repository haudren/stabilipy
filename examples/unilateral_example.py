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

import stabilipy as stab
import numpy as np
import sys

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

  polyhedron = stab.StabilityPolygon(200, dimension=3, radius=1.5, robust_sphere=False)
  polyhedron.contacts = contacts

  polygon = stab.StabilityPolygon(200, dimension=2, radius=1.5)
  polygon.contacts = contacts

  shape = [
              np.array([[-1., 0, 0]]).T,
              np.array([[1., 0, 0]]).T,
              np.array([[0, 1., 0]]).T,
              np.array([[0, -1., 0]]).T,
              np.array([[0, 0., 1]]).T,
              np.array([[0, 0., -1]]).T
          ]

  polytope = [margin*s for s in shape]

  polyhedron.gravity_envelope = polytope
  polyhedron.compute(stab.Mode.iteration, epsilon=2e-3, maxIter=10, solver='qhull',
                     record_anim=False, plot_init=False,
                     plot_step=False, plot_final=True)
  #polyhedron.reset_fig()
  #polyhedron.ax.set_xlabel("x(m)")
  #polyhedron.ax.set_ylabel("y(m)")
  #polyhedron.ax.set_zlabel("z(m)")
  #polyhedron.ax.view_init(elev=elev, azim=azim)
  #polyhedron.ax.set_xlim3d(*xlim)
  #polyhedron.ax.set_ylim3d(*ylim)
  #polyhedron.ax.set_zlim3d(*zlim)
  #polyhedron.plot_contacts()
  #polyhedron.plot_solution()
  #polyhedron.plot_polyhedrons()
  #polyhedron.show()

  #plt.savefig('{}.png'.format(margin))

  #polygon.gravity_envelope = polytope
  #polygon.compute(stab.Mode.best, epsilon=2e-3, maxIter=20, solver='parma',
  #                record_anim=False, plot_init=False,
  #                plot_step=False, plot_final=False)

  #print polyhedron.volume_convex(polyhedron.inner)
  #print 3*polygon.volume_convex(polygon.inner)

print("Margin : {}".format(sys.argv[1]))

main(float(sys.argv[1]))
