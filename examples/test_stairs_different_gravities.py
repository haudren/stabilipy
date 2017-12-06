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
import stabilipy as stab
import numpy as np
import matplotlib.pyplot as plt
import os

np.set_printoptions(linewidth=1000)

from stairs_contacts import pos, normals

def main():
  global pos, normals

  mu = 0.7
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  poly_3d = stab.StabilityPolygon(100, dimension=3)
  poly_3d.contacts = contacts

  poly_2d = stab.StabilityPolygon(100, dimension=2)
  poly_2d.contacts = contacts

  polys = [poly_2d, poly_3d]

  point = np.array([[0.36510907, 0.31419711, 0.73607441]]).T

  for poly in polys:
    poly.addTorqueConstraint(contacts[-4:-2],
                             point,
                             10*np.ones((3, 1)))

  gravities = np.linspace(0, 2, 10)
  shape = [
              np.array([[-1., 0, 0]]).T,
              np.array([[1., 0, 0]]).T,
              np.array([[0, 1., 0]]).T,
              np.array([[0, -1., 0]]).T
          ]

  for gravity in gravities:
    for poly in polys:
      poly.gravity_envelope = [gravity*s for s in shape]
      poly.compute(stab.Mode.best, maxIter=50, epsilon=1e-2, solver='cdd',
                   plot_init=False, plot_final=False, plot_step=False,
                   plot_direction=False, plot_error=False)
      filename = os.path.join('different_gravities', 'g_{}_{}d'.format(gravity, poly.dimension))
      poly.save_polyhedron(filename)
      poly.plot()
    plt.show()
  return poly

if __name__ == '__main__':
  main()
