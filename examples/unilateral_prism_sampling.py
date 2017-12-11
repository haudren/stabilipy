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

  shape = [
              np.array([[-1., 0, 0]]).T,
              np.array([[1., 0, 0]]).T,
              np.array([[0, 1., 0]]).T,
              np.array([[0, -1., 0]]).T,
              #np.array([[1, -1., 0]]).T/np.sqrt(2),
          ]

  polytope = [margin*s for s in shape]

  prisms = stab.PrismIntersection(200, polytope, contacts, radius=1.5)

  points = np.array([1., 1., 3.])*(np.random.random((10**6,3)) - 0.5)
  np.savetxt('random_points.txt', points)
  #points = np.loadtxt('random_points.txt')

  truths, iters = list(zip(*[prisms.sample(point)
                        for point in points]))
  #truths = np.array(truths)

  #prisms.plot()

  print("total iters: {}".format(np.sum(iters)))

  #prisms.threedax.plot(*zip(*points[truths, :]), linestyle="none", marker="*", markerfacecolor="green")
  #prisms.threedax.plot(*zip(*points[~truths, :]), linestyle="none", marker="x", markerfacecolor="red")

  #prisms.show()

print("Margin : {}".format(sys.argv[1]))

main(float(sys.argv[1]))
