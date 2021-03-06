from __future__ import print_function
from builtins import zip
from builtins import map
from builtins import range
from builtins import object
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


from polygon_contacts import contacts

from scipy.spatial import ConvexHull

import matplotlib.pyplot as plt

azim = 48.9035087719
elev = 31.5350877193
xlim = [-0.95389899, 0.95389899]
ylim = [-0.95389899, 0.95389899]
zlim = [-0.95389899, 0.95389899]


class CircularBuffer(object):

  """Holds object and allows to access from both ends"""

  def __init__(self, l):
    """Construct from list

    :l: A list

    """

    self.l = l

  def __getitem__(self, i):
    if not isinstance(i, int):
      raise KeyError("Only integers are accepted in circular buffer")
    if i < len(self.l):
      return self.l[i]
    else:
      return self.l[i % len(self.l)]

def point_seq(n):
  l = [None]*(2*n)
  l[::2] = list(range(n))
  l[1::2] = [-i for i in range(n)]
  return l[1:]

def main(n, alphas):
  if len(alphas) != n:
    raise RuntimeError("Need same number of constraints as contacts")

  if sum(alphas) < 1:
    raise RuntimeError("Bound too low")

  mu = 0.5
  pos, normals = contacts(n)
  cont = [stab.Contact(mu, pi, ni) for pi, ni in zip(pos, normals)]

  polygon = stab.StabilityPolygon(200, dimension=2, radius=1.5)
  polygon.contacts = cont

  for i in range(n):
    polygon.addForceConstraint([polygon.contacts[i]], alphas[i])

  polygon.compute(stab.Mode.best, epsilon=2e-3, maxIter=100, solver='plain',
                  record_anim=False, plot_init=False,
                  plot_step=False, plot_final=False)
  polygon.reset_fig()
  polygon.plot_contacts()

  points = polygon.inner.vertices
  hull = ConvexHull(points)

  vertices = np.hstack((hull.vertices, hull.vertices[0]))
  polygon.ax.plot(points[vertices,0], points[vertices,1], 'r--', lw=2)

  shrink_points = []
  seq = CircularBuffer(point_seq(n))
  alphas = CircularBuffer(alphas)
  points = CircularBuffer([c.r[0:2, :] for c in polygon.contacts])
  for i in range(n):
    a0 = a1 = alphas[i]
    r0 = r1 = alphas[i]*points[i]
    k0 = k1 = i
    ex0 = "{}*p{}".format(alphas[i], i)
    ex1 = "{}*p{}".format(alphas[i], i)

    done0 = done1 = False
    niter = 1
    while not (done0 and done1):
      k0 = i + seq[niter]
      k1 = i - seq[niter]

      if not done0:
        if a0 + alphas[k0] >= 1:
          r0 = r0 + (1-a0)*points[k0]
          ex0 += "+ {}*p{}".format(1-a0, k0)
          done0 = True
        else:
          r0 = r0 + alphas[k0]*points[k0]
          ex0 += "+ {}*p{}".format(alphas[k0], k0)
          a0 += alphas[k0]

      if not done1:
        if a1 + alphas[k1] >= 1:
          r1 = r1 + (1-a1)*points[k1]
          ex1 += " {}*p{}".format(1-a1, k1)
          done1 = True
        else:
          r1 = r1 + alphas[k1]*points[k1]
          ex1 += " {}*p{}".format(alphas[k1], k1)
          a1 += alphas[k1]
      niter += 1

    print(ex0)
    print(ex1)
    shrink_points.append(r0)
    shrink_points.append(r1)
    polygon.ax.plot(r0[0], r0[1], marker='^', markersize=20)
    polygon.ax.plot(r1[0], r1[1], marker='^', markersize=20)

  polygon.show()

print("n-sided gon : {} with limits {}".format(int(sys.argv[1]), list(map(float, sys.argv[2:]))))

main(int(sys.argv[1]), list(map(float, sys.argv[2:])))
