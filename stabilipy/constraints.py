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

from __future__ import absolute_import
from builtins import range
from builtins import object
from enum import Enum, unique
import numpy as np
from .utils import cross_m

@unique
class Constraint(Enum):
  """Constraint types. Can only be inequality or conic."""
  Inequality = 1
  Conic = 2

class TorqueConstraint(object):

  """Constraint to limit torque applied on certain contacts"""

  def __init__(self, indexes, point, ub, lb=None):
    """Default constructor.

    :param indexes: Indexes of the contacts on which the constraint applies
    :param point: Point where the torques are computed (3,1) array
    :param ub: Upper bound (3,1) array
    :param lb: Lower bound, can be None (3,1) array

    """
    self.indexes = indexes
    self.point = point
    self.ub = ub

    if lb is None:
      self.lb = -ub
    else:
      self.lb = lb

    self.size = 6
    self.ctype = Constraint.Inequality

  def compute(self, poly):
    nrVars = poly.nrVars()
    L = np.zeros((6, nrVars))
    for i in self.indexes:
      cont = poly.contacts[i]
      dist = self.point - cont.r
      off = 0
      #Add constraint on every x_i
      for j in range(len(poly.gravity_envelope)):
        L[:, off+3*i:off+3*i+3] = np.vstack([cross_m(dist), -cross_m(dist)])
        off += 3*len(poly.contacts)
    #Filter L, tb to remove zero lines
    tb = np.vstack([self.ub, -self.lb])
    zero_mask = np.all(L == 0, axis=1)

    self.L = L[~zero_mask]
    self.tb = tb[~zero_mask]

    self.size = L.shape[0]

  def matrices(self):
    return self.L, self.tb

class ForceConstraint(object):

  """Constraint to limit force applied on certain contacts"""

  def __init__(self, indexes, limit):
    """Default constructor.

    :param indexes: Indexes of the contacts on which the constraint applies
    :param limit: Maximum force, expressed as percentage of robot weight

    """
    self.indexes = indexes
    self.limit = limit

    self.size = 4
    self.ctype = Constraint.Conic

  def compute(self, poly):
    force_sum = np.zeros((3, poly.nrVars()))
    for index in self.indexes:
      i = 3*index
      off = 0
      for j in range(len(poly.gravity_envelope)):
        force_sum[:, off+i:off+i+3] = np.eye(3)
        off += 3*len(poly.contacts)
    self.g = np.vstack([np.zeros((1, poly.nrVars())), force_sum])
    h_fc = np.zeros((self.size, 1))
    h_fc[0, 0] = self.limit*poly.mass*9.81
    self.h = h_fc

  def matrices(self):
    return self.g, self.h

class DistConstraint(object):

  """Constraint to limit position of the CoM w.r. to an origin.
  The origin will be clamped to dimension of the polygon."""

  def __init__(self, origin, radius):
    """Default constructor.

    :param origin: Origin of the circle/sphere.
    :param radius: Radius of the circle/sphere
    :type origin:  np.array(n, 1)

    """
    self.origin = origin
    self.radius = radius

    self.ctype = Constraint.Conic

  def compute(self, poly):
    self.size = poly.size_z()+1
    S = np.zeros((self.size, poly.nrVars()))
    S[1:, -poly.size_z():] = -np.eye(poly.size_z())
    r = np.zeros((self.size, 1))
    r[0] = self.radius
    r[1:] = -self.origin[:poly.size_z(), :]
    self.S = S
    self.r = r

  def matrices(self):
    return self.S, self.r

class CubeConstraint(object):

  """Constraint to limit position of the CoM to a cuboid centered at an origin"""

  def __init__(self, origin, length):
    """Default constructor. Origin will be clamped to the dimension of the polygon.

    :param origin: Origin of the cuboid.
    :param length: Length of the sides fo the box.
    :type origin: np.array(n, 1).
    """

    self.origin = origin
    self.length = length
    self.ctype = Constraint.Inequality

  def compute(self, poly):
    self.size = 2*poly.size_z()
    C = np.zeros((self.size, poly.nrVars()))
    C[:, -poly.size_z():] = np.vstack([np.eye(poly.size_z()), -np.eye(poly.size_z())])

    d = np.zeros((self.size, 1))
    d[:poly.size_z(), :] = self.origin[:poly.size_z(), :]+self.length
    d[poly.size_z():, :] = -self.origin[:poly.size_z(), :]+self.length

    self.C, self.d = C, d

  def matrices(self):
    return self.C, self.d
