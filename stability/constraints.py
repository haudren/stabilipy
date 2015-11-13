#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from enum import Enum, unique
import numpy as np

def cross_m(vec):
  return np.array([[0, -vec.item(2), vec.item(1)],
                   [vec.item(2), 0, -vec.item(0)],
                   [-vec.item(1), vec.item(0), 0]])

@unique
class Constraint(Enum):
  """Constraint types. Can only be inequality or conic."""
  Inequality = 1
  Conic = 2

class TorqueConstraint(object):

  """Constraint to limit torque applied on certain contacts"""

  def __init__(self, indexes, point, ub, lb=None):
    """Default constructor.

    :indexes: Indexes of the contacts on which the constraint applies
    :point: Point where the torques are computed (3,1) array
    :ub: Upper bound (3,1) array
    :lb: Lower bound, can be None (3,1) array

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

    :indexes: Indexes of the contacts on which the constraint applies
    :limit: Maximum force, expressed as percentage of robot weight

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

  """Constraint to limit position of the CoM w.r. to an origin"""

  def __init__(self, origin, radius):
    """Default constructor.

    :origin: Origin of the circle/sphere. Should be (n, 1)
    Will be clamped to first size_z() dimensions.
    :radius: Radius of the circle/sphere

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
    r[1:] = self.origin[:poly.size_z(), :]
    self.S = S
    self.r = r

  def matrices(self):
    return self.S, self.r
