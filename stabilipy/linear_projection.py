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
from .recursive_projection import RecursiveProjectionProblem
import scipy.optimize as opt
import numpy as np

class LinearProjection(RecursiveProjectionProblem):

  """Recursively compute the projection of a linear set:
    :math:`A x \leq b, C x = d` onto :math:`y = E x + f`"""

  def __init__(self, dimension, A, b, C, d, E=None, f=None):
    """Create a linear projection problem.

       :param dimension: Dimension on which to project
       :param A: Linear inequality matrix
       :param b: Linear inequality RHS
       :param C: Linear equality matrix
       :param d: Linear equality RHS
       :param E: Projection matrix
       :param f: Projection offset

       :type dimension: int
       :type A: np.array(nrineq, dim)
       :type b: np.array(nrineq,)
       :type C: np.array(nreq, dim)
       :type d: np.array(nreq,)
       :type E: np.array(dimension, dim)
       :type f: np.array(dimension,)
       """

    RecursiveProjectionProblem.__init__(self, dimension=dimension)
    self.A = A
    self.b = b
    self.C = C
    self.d = d

    if E is None:
      self.E = np.zeros((dimension, A.shape[1]))
      self.E[:, :dimension] = np.eye(dimension)
    else:
      self.E = E

    if f is None:
      self.f = np.zeros((dimension,))
    else:
      self.f = f

  def solve(self, d):
    #Careful, d is direction but self.d is equality RHS
    res = opt.linprog(-self.E.T.dot(d).squeeze(), self.A, self.b, self.C, self.d, bounds=((None, None)))
    if res.success:
      return (self.E.dot(res.x) + self.f).reshape((self.dimension, 1))
    else:
      return None
