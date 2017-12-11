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
from builtins import map
from builtins import zip
import numpy as np
import stabilipy as stab


#normals = map(stab.normalize,
#              [np.array([[0, 0, 1]]).T,
#               np.array([[0, 0, 1]]).T,
#               np.array([[1, 0, 0]]).T])
#

def main():
  #normals = map(stab.normalize,
  #              [np.array([[-0.7, 0, 1]]).T,
  #               np.array([[0.1, 0.5, 1]]).T,
  #               np.array([[0.5, -0.5, 1]]).T])

  #pos = [np.array([[0, 0, 0]]).T,
  #       np.array([[0, 1, 0]]).T,
  #       np.array([[1, 0, 0.2]]).T]

  #normals = map(stab.normalize,
  #              [np.array([[0., 0., 1.]]).T,
  #               np.array([[0., 0., 1.]]).T,
  #               np.array([[0., 0., 1.]]).T,
  #               np.array([[0., 0., 1.]]).T])

  #pos = [np.array([[0.1, 0.155, -0.105]]).T,
  #       np.array([[-0.07, 0.155, -0.105]]).T,
  #       np.array([[-0.07, 0.055, -0.105]]).T,
  #       np.array([[0.1, 0.055, -0.105]]).T]

  #normals = [np.array([[0, 0, 1]]).T]*3

  normals = list(map(stab.normalize,
                [np.array([[-1, 0.0, 0]]).T,
                 np.array([[1, 0.0, 0]]).T,
                 np.array([[0.0, 1., 0]]).T]))

  pos = [np.array([[1., 0., 0.]]).T,
         np.array([[-1., 0.0, 0.]]).T,
         np.array([[0., -1., 0.]]).T]

  bar = sum(pos)/float(len(pos))
  pos = [p-bar for p in pos]

  mu = 0.7
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  poly = stab.StabilityPolygon(60)
  poly.contacts = contacts

  poly.compute(3e-2, True, False, False, True)

main()
