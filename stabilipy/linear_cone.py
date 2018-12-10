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
import numpy as np
import cdd

from .utils import cross_m

def build_cone(nr_generators, mu, axis, offset_angle=0):
  origin = np.array([[0, 0, 0]])

  angles = [offset_angle+2*np.pi*i/nr_generators for i in range(nr_generators)]
  local_points = [np.array([[mu*np.cos(x), mu*np.sin(x), 1.]]).T for x in angles]
  zaxis = np.array([[0., 0., 1.]]).T
  R = rotate_axis(zaxis, axis)

  points = np.vstack([(R.dot(p)).T for p in local_points])

  all_points = np.vstack((origin, points))

  vertextypes = np.vstack([np.ones((1, 1)), np.zeros((points.shape[0], 1))])

  cdd_points = np.hstack((vertextypes, all_points))

  mat = cdd.Matrix(cdd_points)

  mat.rep_type = cdd.RepType.GENERATOR

  poly = cdd.Polyhedron(mat)
  ineqs = np.array(poly.get_inequalities())[:-1, :]
  return ineqs

def rotate_axis(a, b):
  """Compute rotation as a 3x3 mattrix that rotates unit vector a onto b"""
  assert(abs(np.linalg.norm(a) - 1.0) < 1e-5)
  assert(abs(np.linalg.norm(b) - 1.0) < 1e-5)

  v = np.cross(a, b, axis=0)
  s = np.linalg.norm(v)
  c = np.dot(a.T, b)
  cv = cross_m(v)

  if s != 0:
    R = np.eye(3) + cv + np.dot(cv, cv)*(1-c)/s
  else:
    R = np.eye(3)
  return R
