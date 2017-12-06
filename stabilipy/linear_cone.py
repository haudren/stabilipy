#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import cdd

from utils import cross_m

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
