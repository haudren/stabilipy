#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy.spatial as spatial

def other_points(i, simplex):
  k, l, m = simplex
  if i == k:
    return l, m
  elif i == l:
    return k, m
  elif i == m:
    return k, l
  else:
    raise ValueError("The point {} is not in simplex {}".format(i, simplex))

def cap_hull(ch, n, d):
  above = ch.points.dot(n) + d < 0
  above_i = np.where(above)[0]
  points = [ch.points[~above.squeeze(), :]]

  for i in above_i:
    for simplex in ch.simplices:
      if i in simplex:
        j, k = other_points(i, simplex)
        if not j in above_i:
          alpha = -(d + n.squeeze().dot(ch.points[i, :]))/(n.T.dot(ch.points[j, :] - ch.points[i, :]))
          points.append(ch.points[i]*(1-alpha) + ch.points[j]*alpha)
        if not k in above_i:
          alpha = -(d + n.squeeze().dot(ch.points[i, :]))/(n.T.dot(ch.points[k, :] - ch.points[i, :]))
          points.append(ch.points[i]*(1-alpha) + ch.points[k]*alpha)

  hull = spatial.ConvexHull(np.vstack(points))
  return hull
