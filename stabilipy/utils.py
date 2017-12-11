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
from __future__ import absolute_import
from builtins import zip
import numpy as np
from .backends import CDDBackend, PlainBackend, ParmaBackend

try:
  import shapely.geometry as geom
  HAS_SHAPELY=True
except ImportError:
  HAS_SHAPELY=False


def cross_m(vec):
  return np.array([[0, -vec.item(2), vec.item(1)],
                   [vec.item(2), 0, -vec.item(0)],
                   [-vec.item(1), vec.item(0), 0]])


def normalize(vec):
  return vec/(np.linalg.norm(vec))


def polygon(self, centroid=None):
  """Return the inner approximation as a shapely polygon. Only valid in 2D"""

  if not HAS_SHAPELY:
    raise RuntimeError("You need to have shapely installed to use this method")

  assert(self.dimension == 2)
  if isinstance(self.backend, CDDBackend):
    gen = np.array(self.inner)[:, 1:]
  elif isinstance(self.backend, PlainBackend):
    gen = self.inner.vertices
  elif isinstance(self.backend, ParmaBackend):
    gen = self.inner.vrep()[:, 1:]
  else:
    raise NotImplemented("No polygon method defined for this backend")

  p = np.vstack([gen, gen[0, :]])
  if centroid is None:
    x, y = p[:, 0], p[:, 1]
  else:
    x, y = p[:, 0] + centroid.item(0), p[:, 1] + centroid.item(1)
  return geom.Polygon(list(zip(x, y))).convex_hull
