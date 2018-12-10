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
import numpy as np

normalize = lambda x: x/np.linalg.norm(x)

pos, normals = [], []
#Right foot
p = [[-0.1, 0.4, 0.]]
n = [[0., 0.3, 1.]]

pos.append(p)
normals.append(n)

#Left foot
p = [[0., -0.4, 0.1]]
n = [[0., -0.3, 1.]]

pos.append(p)
normals.append(n)

#Hand
p = [[0.5, 0., 0.8]]
n = [[1., 0., 1.]]

pos.append(p)
normals.append(n)

pos = [np.array(pi).T for pi in pos]
normals = [np.array(normalize(ni)).T for ni in normals]
