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

import numpy as np

pos = []
normals = []

p = [[-0.4722227, -0.24517583, -0.6370031]]
n = [[2.02215104e-04, -3.23903880e-05, 9.99999979e-01]]
#Left
p = [[0, 1, 0]]
n = [[0, 1, 0]]
pos.append(p)
normals.append(n)

#Right
p = [[0, -1, 0]]
n = [[0, -1, 0]]
pos.append(p)
normals.append(n)

#Front
p = [[1, 0, 0]]
n = [[1, 0, 0]]
pos.append(p)
normals.append(n)

#Back
p = [[-1, 0, 0]]
n = [[-1, 0, 0]]
pos.append(p)
normals.append(n)

#Top
p = [[0, 0, 1]]
n = [[0, 0, 1]]
pos.append(p)
normals.append(n)

pos = [np.array(px).T for px in pos]
#for p in pos:
#  p[2, 0] = 0.0
normals = [np.array(nx).T for nx in normals]
