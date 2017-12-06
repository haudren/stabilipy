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

import cdd
import numpy as np

radius = 0.1
points = radius*np.vstack([np.eye(3), -np.eye(3)])
mat_p = cdd.Matrix(np.hstack([np.ones((6, 1)), points]))
mat_p.rep_type = cdd.RepType.GENERATOR

sphere_ineq = np.array(cdd.Polyhedron(mat_p).get_inequalities())
print sphere_ineq
