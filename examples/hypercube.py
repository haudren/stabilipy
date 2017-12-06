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

import stabilipy as stab
import numpy as np

if __name__ == '__main__':

  A = np.vstack((np.eye(6), -np.eye(6)))
  b = np.ones(12,)

  linear_proj = stab.LinearProjection(3, A, b, None, None)
  linear_proj.compute(stab.Mode.precision, solver='cdd', epsilon=1e-3)
