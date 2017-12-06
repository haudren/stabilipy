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

def contacts():
  p = [np.array([[0.21699991],
                 [0.02499986],
                 [0.23000027]]),
       np.array([[0.08199991],
                 [0.02500003],
                 [0.23000024]]),
       np.array([[0.08200009],
                 [0.16500003],
                 [0.22999976]]),
       np.array([[0.21700009],
                 [0.16499986],
                 [0.22999979]])]

  n = [np.array([[-2.05321496e-07],
                 [3.43245159e-06],
                 [1.00000000e+00]]),
       np.array([[-2.05321496e-07],
                 [3.43245159e-06],
                 [1.00000000e+00]]),
       np.array([[-2.05321496e-07],
                 [3.43245159e-06],
                 [1.00000000e+00]]),
       np.array([[-2.05321496e-07],
                 [3.43245159e-06],
                 [1.00000000e+00]])]
  return p, n
