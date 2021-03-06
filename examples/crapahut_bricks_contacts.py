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

points = []
normals = []
points.append([0.550119, 0.0376734, 0.0769174])
normals.append([-0.281987, 0.00328666, 0.959413])
points.append([0.546072, 0.174314, 0.0752598])
normals.append([-0.281987, 0.00328666, 0.959413])
points.append([0.340578, 0.172651, 0.0148676])
normals.append([-0.281987, 0.00328666, 0.959413])
points.append([0.341712, 0.0346581, 0.0156734])
normals.append([-0.281987, 0.00328666, 0.959413])
points.append([-0.0724719, -0.0343614, -0.00169742])
normals.append([-0.0018678, 7.97359e-06, 0.999998])
points.append([-0.0723575, -0.172361, -0.00169611])
normals.append([-0.0018678, 7.97359e-06, 0.999998])
points.append([0.141832, -0.171807, -0.00129605])
normals.append([-0.0018678, 7.97359e-06, 0.999998])
points.append([0.144767, -0.0351274, -0.00129166])
normals.append([-0.0018678, 7.97359e-06, 0.999998])
points.append([0.449843, -0.363954, 0.758192])
normals.append([-0.0157233, -0.0119195, 0.999805])
points.append([0.30986, -0.363778, 0.755992])
normals.append([-0.0157233, -0.0119195, 0.999805])
points.append([0.309774, -0.443773, 0.755037])
normals.append([-0.0157233, -0.0119195, 0.999805])
points.append([0.449757, -0.443949, 0.757237])
normals.append([-0.0157233, -0.0119195, 0.999805])
points = [np.array([pi]).T for pi in points]
normals = [np.array([ni]).T for ni in normals]
