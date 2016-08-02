#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
