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
#p = [[0, 0, 1]]
#n = [[0, 0, 1]]
#pos.append(p)
#normals.append(n)

pos = [np.array(px).T for px in pos]
#for p in pos:
#  p[2, 0] = 0.0
normals = [np.array(nx).T for nx in normals]
