import numpy as np

pos = []
normals = []

p = [[-0.4722227, -0.24517583, -0.6370031]]
n = [[2.02215104e-04, -3.23903880e-05, 9.99999979e-01]]
pos.append(p)
normals.append(n)

p = [[-0.2549828, -0.24587737, -0.63704705]]
n = [[2.02215104e-04, -3.23903880e-05, 9.99999979e-01]]
pos.append(p)
normals.append(n)

p = [[-0.25787751, -0.38255749, -0.63705089]]
n = [[2.02215104e-04, -3.23903880e-05, 9.99999979e-01]]
pos.append(p)
normals.append(n)

p = [[-0.47206733, -0.38317576, -0.6370076]]
n = [[2.02215104e-04, -3.23903880e-05, 9.99999979e-01]]
pos.append(p)
normals.append(n)

#Contact lgripper/handrail
#Left
p = [[0.3651077, 0.33419711, 0.63609439]]
n = [[-3.39491173e-05, 9.99999875e-01, 4.99472000e-04]]
pos.append(p)
normals.append(n)

#Right
#p = [[0.36510907, 0.29419711, 0.63607441]]
#p = [[0.3651077, 0.33419711, 0.63609439]]
#n = [[3.44761855e-05, -9.99999874e-01, -5.00077386e-04]]
#pos.append(p)
#normals.append(n)

#Bottom
#p = [[0.34212609, 0.31418314, 0.66248165]]
#n = [[-6.56636734e-01, -3.99160434e-04, 7.54206895e-01]]
#pos.append(p)
#normals.append(n)
#
##Top
p = [[0.38480749, 0.31420908, 0.61345819]]
n = [[6.56636734e-01, 4.00439950e-04, -7.54206894e-01]]
pos.append(p)
normals.append(n)


pos = [np.array(px).T for px in pos]
#for p in pos:
#  p[2, 0] = 0.0
normals = [np.array(nx).T for nx in normals]
