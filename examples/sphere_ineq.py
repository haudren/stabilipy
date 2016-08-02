import cdd
import numpy as np

radius = 0.1
points = radius*np.vstack([np.eye(3), -np.eye(3)])
mat_p = cdd.Matrix(np.hstack([np.ones((6, 1)), points]))
mat_p.rep_type = cdd.RepType.GENERATOR

sphere_ineq = np.array(cdd.Polyhedron(mat_p).get_inequalities())
print sphere_ineq
