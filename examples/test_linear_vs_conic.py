import stability as stab
import numpy as np
import cdd
import copy
from stability.linear_cone import build_cone
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

np.set_printoptions(linewidth=1000000, precision=2, threshold=1e12)

execfile('stairs_contacts.py')

def main():
  global pos, normals

  #bar = sum(pos)/float(len(pos))
  #pos = [p-bar for p in pos]

  mu = 0.1
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  poly = stab.StabilityPolygon(60, dimension=2,
                               force_lim=1000.)

  poly.contacts = contacts[0:8]

  for c in poly.contacts[0:4]:
    c.r[2] = 0.
    c.n = np.array([[0.4], [0.4], [np.sqrt(1-2*(0.4**2))]])
    #c.n = np.array([[0., 1., 0.]]).T

  for c in poly.contacts[4:]:
    c.r[2] = 0.
    #c.n = np.array([[0., 0., 1.]]).T
    c.n = np.array([[-0.4], [-0.4], [np.sqrt(1-2*(0.4**2))]])

  #poly.reset_fig()
  #poly.plot_contacts()
  #poly.show()

  #poly.make_problem()
  #poly.check_sizes()

  #point = np.array([[0.36510907, 0.31419711, 0.73607441]]).T

  #poly.addTorqueConstraint(contacts[-4:-2],
  #                         point,
  #                         10*np.ones((3, 1)))

  #poly.addTorqueConstraint(contacts[-2:],
  #                         point,
  #                         10*np.ones((3, 1)))

  #bar1 = sum([c.r for c in poly.contacts[0:4]])/4
  #bar2 = sum([c.r for c in poly.contacts[4:8]])/4
  #poly.addDistConstraint(bar1, 1.5)
  #poly.addDistConstraint(bar2, 1.5)

  #poly.make_problem()
  #poly.reset_fig()
  #poly.plot_contacts()
  #poly.ax.plot(point[0], point[1], point[2], 'o', markersize=10)
  #poly.show()

  sol = 'plain'

  poly.compute(stab.Mode.iteration, maxIter=50, epsilon=2e-3,
               solver=sol, plot_error=False, plot_step=False,
               plot_init=False, plot_final=False)

  assert(poly.dimension == 2)
  assert(len(poly.gravity_envelope) == 1)

  A1, A2, t = poly.A1, poly.A2, poly.t

  nr_generators = 4
  friction_cstr = []

  for i, contact in enumerate(poly.contacts):
    cone = build_cone(nr_generators, mu, contact.n, offset_angle=np.pi/4)
    nr_ineq = cone.shape[0]
    cstr = np.zeros((nr_ineq+1, 1+poly.nrVars()))
    cstr[:-1, 0] = cone[:, 0]
    cstr[:-1, 3*i+1:3*i+4] = cone[:, 1:]

    cstr[-1, 0] = 0.
    cstr[-1, 3*i+1:3*i+4] = contact.n.T
    friction_cstr.append(copy.deepcopy(cstr))

  com_max = 2.
  limit_cstr = np.zeros((4, 1+poly.nrVars()))
  limit_cstr[:, 0:1] = com_max*np.vstack((np.ones((2, 1)), np.ones((2, 1))))
  limit_cstr[:, -2:] = np.vstack((-np.identity(2), np.identity(2)))

  nr_forces = 3*len(poly.contacts)
  f_max = 5.*poly.mass*9.81
  limit_force = np.zeros((2*nr_forces, 1+poly.nrVars()))
  limit_force[:, 0:1] = f_max*np.vstack((np.ones((nr_forces, 1)), np.ones((nr_forces, 1))))
  limit_force[:, 1:-2] = np.vstack((-np.identity(nr_forces), np.identity(nr_forces)))

  friction = np.vstack(friction_cstr)
  eq = np.hstack((t, -A1, -A2))

  mat = cdd.Matrix(friction, number_type='fraction')
  #mat = cdd.Matrix(friction)
  mat.rep_type = cdd.RepType.INEQUALITY
  mat.extend(limit_force)
  mat.extend(limit_cstr)
  print np.vstack(friction)
  mat.extend(eq, linear=True)
  cdd_poly = cdd.Polyhedron(mat)
  vertices = np.array(cdd_poly.get_generators())

  print vertices.shape
  if len(vertices.shape) > 1:
    points = vertices[:, -2:]
    hull = ConvexHull(points)

    poly.reset_fig()
    poly.plot_contacts()
    poly.plot_polyhedron(poly.inner, 'blue', 0.5)
    poly.ax.plot(points[hull.vertices, 0], points[hull.vertices, 1], label='cdd')
    poly.show()
  else:
    print "No vertices"

  #if sol == 'plain':
  #  ineq = [l/abs(l[0]) for l in poly.inner.inequalities]
  #  print np.vstack(ineq)
  #elif sol == 'cdd':
  #  poly = cdd.Polyhedron(poly.inner)
  #  print poly.get_inequalities()
  #return poly

if __name__ == '__main__':
  main()
