from __future__ import division
import stability as stab
import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import cdd

from pyparma import Polyhedron as PPL_Poly
from fractions import Fraction

import shapes

np.set_printoptions(linewidth=1000)

execfile('stairs_contacts.py')

def main():
  global pos, normals

  #bar = sum(pos)/float(len(pos))
  #pos = [p-bar for p in pos]

  mu = 0.7
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  poly = stab.StabilityPolygon(100, dimension=3)
  poly.contacts = contacts

  #poly.reset_fig()
  #poly.plot_contacts()
  #poly.show()

  #poly.make_problem()
  #poly.check_sizes()


  #poly.contacts = contacts[0:4]

  #poly.compute(stab.Mode.best, maxIter=50, epsilon=1e-2, solver='cdd',
  #             plot_init=False, plot_final=False, plot_step=False,
  #             plot_direction=False, plot_error=False)

  #p1 = poly.polygon()

  #poly.contacts = contacts

  point = np.array([[0.36510907, 0.31419711, 0.73607441]]).T

  poly.addTorqueConstraint(contacts[-4:-2],
                           point,
                           10*np.ones((3, 1)))

  gravities = np.linspace(0, 2, 10)
  shape = [
              np.array([[-1., 0, 0]]).T,
              np.array([[1., 0, 0]]).T,
              np.array([[0, 1., 0]]).T,
              np.array([[0, -1., 0]]).T
          ]

  for gravity in gravities:
    poly.gravity_envelope = [gravity*s for s in shape]
    poly.compute(stab.Mode.best, maxIter=50, epsilon=1e-2, solver='cdd',
                 plot_init=False, plot_final=True, plot_step=False,
                 plot_direction=False, plot_error=False)

  #p0 = poly.polygon()
  #interp = shapes.PolygonInterpolator(p0, p1)

  #poly.addTorqueConstraint(contacts[-2:],
  #                         point,
  #                         10*np.ones((3, 1)))
  #poly.reset_fig()
  #poly.plot_contacts()

  #f = open('top_lel.txt', 'w')
  #for i in range(101):
  #  poly.forces_constraints = []
  #  for k in range(4, 8):
  #    poly.addForceConstraint(contacts[k:k+1], (100-i)*0.01)

  #  poly.compute(stab.Mode.best, maxIter=50, epsilon=1e-2, solver='cdd',
  #               plot_init=False, plot_final=False, plot_step=False,
  #               plot_direction=False, plot_error=False)

  #  cur_p = poly.polygon()
  #  cur_area = cur_p.area
  #  max_j = 1.
  #  min_j = 0.
  #  j = 0.5

  #  pi = interp.fast_interpolate(j)
  #  interp_area = pi.area

  #  nr_iter = 0
  #  while abs(cur_area - interp_area) > 1e-3 and nr_iter < 200:
  #    if cur_area > interp_area:
  #      max_j = j
  #      j = (min_j+j)/2
  #    else:
  #      min_j = j
  #      j = (max_j+j)/2
  #    pi = interp.fast_interpolate(j)
  #    interp_area = pi.area
  #    #print cur_area - interp_area, max_j, j, min_j
  #    nr_iter += 1
  #  f.write("{},{}\n".format((100-i)*0.01, j))
  #  #poly.plot()
  #  #poly.ax.plot(*pi.exterior.coords.xy)
  #  #poly.show()
  #f.close()

  ##A = np.vstack([p.T for p in poly.points])
  #print zip(*A)
  #qhull = ConvexHull(A)
  #coords, tri = [c for c in qhull.points.T], qhull.simplices

  #dirs = np.vstack(poly.directions[:-1])
  #offsets = np.vstack(poly.offsets[:-1])
  #mat = np.hstack([offsets, -dirs])
  #out = cdd.Matrix(mat, number_type='fraction')
  #print np.array(out)
  #out.rep_type = cdd.RepType.INEQUALITY
  #polyhedron = cdd.Polyhedron(out)
  #gen = np.array(polyhedron.get_generators())

  #f = np.vectorize(lambda x: Fraction(str(x)).limit_denominator(1000))
  #ppl_poly = PPL_Poly(hrep=f(mat))
  #p3 = ppl_poly.vrep()[:, 1:]
  #qhull = ConvexHull(p3)
  #coords3, tri3 = [c for c in qhull.points.T], qhull.simplices

  #print np.array(out).shape
  #print np.array(poly.outer).shape
  ##print np.array(out) - np.vstack([np.array(poly.outer), np.zeros((1, 4))])

  #p2 = gen[:, 1:]
  #qhull = ConvexHull(p2)
  #coords2, tri2 = [c for c in qhull.points.T], qhull.simplices

  #plt.close("all")
  #fig = plt.figure()
  #ax = fig.add_subplot('111', aspect='equal', projection='3d')
  #x, y, z = zip(*A)
  #ax.scatter(x, y, z)
  #x, y, z = zip(*p3)
  #ax.scatter(x, y, z, color='b', marker='^', s=120.)

  #ax.plot_trisurf(*coords, triangles=tri, color='r', alpha=0.5)

  #ax.plot_trisurf(*coords2, triangles=tri2, color='b', alpha=0.1)

  #ax.plot_trisurf(*coords3, triangles=tri3, color='g', alpha=0.1)

  #plt.show()

  return poly

if __name__ == '__main__':
  main()
