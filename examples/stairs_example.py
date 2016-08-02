import stability as stab
import numpy as np
import cdd
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

execfile('stairs_contacts.py')

def main():
  global pos, normals

  #bar = sum(pos)/float(len(pos))
  #pos = [p-bar for p in pos]

  offset = np.array([[0.0], [0.], [0.]])
  mu = 0.5
  contacts = [stab.Contact(mu, offset+p, n) for p, n in zip(pos, normals)]

  poly = stab.StabilityPolygon(60, dimension=2, radius=100.)
  poly.printer.verbosity = stab.Verbosity.none
  poly.contacts = contacts

  point = offset+np.array([[0.36510907, 0.31419711, 0.73607441]]).T

  poly.make_problem()
  poly.reset_fig()
  poly.plot_contacts()
  poly.ax.plot(point[0], point[1], point[2], 'o', markersize=10)
  poly.show()

  sol = 'cdd'

  poly.compute(stab.Mode.iteration, maxIter=50, epsilon=2e-3,
               solver=sol, plot_error=False, plot_step=False,
               plot_init=False, plot_final=False)

  unconstrained = poly.polyhedron()

  bar1 = sum([c.r for c in poly.contacts[0:4]])/4
  bar2 = sum([c.r for c in poly.contacts[4:8]])/4

  #Foot
  poly.addDistConstraint(bar1, 1.5)
  #Hand
  poly.addDistConstraint(bar2, 1.5)

  poly.compute(stab.Mode.iteration, maxIter=50, epsilon=2e-3,
               solver=sol, plot_error=False, plot_step=False,
               plot_init=False, plot_final=False)

  dist_constrained = poly.polyhedron()

  poly.addTorqueConstraint(contacts[-4:-2],
                           point,
                           1*np.ones((3, 1)))

  poly.addTorqueConstraint(contacts[-2:],
                           point,
                           1*np.ones((3, 1)))

  poly.compute(stab.Mode.iteration, maxIter=50, epsilon=2e-3,
               solver=sol, plot_error=False, plot_step=False,
               plot_init=False, plot_final=False)

  all_constrained = poly.polyhedron()

  def envelope(points):
    hull = ConvexHull(points)
    exterior = hull.points[hull.vertices, :]
    exterior = np.vstack([exterior, hull.points[hull.vertices[0], :]])
    return exterior

  poly.reset_fig()
  poly.plot_contacts()

  unconstrained = envelope(unconstrained)
  poly.ax.plot(unconstrained[:, 0], unconstrained[:, 1], label="Unconstrained")
  dist_constrained = envelope(dist_constrained)
  poly.ax.plot(dist_constrained[:, 0], dist_constrained[:, 1], label="Dist constraint")
  all_constrained = envelope(all_constrained)
  poly.ax.plot(all_constrained[:, 0], all_constrained[:, 1], label="Dist + force constrained")
  plt.legend()
  plt.show()

  if sol == 'plain':
    ineq = [l/abs(l[0]) for l in poly.inner.inequalities]
    print np.vstack(ineq)
  elif sol == 'cdd':
    poly = cdd.Polyhedron(poly.inner)
    print poly.get_inequalities()
  return poly

if __name__ == '__main__':
  main()
