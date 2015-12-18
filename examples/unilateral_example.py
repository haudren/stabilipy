import stability as stab
import numpy as np
import sys

from unilateral_contacts import pos, normals

import matplotlib.pyplot as plt

azim = 48.9035087719
elev = 31.5350877193
xlim = [-0.95389899, 0.95389899]
ylim = [-0.95389899, 0.95389899]
zlim = [-0.95389899, 0.95389899]

def main(margin):
  mu = 0.5
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  contacts[2].mu = 0.5

  polyhedron = stab.StabilityPolygon(200, dimension=3, radius=1.5)
  polyhedron.contacts = contacts

  polygon = stab.StabilityPolygon(200, dimension=2, radius=1.5)
  polygon.contacts = contacts

  shape = [
              np.array([[-1., 0, 0]]).T,
              np.array([[1., 0, 0]]).T,
              np.array([[0, 1., 0]]).T,
              np.array([[0, -1., 0]]).T
          ]

  polytope = [margin*s for s in shape]

  polyhedron.gravity_envelope = polytope
  polyhedron.compute(stab.Mode.best, epsilon=2e-3, maxIter=50, solver='parma',
                     record_anim=False, plot_init=False,
                     plot_step=False, plot_final=False)
  polyhedron.reset_fig()
  polyhedron.ax.set_xlabel("x(m)")
  polyhedron.ax.set_ylabel("y(m)")
  polyhedron.ax.set_zlabel("z(m)")
  polyhedron.ax.view_init(elev=elev, azim=azim)
  polyhedron.ax.set_xlim3d(*xlim)
  polyhedron.ax.set_ylim3d(*ylim)
  polyhedron.ax.set_zlim3d(*zlim)
  polyhedron.plot_contacts()
  polyhedron.plot_solution()
  polyhedron.plot_polyhedrons()
  #polyhedron.show()

  plt.savefig('{}.png'.format(margin))

  polygon.gravity_envelope = polytope
  polygon.compute(stab.Mode.best, epsilon=2e-3, maxIter=20, solver='parma',
                  record_anim=False, plot_init=False,
                  plot_step=False, plot_final=False)

  print polyhedron.volume_convex(polyhedron.inner)
  print 3*polygon.volume_convex(polygon.inner)

print "Margin : {}".format(sys.argv[1])

main(float(sys.argv[1]))
