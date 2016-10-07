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

  shape = [
              np.array([[-1., 0, 0]]).T,
              np.array([[1., 0, 0]]).T,
              np.array([[0, 1., 0]]).T,
              np.array([[0, -1., 0]]).T,
              np.array([[1, -1., 0]]).T/np.sqrt(2),
          ]

  polytope = [margin*s for s in shape]

  prisms = stab.PrismIntersection(200, polytope, contacts, radius=1.5)

  polyhedron.gravity_envelope = polytope
  polyhedron.compute(stab.Mode.best, epsilon=1e-3, maxIter=10, solver='qhull', plot_final=False)

  points = 3*(np.random.random((10000,3)) - 0.5)
  truths = np.array([prisms.sample(point)
                      for point in points])


  prisms.plot()
  #polyhedron.plot()

  prisms.threedax.plot(*zip(*points[truths, :]), linestyle="none", marker="*", markerfacecolor="green")
  prisms.threedax.plot(*zip(*points[~truths, :]), linestyle="none", marker="x", markerfacecolor="red")

  prisms.show()

  #polyhedron.show()

print "Margin : {}".format(sys.argv[1])

main(float(sys.argv[1]))
