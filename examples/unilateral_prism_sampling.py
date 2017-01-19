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

  shape = [
              np.array([[-1., 0, 0]]).T,
              np.array([[1., 0, 0]]).T,
              np.array([[0, 1., 0]]).T,
              np.array([[0, -1., 0]]).T,
              #np.array([[1, -1., 0]]).T/np.sqrt(2),
          ]

  polytope = [margin*s for s in shape]

  prisms = stab.PrismIntersection(200, polytope, contacts, radius=1.5)

  points = np.array([1., 1., 3.])*(np.random.random((10**6,3)) - 0.5)
  np.savetxt('random_points.txt', points)
  #points = np.loadtxt('random_points.txt')

  truths, iters = zip(*[prisms.sample(point)
                        for point in points])
  #truths = np.array(truths)

  #prisms.plot()

  print("total iters: {}".format(np.sum(iters)))

  #prisms.threedax.plot(*zip(*points[truths, :]), linestyle="none", marker="*", markerfacecolor="green")
  #prisms.threedax.plot(*zip(*points[~truths, :]), linestyle="none", marker="x", markerfacecolor="red")

  #prisms.show()

print "Margin : {}".format(sys.argv[1])

main(float(sys.argv[1]))
