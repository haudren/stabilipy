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

def main():
  mu = 0.5
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  contacts[2].mu = 0.5

  polyhedron = stab.StabilityPolygon(200, dimension=2, radius=1.5)
  polyhedron.contacts = contacts

  polyhedron.select_solver('plain')
  polyhedron.make_problem()
  polyhedron.init_algo()
  polyhedron.build_polys()

  points = np.random.random((100,2)) - 0.5

  truths = np.array([polyhedron.sample(point, plot_final=False, plot_step=False)
            for point in points])

  polyhedron.plot()

  polyhedron.ax.plot(*zip(*points[truths, :]), linestyle="none", marker="*", markerfacecolor="green")
  polyhedron.ax.plot(*zip(*points[~truths, :]), linestyle="none", marker="x", markerfacecolor="red")

  polyhedron.show()

if __name__ == "__main__":
  main()
