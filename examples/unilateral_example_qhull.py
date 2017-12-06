import stabilipy as stab
import numpy as np
import sys

#from unilateral_contacts import pos, normals
from box_contacts import pos, normals

import matplotlib.pyplot as plt

azim = 48.9035087719
elev = 31.5350877193
xlim = [-0.95389899, 0.95389899]
ylim = [-0.95389899, 0.95389899]
zlim = [-0.95389899, 0.95389899]

def main(margin, solver):
  mu = 0.7
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  contacts[2].mu = 0.5

  polyhedron = stab.StabilityPolygon(200, dimension=3, radius=1.)
  polyhedron.contacts = contacts

  shape = [
              np.array([[-1., 0, 0]]).T,
              np.array([[1., 0, 0]]).T,
              np.array([[0, 1., 0]]).T,
              np.array([[0, -1., 0]]).T,
              #np.array([[0, 0., 1]]).T,
              #np.array([[0, 0., -1]]).T
          ]

  polytope = [margin*s for s in shape]

  polyhedron.gravity_envelope = polytope
  polyhedron.compute(stab.Mode.iteration, epsilon=2e-3, maxIter=1000, solver=solver,
                     record_anim=False, plot_init=False,
                     plot_step=False, plot_final=False, plot_direction=False,
                     plot_error=True)

print("Margin : {}".format(sys.argv[1]))

if len(sys.argv) >= 3:
  solver = sys.argv[2]
else:
  solver = 'qhull'

main(float(sys.argv[1]), solver)
