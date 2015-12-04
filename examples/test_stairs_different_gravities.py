from __future__ import division
import stability as stab
import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(linewidth=1000)

from stairs_contacts import pos, normals

def main():
  global pos, normals

  mu = 0.7
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  poly_3d = stab.StabilityPolygon(100, dimension=3)
  poly_3d.contacts = contacts

  poly_2d = stab.StabilityPolygon(100, dimension=2)
  poly_2d.contacts = contacts

  polys = [poly_2d, poly_3d]

  point = np.array([[0.36510907, 0.31419711, 0.73607441]]).T

  for poly in polys:
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
    for poly in polys:
      poly.gravity_envelope = [gravity*s for s in shape]
      poly.compute(stab.Mode.best, maxIter=50, epsilon=1e-2, solver='cdd',
                   plot_init=False, plot_final=False, plot_step=False,
                   plot_direction=False, plot_error=False)
      poly.plot()
    plt.show()
  return poly

if __name__ == '__main__':
  main()
