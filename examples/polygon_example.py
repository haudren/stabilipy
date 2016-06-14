import stability as stab
import numpy as np
import sys

from polygon_contacts import contacts

from scipy.spatial import ConvexHull

import matplotlib.pyplot as plt

azim = 48.9035087719
elev = 31.5350877193
xlim = [-0.95389899, 0.95389899]
ylim = [-0.95389899, 0.95389899]
zlim = [-0.95389899, 0.95389899]

def main(n, alpha):
  if alpha < 1/float(n):
    raise RuntimeError("Unfeasible")

  mu = 0.5
  pos, normals = contacts(n)
  cont = [stab.Contact(mu, pi, ni) for pi, ni in zip(pos, normals)]

  polygon = stab.StabilityPolygon(200, dimension=2, radius=1.5)
  polygon.contacts = cont

  for i in range(n):
    polygon.addForceConstraint([polygon.contacts[i]], alpha)

  polygon.compute(stab.Mode.best, epsilon=2e-3, maxIter=100, solver='plain',
                  record_anim=False, plot_init=False,
                  plot_step=False, plot_final=False)
  polygon.reset_fig()
  polygon.plot_contacts()

  points = polygon.inner.vertices
  hull = ConvexHull(points)

  polygon.ax.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)

  for i in range(n):
    p0 = polygon.contacts[i].r[0:2, :]
    p1 = polygon.contacts[i-1].r[0:2, :]
    if i < n-1:
      p2 = polygon.contacts[i+1].r[0:2, :]
    else:
      p2 = polygon.contacts[i-(n-1)].r[0:2, :]

    p3 = polygon.contacts[i-2].r[0:2, :]

    if i < n-2:
      p4 = polygon.contacts[i+2].r[0:2, :]
    else:
      p4 = polygon.contacts[i-(n-2)].r[0:2, :]

    p5 = polygon.contacts[i-3].r[0:2, :]

    if i < n-3:
      p6 = polygon.contacts[i+3].r[0:2, :]
    else:
      p6 = polygon.contacts[i-(n-3)].r[0:2, :]


    if alpha > 0.5: #1/2
      r1 = alpha*p0 + (1-alpha)*p1
      r2 = alpha*p0 + (1-alpha)*p2
    elif alpha > 0.33: # 1/3
      r1 = alpha*p0 + alpha*p1 + (1-2*alpha)*p2
      r2 = alpha*p0 + alpha*p2 + (1-2*alpha)*p1
    elif alpha > 0.25: # 1/4
      r1 = alpha*p0 + alpha*p1 + alpha*p2 + (1 - 3*alpha)*p3
      r2 = alpha*p0 + alpha*p2 + alpha*p1 + (1 - 3*alpha)*p4
    elif alpha > 0.2: # 1/5
      r1 = alpha*p0 + alpha*p1 + alpha*p2 + alpha*p3 + (1-4*alpha)*p4
      r2 = alpha*p0 + alpha*p2 + alpha*p1 + alpha*p4 + (1-4*alpha)*p3
    elif alpha > 0.16: #1/6
      r1 = alpha*p0 + alpha*p1 + alpha*p2 + alpha*p3 + alpha*p4 + (1-5*alpha)*p6
      r2 = alpha*p0 + alpha*p2 + alpha*p1 + alpha*p4 + alpha*p3 + (1-5*alpha)*p5
    else:
      raise RuntimeError("Bound too low")

    polygon.ax.plot(r1[0], r1[1], marker='^', markersize=20)
    polygon.ax.plot(r2[0], r2[1], marker='^', markersize=20)

  polygon.show()

print("n-sided gon : {} with limit {}".format(int(sys.argv[1]), float(sys.argv[2])))

main(int(sys.argv[1]), float(sys.argv[2]))
