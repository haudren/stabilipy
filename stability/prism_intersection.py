import stability
import numpy as np
from scipy.spatial import ConvexHull, HalfspaceIntersection

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def _make_polygon(mass, gravity, contacts, radius):
  polygon = stability.StabilityPolygon(mass, dimension=2, radius=radius)
  polygon.contacts = contacts
  polygon.gravity_envelope = gravity
  return polygon

def _intersect(hulls, feas_point=None):
  all_eq = np.vstack([h.equations for h in hulls])

  if feas_point is None:
    points = np.vstack([h.points[h.vertices, :] for h in hulls])
    feas_point = np.sum(points, axis=0)/points.shape[0]

  intersected = ConvexHull(HalfspaceIntersection(all_eq, feas_point).intersections)
  return intersected

class PrismIntersection():

  def __init__(self, mass, polytope, contacts, radius=1.5):
    self.mass = mass
    self.polytope = polytope
    self.contacts = contacts
    self.radius = radius

    self.initialized = False

    self.polygons = [_make_polygon(mass, s, contacts, self.radius)
                      for s in self.polytope]
    self.prisms = None
    self.figure = None
    self.axes = []
    self.threedax = None

  def initialize(self):
    for polygon in self.polygons:
      polygon.select_solver('plain')
      polygon.make_problem()
      polygon.init_algo()
      polygon.build_polys()

    self.initialized = True

  def compute(self, mode, epsilon, maxIter):
    for polygon in self.polygons:
      polygon.compute(mode, epsilon=epsilon, maxIter=maxIter, solver='plain',
          record_anim=False, plot_init=False, plot_step=False, plot_final=False)
    self.initialized = True

  def sample(self, point):
    if not self.initialized:
      self.initialize()

    assert(point.shape == (3,))

    res = []
    for s, polygon in zip(self.polytope, self.polygons):
      trans_p = point[:2] + point.item(2)*np.array([s[0,0], s[1,0]])/(s[2]+9.81)
      res.append(polygon.sample(trans_p, plot_final=False))

    return all(res)

  def polyhedron(self):
    self.prisms = []
    for s, polygon in zip(self.polytope, self.polygons):
      points = ConvexHull(polygon.polyhedron()).points
      origin = np.array([[-s[0]/(s[2]+9.81), -s[1]/(s[2]+9.81), 1]])
      top = np.hstack((points, np.zeros((points.shape[0], 1)))) + self.radius*origin
      bot = np.hstack((points, np.zeros((points.shape[0], 1)))) - self.radius*origin
      self.prisms.append(ConvexHull(np.vstack((top, bot))))

    return _intersect(self.prisms)

  def plot(self):
    self.figure = plt.figure()
    nraxes = len(self.polytope)
    nrows = 2
    ncols = nraxes // nrows + 1
    if(nraxes % nrows != 0):
      ncols += 1

    self.axes = [plt.subplot2grid((nrows, ncols), (i % nrows, i // nrows), projection="3d")
        for i in range(nraxes)]

    self.threedax = plt.subplot2grid((nrows, ncols), (0, ncols-1), rowspan=3, projection="3d")

    for polygon, ax in zip(self.polygons, self.axes):
      polygon.ax = ax
      polygon.plot_polyhedrons()

    polyhedron = self.polyhedron()

    coords = [c for c in polyhedron.points.T]
    surf = self.threedax.plot_trisurf(*coords, triangles=polyhedron.simplices, color="red", alpha=0.3, shade=True)
    surf.set_edgecolor('red')

    self.set_axes_properties()

  def set_axes_properties(self):
    for s, ax in zip(self.polytope, self.axes):
      ax.set_aspect("equal")
      ax.set_xlabel("x(m)")
      ax.set_ylabel("y(m)")
      ax.set_title("s : {}".format(s.T))

  def show(self):
    plt.show()
