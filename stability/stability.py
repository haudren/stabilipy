import numpy as np
from cvxopt import matrix, solvers
import cdd
import shapely.geometry as geom

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def cross_m(vec):
  return np.array([[-0, -vec.item(2), vec.item(1)],
                   [vec.item(2), 0, -vec.item(0)],
                   [-vec.item(1), vec.item(0), 0]])

def normalize(vec):
  return vec/(np.linalg.norm(vec))

def segments(p):
      return zip(p, p[1:] + [p[0]])

def area(hrep):
  points = [(p.item(0), p.item(1)) for p in hrep.generators]
  return 0.5 * abs(sum(x0*y1 - x1*y0
                       for ((x0, y0), (x1, y1)) in segments(points)))

def area_convex(hrep):
  gen = np.array(cdd.Polyhedron(hrep).get_generators())
  #If the polygon is empty or degenerate, return 0
  if gen.shape[0] < 3:
    return 0
  p = np.vstack([gen, gen[0, :]])
  x, y = p[:, 1], p[:, 2]
  poly = geom.Polygon(zip(x, y))
  return poly.convex_hull.area

def convex_hull(x_p, y_p):
  poly = geom.Polygon(zip(x_p, y_p))
  x_c, y_c = poly.convex_hull.exterior.xy
  return x_c, y_c

def convexify(points):
  A = np.hstack(points)
  x_p, y_p = A[0, :], A[1, :]
  x_c, y_c = convex_hull(x_p, y_p)
  p_c = [np.array([[x], [y]]) for x, y in zip(x_c, y_c)]
  return p_c

def convexify_h(hrep):
  A = hrep.generators
  x_p, y_p = A[:, 0], A[:, 1]
  x_c, y_c = convex_hull(x_p, y_p)
  return Vrep(np.vstack([x_c, y_c]).T)

#A contact should contain :
# - r : position world
# - n : normal to teh surface
# - mu : friction coefficient

class Contact():
  def __init__(self, mu, pos, normal):
    self.n = normal
    self.r = pos
    self.mu = mu

  def check_size(self):
    assert self.n.shape == (3, 1),\
        "Normal is {} shape instead of (3,1))".format(self.n.shape)
    assert self.r.shape == (3, 1),\
        "Position is {} shape instead of (3,1))".format(self.r.shape)

  def __str__(self):
    lines = ["Contact at :", str(self.r.T), "With normal :", str(self.n.T)]
    return "\n".join(lines)

#Algorithm to compute stability polygon according to
#Bretl et al. "Testing static equilibrium of legged robots"
# You need to first set some contacts and a robot mass
# Then call compute with desired precision.

class StabilityPolygon():
  def __init__(self, robotMass, gravity=-9.81):
    solvers.options['show_progress'] = False
    self.contacts = []
    self.mass = robotMass
    self.gravity = np.array([[0, 0, gravity]]).T
    self.proj = np.array([[1, 0, 0], [0, 1, 0]])
    self.inner = []
    self.outer = []

  def nrVars(self):
    return self.size_x() + self.size_z()

  def size_x(self):
    return 3*len(self.contacts)

  def size_z(self):
    return 2

  def addContact(self, contact):
    self.contacts.append(contact)

  def reset(self):
    self.contacts = []

  def make_problem(self):
    A_s = []
    B_s = []
    u_s = []

    for c in self.contacts:
      A_s.append(np.vstack([np.eye(3), cross_m(c.r)]))
      B_s.append(np.eye(3) - np.dot(c.n, c.n.T))
      u_s.append(c.mu*c.n)

    self.A1 = np.hstack(A_s)
    self.A2 = np.vstack([np.zeros((3, 2)),
                         -cross_m(self.mass*self.gravity).dot(self.proj.T)])
    self.t = np.vstack([-self.mass*self.gravity, np.array([[0], [0], [0]])])

    self.B_s = B_s
    self.u_s = u_s

  def check_sizes(self):
    assert(self.A1.shape[1] + self.A2.shape[1] == self.nrVars())

  def solve(self, a):
    self.sol = self.socp(a, self.A1, self.A2, self.t, self.B_s, self.u_s)
    if self.sol['status'] == 'optimal':
      vec = np.array(self.sol['x'])
      self.com = vec[-2:]
      self.forces = vec[:-2].reshape((3, len(self.contacts))).T
      return True
    return False

  def socp(self, a, A1, A2, t, B_s, u_s):
    #Max a^T z ~ min -a^T z
    c = matrix(np.vstack([np.zeros((self.size_x(), 1)), -a]))

    A = matrix(np.hstack([A1, A2]))

    #For G : compute all cones
    G = []
    H = []
    for i, (b, u) in enumerate(zip(B_s, u_s)):
      block = np.hstack([-np.vstack([u.T, b]),
                        np.zeros((4, 2))])
      g = np.hstack([np.zeros((4, 3*i)),
                     block,
                     np.zeros((4, 3*(len(B_s)-1-i)))])
      G.append(g)
      H.append(np.zeros((4, 1)))

    sol = solvers.socp(c, Gq=map(matrix, G), hq=map(matrix, H),
                       A=A, b=matrix(t))
    return sol

  def init_algo(self):
    self.make_problem()
    self.directions = []
    self.points = []
    self.offsets = []
    self.inner = None
    self.outer = None
    #Search in three random directions
    directions = map(normalize, [np.array([[0, 1]]).T,
                                 np.array([[1, 0]]).T,
                                 np.array([[-1, -1]]).T])
    for d in directions:
      self.step(d)

  def step(self, d):
    if self.solve(d):
      self.directions.append(d.T)
      self.points.append(self.com)
      self.offsets.append(d.T.dot(self.com))
    else:
      raise Exception("Failed to step in direction {}".format(d))

  def build_polys(self):
    if self.outer is None:
      A = np.vstack(self.directions)
      b = np.vstack(self.offsets)
      self.outer = cdd.Matrix(np.hstack([b, -A]))
      self.outer.rep_type = cdd.RepType.INEQUALITY
    else:
      self.outer.extend(np.hstack((self.offsets[-1], -self.directions[-1])))
      self.outer.canonicalize()

    if self.inner is None:
      A = np.vstack([p.T for p in self.points])
      b = np.ones((len(self.points), 1))
      self.inner = cdd.Matrix(np.hstack((b, A)))
      self.inner.rep_type = cdd.RepType.GENERATOR
    else:
      self.inner.extend(np.hstack(([[1]], self.points[-1].T)))
      self.inner.canonicalize()

  def find_direction(self):
    self.build_polys()
    out_area = area_convex(self.outer)

    areas = []
    ineq = np.array(cdd.Polyhedron(self.inner).get_inequalities())

    for line in ineq:
      A_e = self.outer.copy()
      A_e.extend(cdd.Matrix(line.reshape(1, line.size)))
      A_e.canonicalize()
      areas.append(out_area - area_convex(A_e))

    i, a = max(enumerate(areas), key=lambda x: x[1])
    return -ineq[i, 1:]

  def next_edge(self, plot=False, record=False, number=0):
    d = normalize(self.find_direction().reshape((2, 1)))
    self.step(d)

    if plot or record:
      self.plot()
      if plot:
        self.show()
      if record:
        filename = "stability_{0:04d}".format(number)
        plt.savefig(filename)
        plt.close()

  def iterBound(self, nr_edges_init, error_init, prec):
    c = float(343)/float(243)
    return nr_edges_init*(np.sqrt(c*error_init/prec) - 1)

  def compute(self, epsilon=1e-4, plot_init=False,
              plot_step=False,
              record_anim=False,
              plot_final=True):
    self.make_problem()
    self.init_algo()
    self.build_polys()

    if plot_init:
      self.plot()
      self.show()

    error = area_convex(self.outer) - area_convex(self.inner)

    iterBound = self.iterBound(3, error, epsilon)

    print "This should take {} iterations".format(np.ceil(iterBound))

    nrSteps = 0
    while(error > epsilon):
      print error
      self.next_edge(plot_step, record_anim, nrSteps)
      error = area_convex(self.outer) - area_convex(self.inner)
      nrSteps += 1

    print "NrIter : {} | Remainder : {}".format(nrSteps, error)

    if plot_final:
      self.plot()
      self.show()

  def polygon(self, centroid=None):
    gen = np.array(cdd.Polyhedron(self.inner).get_generators())
    p = np.vstack([gen, gen[0, :]])
    if centroid is None:
      x, y = p[:, 1], p[:, 2]
    else:
      x, y = p[:, 1] + centroid.item(0), p[:, 2] + centroid.item(1)
    return geom.Polygon(zip(x, y))

  def plot(self):
    fig = plt.figure()
    self.ax = fig.add_subplot('111', aspect='equal', projection='3d')
    self.ax.set_xlim([-0.5, 1])
    self.ax.set_ylim([-0.5, 1])
    self.ax.set_zlim([-0.5, 1])

    self.plot_contacts()
    self.plot_solution()
    self.plot_polygons()

  def show(self):
    plt.show()

  def plot_contacts(self):
    X, Y, Z, U, V, W = [], [], [], [], [], []
    positions = np.hstack([c.r for c in self.contacts])
    normals = np.hstack([c.n for c in self.contacts])
    X, Y, Z = positions[0, :], positions[1, :], positions[2, :]
    U, V, W = normals[0, :], normals[1, :], normals[2, :]
    self.ax.quiver(X, Y, Z, U, V, W, color='black', linestyles='dashed')

  def plot_solution(self):
    com_pos = self.com
    forces = self.forces

    x, y = com_pos.item(0), com_pos.item(1)
    self.ax.plot([x], [y], linestyle="none",
                 marker='o', alpha=0.6,
                 markersize=10, markerfacecolor='black')

    positions = np.hstack([c.r for c in self.contacts])

    for pos, force in zip(positions.T, forces.T):
      X, Y, Z = pos[0], pos[1], pos[2]
      U, V, W = -force[0], -force[1], -force[2]
      l = np.linalg.norm(force)/(self.mass*abs(self.gravity.item(2)))
      self.ax.quiver(X, Y, Z, U, V, W, color='r', length=l)

  def plot_direction(self, d):
    self.ax.plot([0, d.item(0)], [0, d.item(1)], marker='d')

  def plot_polygons(self):
    self.plot_polygon(self.inner, 'x')
    self.plot_polygon(self.outer, 'o')

  def plot_polygon(self, poly, m):
    gen = np.array(cdd.Polyhedron(poly).get_generators())
    p = np.vstack([gen, gen[0, :]])
    x, y = p[:, 1], p[:, 2]
    poly = geom.Polygon(zip(x, y))
    x, y = poly.convex_hull.exterior.xy
    self.ax.plot(x, y, marker=m)
