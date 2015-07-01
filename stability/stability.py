import numpy as np
from cvxopt import matrix, solvers
import cdd
import shapely.geometry as geom

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

from scipy.linalg import block_diag

from collections import namedtuple

from CGAL.CGAL_Polyhedron_3 import Polyhedron_3
from CGAL.CGAL_Kernel import Point_3, Tetrahedron_3
from CGAL.CGAL_Convex_hull_3 import convex_hull_3

def cross_m(vec):
  return np.array([[0, -vec.item(2), vec.item(1)],
                   [vec.item(2), 0, -vec.item(0)],
                   [-vec.item(1), vec.item(0), 0]])

def normalize(vec):
  return vec/(np.linalg.norm(vec))

def area_convex(hrep):
  try:
    gen = np.array(cdd.Polyhedron(hrep).get_generators())
  except RuntimeError:  # Whenever numerical inconsistency is found
    return 0
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

def add_cgal(p1, p2):
  return Point_3(p1.x() + p2.x(),
                 p1.y() + p2.y(),
                 p1.z() + p2.z())

def mult_cgal(p, a):
  return Point_3(a*p.x(), a*p.y(), a*p.z())

def as_list(p):
  return [p.x(), p.y(), p.z()]

def convexify_polyhedron(hrep):
  gen = np.array(cdd.Polyhedron(hrep).get_generators())
  #If the polygon is empty or degenerate, return 0
  if gen.shape[0] < 3:
    return 0

  points = [Point_3(x, y, z) for x, y, z in gen[:, 1:]]
  poly = Polyhedron_3()
  convex_hull_3(points, poly)
  lines = [np.array([as_list(p)]) for p in poly.points()]
  return np.vstack(lines)

def tetrahedron_from_facet(facet, center):
  points = []
  he = facet.halfedge()
  points.append(he.vertex().point())
  h = he.next()
  for i in range(facet.facet_degree()-1):
    points.append(h.vertex().point())
    h = h.next()
  points.append(center)
  return Tetrahedron_3(*points)

def volume_convex(hrep):
  gen = np.array(cdd.Polyhedron(hrep).get_generators())
  #If the polygon is empty or degenerate, return 0
  if gen.shape[0] < 3:
    return 0
  #p = np.vstack([gen, gen[0, :]])
  p = gen

  points = [Point_3(x, y, z) for x, y, z in p[:, 1:]]
  poly = Polyhedron_3()
  convex_hull_3(points, poly)
  points = list(poly.points())

  center = mult_cgal(reduce(add_cgal, points),
                     1/float(len(points)))

  tetrahedrons = [tetrahedron_from_facet(f, center) for f in poly.facets()]

  return sum([abs(t.volume()) for t in tetrahedrons])

TorqueConstraint = namedtuple('TorqueConstraint',
                              ['indexes', 'point', 'limit'])

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

class SteppingException(Exception):
  def __init__(self, m):
    self.message = m

  def __str__(self):
    return self.message

#Algorithm to compute stability polygon according to
#Bretl et al. "Testing static equilibrium of legged robots"
# You need to first set some contacts and a robot mass
# Then call compute with desired precision.

class StabilityPolygon():
  def __init__(self, robotMass, gravity=-9.81):
    solvers.options['show_progress'] = False
    self.contacts = []
    self.torque_constraints = []
    self.mass = robotMass
    self.gravity = np.array([[0, 0, gravity]]).T
    self.gravity_envelope = [
        np.array([[-0.15, 0, 0]]).T,
        np.array([[0.15, 0, 0]]).T,
        np.array([[0, 0.15, 0]]).T,
        np.array([[0, -0.15, 0]]).T
                            ]
    self.radius = 2.0
    self.force_lim = 1.0
    self.proj = np.array([[1, 0, 0], [0, 1, 0]])
    self.inner = []
    self.outer = []

  def nrVars(self):
    return self.size_x() + self.size_z()

  def size_tb(self):
    if self._size_tb is None:
      return 6*len(self.torque_constraints)
    else:
      return self._size_tb

  def size_x(self):
    return 3*len(self.contacts)*len(self.gravity_envelope)

  def size_z(self):
    return 3

  def addContact(self, contact):
    self.contacts.append(contact)

  def addTorqueConstraint(self, contacts, point, limit):
    #Add a limit on torque at a point over contacts
    indexes = []
    for c in contacts:
      indexes.append(self.contacts.index(c))

    self.torque_constraints.append(TorqueConstraint(indexes, point, limit))

  def reset(self):
    self.contacts = []
    self.torque_constraints = []

  def make_problem(self):
    A_s = []
    B_s = []
    u_s = []
    L_s = []
    tb_s = []

    self._size_tb = 0

    for tc in self.torque_constraints:
      L = np.zeros((6, self.nrVars()))
      for i in tc.indexes:
        cont = self.contacts[i]
        dist = tc.point - cont.r
        L[:, 3*i:3*i+3] = np.vstack([cross_m(dist), -cross_m(dist)])
      #Filter L, tb to remove zero lines
      tb = np.vstack([tc.limit, tc.limit])
      zero_mask = np.all(L == 0, axis=1)

      L = L[~zero_mask]
      tb = tb[~zero_mask]

      L_s.append(L)
      tb_s.append(tb)

      self._size_tb += L.shape[0]

    for c in self.contacts:
      A_s.append(np.vstack([np.eye(3), cross_m(c.r)]))
      B_s.append(np.eye(3) - np.dot(c.n, c.n.T))
      u_s.append(c.mu*c.n)

    self.A1 = np.hstack(A_s)
    self.A2 = self.computeA2(self.gravity)

    self.t = self.computeT(self.gravity)

    self.B_s = B_s
    self.u_s = u_s

    self.L_s = L_s
    self.tb_s = tb_s

  def computeA2(self, gravity):
    return np.vstack([np.zeros((3, self.size_z())),
                      -cross_m(self.mass*gravity)])

  def computeT(self, gravity):
    return np.vstack([-self.mass*gravity, np.array([[0], [0], [0]])])

  def check_sizes(self):
    assert(self.A1.shape[1] + self.A2.shape[1] == self.nrVars())
    assert(self.A1.shape == (6, self.size_x()))
    assert(self.A2.shape == (6, self.size_z()))
    assert(self.t.shape == (6, 1))

  def solve(self, a):
    self.sol = self.block_socp(a, self.A1, self.A2, self.t, self.B_s, self.u_s)
    if self.sol['status'] == 'optimal':
      vec = np.array(self.sol['x'])
      self.com = vec[-self.size_z():]
      self.forces = vec[:-self.size_z()].reshape((len(self.contacts)*len(self.gravity_envelope), 3)).T
      return True
    return False

  #Compute B as diag(B_s), resulting in only one cone constraint
  def block_socp(self, a, A1, A2, t, B_s, u_s):
    dims = {
        'l': self.size_tb() + 2*self.size_x(),  # Pure inequality constraints
            # Com cone is now 3d, Size of the cones: x,y,z+1
        'q': [4]+[4]*len(self.contacts)*len(self.gravity_envelope),
        's': []  # No sd cones
            }

    size_cones = self.size_x()*4 // 3

    #Max a^T z ~ min -a^T z
    c = matrix(np.vstack([np.zeros((self.size_x(), 1)), -a]))

    A1_diag = block_diag(*([A1]*len(self.gravity_envelope)))
    A2 = np.vstack([self.computeA2(self.gravity+e)
                    for e in self.gravity_envelope])

    T = np.vstack([self.computeT(self.gravity+e)
                   for e in self.gravity_envelope])

    A = matrix(np.hstack([A1_diag, A2]))

    g_s = []
    h_s = []

    if self.L_s:
      g_s.append(np.vstack(self.L_s))
      h_s.append(np.vstack(self.tb_s))

    g_force = np.vstack([np.eye(self.size_x()), -np.eye(self.size_x())])
    g_s.append(np.hstack([g_force, np.zeros((2*self.size_x(), self.size_z()))]))

    h_s.append(self.force_lim*self.mass*9.81*np.ones((2*self.size_x(), 1)))

    #This com cone should prevent com from going to infinity : ||com|| =< max
    size_com_cone = self.size_z()+1
    com_cone = np.zeros((self.size_z()+1, self.nrVars()))
    com_cone[1, -3] = -1.0  # com_x
    com_cone[2, -2] = -1.0  # com_y
    com_cone[3, -1] = -1.0  # com_z

    g_s.append(com_cone)

    h_com_cone = np.zeros((size_com_cone, 1))
    h_com_cone[0, 0] = self.radius
    h_s.append(h_com_cone)

    #B = diag{[u_i b_i.T].T}
    blocks = [-np.vstack([u.T, B]) for u, B in zip(u_s, B_s)]*len(self.gravity_envelope)
    block = block_diag(*blocks)

    g_contacts = np.hstack([block, np.zeros((size_cones, self.size_z()))])
    g_s.append(g_contacts)
    h_cones = np.zeros((size_cones, 1))
    h_s.append(h_cones)

    for gi, hi in zip(g_s, h_s):
      print gi.shape, hi.shape

    g = np.vstack(g_s)
    h = np.vstack(h_s)

    sol = solvers.conelp(c, G=matrix(g), h=matrix(h),
                         A=A, b=matrix(T), dims=dims)
    return sol

  def socp(self, a, A1, A2, t, B_s, u_s):
    #Max a^T z ~ min -a^T z
    c = matrix(np.vstack([np.zeros((self.size_x(), 1)), -a]))

    A = matrix(np.hstack([A1, A2]))

    #Compute com friction cone
    g_com = np.zeros((self.size_z()+1, self.nrVars()))
    g_com[1, -2] = -1.0  # com_x
    g_com[2, -1] = -1.0  # com_y
    h_com = np.zeros((3, 1))
    h_com[0, 0] = 2.0
    G = [g_com]
    H = [h_com]

    #For G : compute all cones
    for i, (b, u) in enumerate(zip(B_s, u_s)):
      block = -np.vstack([u.T, b])
      g = np.hstack([np.zeros((4, 3*i)),
                     block,
                     np.zeros((4, 3*(len(self.contacts)-1-i))),
                     np.zeros((4, self.size_z()))])
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
    #Search in "random" directions
    directions = map(normalize, [np.array([[0, 1, 1]]).T,
                                 np.array([[1, 0, -1]]).T,
                                 np.array([[-1, -1, 0]]).T,
                                 np.array([[1, -1, 1]]).T])

    directions = map(normalize, [np.array([[1, 0, 0]]).T,
                                 np.array([[-1, 0, 0]]).T,
                                 np.array([[0, 1, 0]]).T,
                                 np.array([[0, -1, 0]]).T,
                                 np.array([[0, 0, 1]]).T,
                                 np.array([[0, 0, -1]]).T])

    rdirs = []

    for d in directions:
      try:
        self.step(d)
      except SteppingException as e:
        rdir = np.random.random((2, 1))
        print str(e), " Will try in a random direction {}".format(rdir.T)
        rdirs.append(rdir)

    for d in rdirs:
      self.step(d)

    assert len(self.points) >= 3, "Not enough points to form a triangle"

  def step(self, d):
    if self.solve(d):
      self.directions.append(d.T)
      self.points.append(self.com)
      self.offsets.append(d.T.dot(self.com))
    else:
      m = ["Failed to step in direction {}".format(d.T),
           "Terminated in {} state".format(self.sol['status'])]
      raise SteppingException('\n'.join(m))

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

    volumes = []

    try:
      ineq = np.array(cdd.Polyhedron(self.inner).get_inequalities())
    except RuntimeError:
      raise SteppingException('Numerical inconsistency found')

    for line in ineq:
      A_e = self.outer.copy()
      A_e.extend(cdd.Matrix(-line.reshape(1, line.size)))
      A_e.canonicalize()
      volumes.append(volume_convex(A_e))

    i, a = max(enumerate(volumes), key=lambda x: x[1])
    return -ineq[i, 1:]

  def next_edge(self, plot=False, record=False, fname_poly=False, number=0):
    d = normalize(self.find_direction().reshape((self.size_z(), 1)))
    self.step(d)

    if plot or record:
      self.plot()
      if plot:
        self.show()
      if record:
        filename = "stability_{0:04d}".format(number)
        plt.savefig(filename)
        plt.close()

    if fname_poly is not None:
      self.save_polyhedron(fname_poly+'_inner_{0:04d}'.format(number))
      self.save_outer(fname_poly+'_outer_{0:04d}'.format(number))

  def iterBound(self, nr_edges_init, error_init, prec):
    c = float(343)/float(243)
    return nr_edges_init*(np.sqrt(c*error_init/prec) - 1)

  def compute(self, epsilon=1e-4, plot_init=False,
              plot_step=False,
              record_anim=False,
              plot_final=True,
              fname_polys=None):
    self.make_problem()
    self.init_algo()
    self.build_polys()

    failure = False

    if plot_init:
      self.plot()
      self.show()

    error = volume_convex(self.outer) - volume_convex(self.inner)

    iterBound = self.iterBound(3, error, epsilon)

    print "This should take {} iterations".format(np.ceil(iterBound))

    nrSteps = 0
    while(error > epsilon):
      try:
        self.next_edge(plot_step, record_anim, fname_polys, nrSteps)
      except SteppingException as e:
        print "Failure detected... Aborting"
        print e.message
        failure = True
        break
      error = volume_convex(self.outer) - volume_convex(self.inner)
      nrSteps += 1

    print "NrIter : {} | Remainder : {}".format(nrSteps, error)

    if plot_final and not failure:
      self.plot()
      self.show()

  def polygon(self, centroid=None):
    gen = np.array(cdd.Polyhedron(self.inner).get_generators())
    p = np.vstack([gen, gen[0, :]])
    if centroid is None:
      x, y = p[:, 1], p[:, 2]
    else:
      x, y = p[:, 1] + centroid.item(0), p[:, 2] + centroid.item(1)
    return geom.Polygon(zip(x, y)).convex_hull

  def polyhedron(self):
    p = convexify_polyhedron(self.inner)
    return p

  def save_polyhedron(self, fname):
    np.savetxt(fname, self.polyhedron())

  def save_outer(self, fname):
    np.savetxt(fname, convexify_polyhedron(self.outer))

  def reset_fig(self):
    fig = plt.figure()
    self.ax = fig.add_subplot('111', aspect='equal', projection='3d')
    tup = [-1.1*self.radius, 1.1*self.radius]
    self.ax.set_xlim(tup)
    self.ax.set_ylim(tup)
    self.ax.set_zlim(tup)

    self.ax.elev = 30.
    #self.ax.dist = 20.

  def plot(self):
    self.reset_fig()
    self.plot_contacts()
    self.plot_solution()
    self.plot_polyhedrons()
    self.plot_sphere(self.radius, 'b')

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

    x, y, z = com_pos.item(0), com_pos.item(1), com_pos.item(2)
    self.ax.plot([x], [y], [z], linestyle="none",
                 marker='o', alpha=0.6,
                 markersize=10, markerfacecolor='black')

    positions = np.hstack([c.r for c in self.contacts*len(self.gravity_envelope)])

    for pos, force in zip(positions.T, forces.T):
      X, Y, Z = pos[0], pos[1], pos[2]
      U, V, W = -force[0], -force[1], -force[2]
      l = np.linalg.norm(force)/(self.mass*abs(self.gravity.item(2)))
      self.ax.quiver(X, Y, Z, U, V, W, color='r', length=l)

  def plot_direction(self, d):
    self.ax.plot([0, d.item(0)], [0, d.item(1)], marker='d')

  def plot_polyhedrons(self):
    self.plot_polyhedron(self.inner, 'red')
    self.plot_polyhedron(self.outer, 'blue')

  def plot_polyhedron(self, poly, m):
    p = convexify_polyhedron(poly)
    x, y, z = p[:, 0], p[:, 1], p[:, 2]
    self.ax.scatter(x, y, z, color=m)

  def plot_sphere(self, radius, color):
    u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
    r = radius
    x = r*np.cos(u)*np.sin(v)
    y = r*np.sin(u)*np.sin(v)
    z = r*np.cos(v)
    self.ax.plot_wireframe(x, y, z, color=color, alpha=0.2)
