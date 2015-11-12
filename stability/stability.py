import numpy as np
from cvxopt import matrix, solvers
import cdd
import shapely.geometry as geom
import sys

import hashlib

from fractions import Fraction

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

from scipy.linalg import block_diag

from collections import namedtuple
from enum import Enum, unique

from CGAL.CGAL_Polyhedron_3 import Polyhedron_3
from CGAL.CGAL_Kernel import Point_3, Tetrahedron_3
from CGAL.CGAL_Convex_hull_3 import convex_hull_3

from qhull_sch import convexVolume, convexHull

from scipy.spatial import ConvexHull
from scipy.spatial.qhull import QhullError

import pyparma
from pyparma.utils import fractionize, floatize

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
  return scipy_convexify_polyhedron(hrep)

def qhull_convexify_polyhedron(hrep):
  gen = np.array(cdd.Polyhedron(hrep).get_generators())
  #If the polygon is empty or degenerate, return 0
  if gen.shape[0] < 3:
    return 0

  points = gen[:, 1:].tolist()
  pout = convexHull(points)
  return np.array(pout)

def scipy_convexify_polyhedron(hrep):
  points = np.array(cdd.Polyhedron(hrep).get_generators())[:, 1:]
  ch = ConvexHull(points)
  return ch.points

def scipy_cdd_triangulate_polyhedron(hrep):
  points = np.array(cdd.Polyhedron(hrep).get_generators())[:, 1:]
  ch = ConvexHull(points)
  return [c for c in ch.points.T], ch.simplices

def scipy_parma_triangulate_polyhedron(poly):
  points = floatize(poly.vrep()[:, 1:])
  ch = ConvexHull(points)
  return [c for c in ch.points.T], ch.simplices

def cgal_convexify_polyhedron(hrep):
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

def scipy_parma_volume_convex(poly):
  points = floatize(poly.vrep()[:, 1:])
  try:
    ch = ConvexHull(points)
  except QhullError:
    return 0
  return ch.volume

def scipy_cdd_volume_convex(hrep):
  try:
    points = np.array(cdd.Polyhedron(hrep).get_generators())[:, 1:]
  except RuntimeError:
    return 0

  if points.shape[0] < points.shape[1]+1:
    return 0

  try:
    ch = ConvexHull(points)
  except QhullError:
    return 0
  return ch.volume

def qhull_volume_convex(hrep):
  gen = np.array(cdd.Polyhedron(hrep).get_generators())
  #If the polygon is empty or degenerate, return 0
  if gen.shape[0] < 3:
    return 0
  points = gen[:, 1:].tolist()
  return convexVolume(points)

def cgal_volume_convex(hrep):
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
                              ['indexes', 'point', 'ub', 'lb'])

ForceConstraint = namedtuple('ForceConstraint',
                             ['indexes', 'limit'])

DistConstraint = namedtuple('DistConstraint',
                            ['origin', 'radius'])

@unique
class Mode(Enum):
  precision = 1
  iteration = 2
  best = 3

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
  def __init__(self, robotMass, dimension=3, gravity=-9.81):
    solvers.options['show_progress'] = False
    self.contacts = []
    self.torque_constraints = []
    self.force_constraints = []
    self.dist_constraints = []
    self.mass = robotMass
    self.gravity = np.array([[0, 0, gravity]]).T
    self.dimension = dimension

    self.volume_dic = {}
    self.vrep_dic = {}

    if dimension == 3:
      shape = [
                  np.array([[-1., 0, 0]]).T,
                  np.array([[1., 0, 0]]).T,
                  np.array([[0, 1., 0]]).T,
                  np.array([[0, -1., 0]]).T
              ]
      self.gravity_envelope = [1.45*s for s in shape]
      self.proj = np.eye(3)

    elif dimension == 2:
      self.gravity_envelope = [
          np.array([[0.0, 0.0, 0.0]]).T
                              ]
      self.proj = np.array([[1, 0, 0], [0, 1, 0]])
    else:
      raise ValueError("Dimension can only be 2 or 3")

    self.radius = 2.
    self.force_lim = 1.0
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
    return self.dimension

  def addContact(self, contact):
    self.contacts.append(contact)

  def addTorqueConstraint(self, contacts, point, ub, lb=None):
    #Add a limit on torque at a point over contacts
    indexes = []
    for c in contacts:
      indexes.append(self.contacts.index(c))

    if lb is None:
      lb = -ub

    self.torque_constraints.append(TorqueConstraint(indexes, point, ub, lb))

  def addForceConstraint(self, contacts, limit):
    """Limit the sum of forces applied on contacts"""
    indexes = []
    for c in contacts:
      indexes.append(self.contacts.index(c))

    self.force_constraints.append(ForceConstraint(indexes, limit))

  def addDistConstraint(self, origin, radius):
    """Limit the CoM to || CoM - origin || < radius"""
    self.dist_constraints.append(DistConstraint(origin, radius))

  def reset(self):
    self.contacts = []
    self.torque_constraints = []
    self.force_constraints = []
    self.dist_constraints = []
    self.volume_dic = {}
    self.vrep_dic = {}

  def make_problem(self):
    A_s = []
    B_s = []
    u_s = []
    L_s = []
    tb_s = []
    S_s = []
    r_s = []

    self._size_tb = 0

    for tc in self.torque_constraints:
      L = np.zeros((6, self.nrVars()))
      for i in tc.indexes:
        cont = self.contacts[i]
        dist = tc.point - cont.r
        off = 0
        #Add constraint on every x_i
        for j in range(len(self.gravity_envelope)):
          L[:, off+3*i:off+3*i+3] = np.vstack([cross_m(dist), -cross_m(dist)])
          off += 3*len(self.contacts)
      #Filter L, tb to remove zero lines
      tb = np.vstack([tc.ub, -tc.lb])
      zero_mask = np.all(L == 0, axis=1)

      L = L[~zero_mask]
      tb = tb[~zero_mask]

      L_s.append(L)
      tb_s.append(tb)

      self._size_tb += L.shape[0]

    for dc in self.dist_constraints:
      S = np.zeros((self.size_z()+1, self.nrVars()))
      S[1:, -self.size_z():] = -np.eye(self.size_z())
      r = np.zeros((self.size_z()+1, 1))
      r[0] = dc.radius
      r[1:] = dc.origin

      S_s.append(S)
      r_s.append(r)

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

    self.S_s = S_s
    self.r_s = r_s

  def computeA2(self, gravity):
    return np.vstack([np.zeros((3, self.size_z())),
                      -cross_m(self.mass*gravity).dot(self.proj.T)])

  def computeT(self, gravity):
    return np.vstack([-self.mass*gravity, np.array([[0], [0], [0]])])

  def check_sizes(self):
    assert(self.A1.shape[1]*len(self.gravity_envelope)
           + self.A2.shape[1]
           == self.nrVars())
    assert(self.A1.shape == (6, self.size_x() // len(self.gravity_envelope)))
    assert(self.A2.shape == (6, self.size_z()))
    assert(self.t.shape == (6, 1))

  def solve(self, a):
    self.sol = self.block_socp(a, self.A1, self.A2, self.t, self.B_s, self.u_s)
    if self.sol['status'] == 'optimal':
      vec = np.array(self.sol['x'])
      self.com = vec[-self.size_z():]
      nrForces = len(self.contacts)*len(self.gravity_envelope)
      self.forces = vec[:-self.size_z()].reshape((nrForces, 3)).T
      return True
    return False

  def sample(self, p):
    self.sol = self.check_point(p, self.A1, self.B_s, self.u_s)
    if self.sol['status'] == 'optimal':
      vec = np.array(self.sol['x'])
      self.com = p
      nrForces = len(self.contacts)*len(self.gravity_envelope)
      self.forces = vec.reshape((nrForces, 3)).T
      return True
    return False

  #Compute B as diag(B_s), resulting in only one cone constraint
  def block_socp(self, a, A1, A2, t, B_s, u_s):
    dims = {
        'l': self.size_tb() + 2*self.size_x(),  # Pure inequality constraints
            # Com cone is now 3d, Size of the cones: x,y,z+1
        'q': [self.size_z()+1]*(1+len(self.dist_constraints))+[4]*len(self.contacts)*len(self.gravity_envelope)+[4]*len(self.force_constraints),
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
    com_cone[1:, -self.size_z():] = -np.eye(self.size_z())

    g_s.append(com_cone)

    h_com_cone = np.zeros((size_com_cone, 1))
    h_com_cone[0, 0] = self.radius
    h_s.append(h_com_cone)

    #These are additional dist constraints
    g_s.append(np.vstack(self.S_s))
    h_s.append(np.vstack(self.r_s))

    #B = diag{[u_i b_i.T].T}
    blocks = [-np.vstack([u.T, B]) for u, B in zip(u_s, B_s)]*len(self.gravity_envelope)
    block = block_diag(*blocks)

    g_contacts = np.hstack([block, np.zeros((size_cones, self.size_z()))])
    g_s.append(g_contacts)
    h_cones = np.zeros((size_cones, 1))
    h_s.append(h_cones)

    #Force constraint ||\sum_i f_i || < lim
    for fc in self.force_constraints:
      force_sum = np.zeros((3, self.nrVars()))
      for index in fc.indexes:
        i = 3*index
        off = 0
        for j in range(len(self.gravity_envelope)):
          force_sum[:, off+i:off+i+3] = np.eye(3)
          off += 3*len(self.contacts)
      g_fc = np.vstack([np.zeros((1, self.nrVars())), force_sum])
      h_fc = np.zeros((4, 1))
      h_fc[0, 0] = fc.limit*self.mass*9.81
      g_s.append(g_fc)
      h_s.append(h_fc)

    g = np.vstack(g_s)
    h = np.vstack(h_s)

    sol = solvers.conelp(c, G=matrix(g), h=matrix(h),
                         A=A, b=matrix(T), dims=dims)
    return sol

  def socp(self, a, A1, A2, t, B_s, u_s):
    #Max a^T z ~ min -a^T z
    c = matrix(np.vstack([np.zeros((self.size_x(), 1)), -a]))

    A1_diag = block_diag(*([A1]*len(self.gravity_envelope)))
    A2 = np.vstack([self.computeA2(self.gravity+e)
                    for e in self.gravity_envelope])

    T = np.vstack([self.computeT(self.gravity+e)
                   for e in self.gravity_envelope])

    A = matrix(np.hstack([A1_diag, A2]))

    g_s, h_s = [], []

    if self.L_s:
      g_s.append(np.vstack(self.L_s))
      h_s.append(np.vstack(self.tb_s))

    g_force = np.vstack([np.eye(self.size_x()), -np.eye(self.size_x())])
    g_s.append(np.hstack([g_force, np.zeros((2*self.size_x(), self.size_z()))]))

    h_s.append(self.force_lim*self.mass*9.81*np.ones((2*self.size_x(), 1)))

    gl = np.vstack(g_s)
    hl = np.vstack(h_s)

    #Compute com friction cone
    g_com = np.zeros((self.size_z()+1, self.nrVars()))
    g_com[1:, -3:] = -np.eye(self.size_z())
    h_com = np.zeros((self.size_z()+1, 1))
    h_com[0, 0] = self.radius
    G = [g_com]
    H = [h_com]

    #For G : compute all cones
    for i, (b, u) in enumerate(zip(B_s, u_s)*len(self.gravity_envelope)):
      block = -np.vstack([u.T, b])
      g = np.hstack([np.zeros((4, 3*i)),
                     block,
                     np.zeros((4, 3*(len(self.contacts)*len(self.gravity_envelope)-1-i))),
                     np.zeros((4, self.size_z()))])
      G.append(g)
      H.append(np.zeros((4, 1)))

    sol = solvers.socp(c, Gl=matrix(gl), hl=matrix(hl),
                       Gq=map(matrix, G), hq=map(matrix, H),
                       A=A, b=matrix(T))
    return sol

  #Check if one point is stable or not. Based on block_socp
  #TODO: refactor away this, block_socp and regular socp
  def check_point(self, p, A1, B_s, u_s):
    dims = {
        'l': self.size_tb() + 2*self.size_x(),  # Pure inequality constraints
            # No com cone
        'q': [4]*len(self.contacts)*len(self.gravity_envelope),
        's': []  # No sd cones
            }

    size_cones = self.size_x()*4 // 3

    #Min x ~ who cares, we just want to know if there is a solution
    c = matrix(np.ones((self.size_x(), 1)))

    A1_diag = block_diag(*([A1]*len(self.gravity_envelope)))
    A2 = np.vstack([self.computeA2(self.gravity+e)
                    for e in self.gravity_envelope])

    T = np.vstack([self.computeT(self.gravity+e)
                   for e in self.gravity_envelope]) - A2.dot(p)

    A = matrix(A1_diag)

    g_s = []
    h_s = []

    if self.L_s:
      g_s.append(np.vstack(self.L_s))
      h_s.append(np.vstack(self.tb_s))

    g_force = np.vstack([np.eye(self.size_x()), -np.eye(self.size_x())])
    g_s.append(g_force)

    h_s.append(self.force_lim*self.mass*9.81*np.ones((2*self.size_x(), 1)))

    #B = diag{[u_i b_i.T].T}
    blocks = [-np.vstack([u.T, B]) for u, B in zip(u_s, B_s)]*len(self.gravity_envelope)
    block = block_diag(*blocks)

    g_s.append(block)
    h_cones = np.zeros((size_cones, 1))
    h_s.append(h_cones)

    g = np.vstack(g_s)
    h = np.vstack(h_s)

    sol = solvers.conelp(c, G=matrix(g), h=matrix(h),
                         A=A, b=matrix(T), dims=dims)
    return sol

  def init_algo(self):
    self.make_problem()
    self.directions = []
    self.points = []
    self.offsets = []
    self.inner = None
    self.outer = None
    #Search in "random" directions
    if self.dimension == 3:
      directions = map(normalize, [np.array([[1, 0, 0]]).T,
                                   np.array([[-1, 0, 0]]).T,
                                   np.array([[0, 1, 0]]).T,
                                   np.array([[0, -1, 0]]).T,
                                   np.array([[0, 0, 1]]).T,
                                   np.array([[0, 0, -1]]).T])
    elif self.dimension == 2:
      directions = map(normalize, [np.array([[1, 0]]).T,
                                   np.array([[-1, 0]]).T,
                                   np.array([[0, 1]]).T,
                                   np.array([[0, -1]]).T])

    rdirs = []

    for d in directions:
      try:
        self.step(d)
      except SteppingException as e:
        rdir = np.random.random((self.size_z(), 1))
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

    self.invalidate_vreps()

  def invalidate_vreps(self):
    offset = self.offsets[-1]
    direction = self.directions[-1]
    keys = []
    for key, vrep in self.vrep_dic.iteritems():
      valid = all(((offset+direction.dot(p[1:].T)) > 0 for p in vrep))
      if not valid:
        keys.append(key)

    #print "Invalidating {} keys out of {}".format(len(keys),
    #                                              len(self.vrep_dic.keys()))
    #Keys should always be present in both dictionaries !
    for key in keys:
      del self.volume_dic[key]
      del self.vrep_dic[key]

  def build_parma_polys(self):
    if self.outer is None:
      A = np.vstack(self.directions)
      b = np.vstack(self.offsets)
      self.outer = pyparma.Polyhedron(hrep=fractionize(np.hstack([b, -A])))
    else:
      self.outer.add_ineq(fractionize(np.hstack((self.offsets[-1],
                                                 -self.directions[-1]))))
    if self.inner is None:
      A = np.vstack([p.T for p in self.points])
      b = np.ones((len(self.points), 1))
      self.inner = pyparma.Polyhedron(vrep=np.hstack([b, fractionize(A)]))
    else:
      self.inner.add_generator(np.hstack(([[1]],
                                          fractionize(self.points[-1].T))))

  def build_cdd_polys(self):
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

  def find_parma_direction(self, plot=False):
    self.build_parma_polys()
    volumes = []

    ineq = self.inner.hrep()

    for line in ineq:
      key = hashlib.sha1(line).hexdigest()
      if key in self.volume_dic:
        volumes.append(self.volume_dic[key])
      else:
        A_e = self.outer.copy()
        A_e.add_ineq(-line)
        vol = self.volume_convex(A_e)
        self.volume_dic[key] = vol
        volumes.append(vol)
        self.vrep_dic[key] = A_e.vrep()

      if plot:
        self.reset_fig()
        self.plot_polyhedrons()
        self.plot_polyhedron(A_e, 'm', 0.5)
        self.show()

    i, a = max(enumerate(volumes), key=lambda x: x[1])
    return floatize(-ineq[i, 1:])

  def find_cdd_direction(self, plot=False):
    self.build_cdd_polys()

    volumes = []

    try:
      ineq = np.array(cdd.Polyhedron(self.inner).get_inequalities())
    except RuntimeError:
      raise SteppingException('Numerical inconsistency found')

    for line in ineq:
      A_e = self.outer.copy()
      A_e.extend(cdd.Matrix(-line.reshape(1, line.size)))
      A_e.canonicalize()

      if plot:
        self.reset_fig()
        self.plot_polyhedrons()
        self.plot_polyhedron(A_e, 'm', 0.5)
        self.show()

      volumes.append(self.volume_convex(A_e))

    i, a = max(enumerate(volumes), key=lambda x: x[1])
    return -ineq[i, 1:]

  def next_edge(self, plot=False, plot_direction=False,
                record=False, fname_poly=False, number=0):
    d = normalize(self.find_direction(plot_direction).reshape((self.size_z(), 1)))
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

  def select_solver(self, solver):
    if solver == 'cdd':
      self.volume_convex = scipy_cdd_volume_convex
      self.build_polys = self.build_cdd_polys
      self.find_direction = self.find_cdd_direction
      self.triangulate_polyhedron = scipy_cdd_triangulate_polyhedron
    elif solver == 'parma':
      self.volume_convex = scipy_parma_volume_convex
      self.build_polys = self.build_parma_polys
      self.find_direction = self.find_parma_direction
      self.triangulate_polyhedron = scipy_parma_triangulate_polyhedron
    else:
      raise ValueError("Only 'cdd' or 'parma' solvers are available")

  def compute(self, mode, maxIter=100, epsilon=1e-4,
              solver='cdd',
              plot_error=False,
              plot_init=False,
              plot_step=False,
              plot_direction=False,
              record_anim=False,
              plot_final=True,
              fname_polys=None):

    self.select_solver(solver)
    self.make_problem()
    self.init_algo()
    self.build_polys()
    failure = False

    if plot_init:
      self.plot()
      self.show()

    error = self.volume_convex(self.outer) - self.volume_convex(self.inner)
    nrSteps = 0

    if plot_error:
      self.fig_error = plt.figure()
      self.ax_error = self.fig_error.add_subplot(111)
      self.line_error, = self.ax_error.plot([nrSteps], [error], 'r-')

    if(mode is Mode.precision):
      iterBound = self.iterBound(len(self.points), error, epsilon)
      print "Reaching {} should take {} iterations".format(epsilon,
                                                           np.ceil(iterBound))
      stop_condition = lambda: error > epsilon
    elif(mode is Mode.iteration):
      print "This will take {} iterations".format(maxIter)
      stop_condition = lambda: nrSteps < maxIter
    elif(mode is Mode.best):
      print "Will try to reach {} under {} iterations".format(epsilon, maxIter)
      stop_condition = lambda: nrSteps < maxIter and error > epsilon
    else:
      raise ValueError("Unknown mode, please use a value supplied in enum")

    while(stop_condition()):
      try:
        self.next_edge(plot_step, plot_direction, record_anim,
                       fname_polys, nrSteps)
      except SteppingException as e:
        print "Failure detected... Aborting"
        print e.message
        failure = True
        break
      error = self.volume_convex(self.outer) - self.volume_convex(self.inner)
      print error
      sys.stdout.flush()

      nrSteps += 1
      if plot_error:
        self.line_error.set_xdata(np.append(self.line_error.get_xdata(),
                                            nrSteps))
        self.line_error.set_ydata(np.append(self.line_error.get_ydata(),
                                            error))
        self.ax_error.relim()
        self.ax_error.autoscale_view()
        self.fig_error.canvas.draw()
        plt.pause(0.01)

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
    self.plot_sphere(np.zeros((3, 1)), self.radius, 'b')

    for dc in self.dist_constraints:
      self.plot_sphere(dc.origin, dc.radius, 'b')

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

    self.ax.plot(*com_pos, linestyle="none",
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
    self.plot_polyhedron(self.inner, 'red', 0.5)
    self.plot_polyhedron(self.outer, 'blue', 0.1)

  def plot_polyhedron(self, poly, c, a):
    coords, tri = self.triangulate_polyhedron(poly)
    if len(coords) == 3:
      self.ax.plot_trisurf(*coords, triangles=tri, color=c, alpha=a)
    else:
      self.ax.plot(*coords, linestyle='+', color=c, alpha=1.0)
      self.ax.triplot(*coords, triangles=tri, color=c, alpha=a)

  def plot_sphere(self, origin, radius, color):
    if self.size_z() == 3:
      u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
      r = radius
      x = r*np.cos(u)*np.sin(v) + origin[0]
      y = r*np.sin(u)*np.sin(v) + origin[1]
      z = r*np.cos(v) + origin[2]
      self.ax.plot_wireframe(x, y, z, color=color, alpha=0.2)
    else:
      u = np.mgrid[0:2*np.pi:40j]
      r = radius
      x = r*np.cos(u) + origin[0]
      y = r*np.sin(u) + origin[1]
      self.ax.plot(x, y, color=color, alpha=0.2)
