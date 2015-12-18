import numpy as np
from cvxopt import matrix, solvers
import cdd
import shapely.geometry as geom
import sys

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

from scipy.linalg import block_diag

from enum import Enum, unique

from constraints import TorqueConstraint, DistConstraint, ForceConstraint, CubeConstraint

from backends import CDDBackend, ParmaBackend, PlainBackend

from functools import partial

from utils import cross_m, normalize

from geomengines import convexify_polyhedron

from linear_cone import rotate_axis

@unique
class Mode(Enum):

  """Enum that shows the three modes of functionment."""

  precision = 1
  iteration = 2
  best = 3

class Contact():

  """Class representing a contact as a single point"""

  def __init__(self, mu, pos, normal):
    """Default constructor.

    :param r: position in the world frame
    :param n: associated normal
    :param mu: friction coefficient
    :type n: np.array(3,1)
    :type r: np.array(3,1)
    :type mu: double"""

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

class StabilityPolygon():
  """Algorithm to compute stability polygon according to
  Bretl et al. \"Testing static equilibrium of legged robots\".
  You need to first set some contacts and a robot mass
  Then call compute with desired precision.
  In 2D, computes a static stability polygon without discretizing
  cones. In 3D, computes a robust static stability polyhedron."""

  def __init__(self, robotMass, dimension=3, gravity=-9.81,
               radius=2., force_lim=1.):
    """The default constructor to build a polygon/polyhedron.

    :param robotMass: Mass of the robot
    :param dimension: Number of dimensions.
    :param gravity: Intensity of gravity given along upwards z-axis.
    :param radius: Radius of the CoM limitation constraint.
    :param force_lim: Maximum force, expressed as a factor of the robot's weight.
    :type dimension: 2,3
    :type gravity: double
    """
    solvers.options['show_progress'] = False
    self.contacts = []

    self.clearConstraints()
    self.clearAlgo()

    self.mass = robotMass
    self.gravity = np.array([[0, 0, gravity]]).T
    self.dimension = dimension
    self.radius = radius
    self.force_lim = force_lim

    if dimension == 3:
      shape = [
                  np.array([[-1., 0, 0]]).T,
                  np.array([[1., 0, 0]]).T,
                  np.array([[0, 1., 0]]).T,
                  np.array([[0, -1., 0]]).T
              ]
      self.gravity_envelope = [0.45*s for s in shape]
      self.proj = np.eye(3)

    elif dimension == 2:
      self.gravity_envelope = [
          np.array([[0.0, 0.0, 0.0]]).T
                              ]
      self.proj = np.array([[1, 0, 0], [0, 1, 0]])
    else:
      raise ValueError("Dimension can only be 2 or 3")

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
    """Add a limit on torque at a point over a set of contacts"""
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

  def addCubeConstraint(self, origin, length):
    """Limit the CoM to | CoM - origin | < length"""
    self.cube_constraints.append(CubeConstraint(origin, length))

  def reset(self):
    """Remove all contacts, constraints and resets inner state"""
    self.contacts = []
    self.clearAlgo()
    self.clearConstraints()

  def clearConstraints(self):
    """Remove all constraints"""
    self.torque_constraints = []
    self.force_constraints = []
    self.dist_constraints = []
    self.cube_constraints = []

  def clearAlgo(self):
    """Resets internal state"""
    self.volume_dic = {}
    self.vrep_dic = {}

    self.directions = []
    self.points = []
    self.offsets = []

    self.inner = None
    self.outer = None

  def make_problem(self):
    """Compute all matrices necessary to solving the problems. Only needs to be
    called once, because the problem shape never changes. This adds global dist
    constraint that should prevent CoM from going to infinity : ||com|| =< max.
    However, make sure you remove it between calls to compute or to set it to
    None when creating the polygon."""

    A_s = []
    B_s = []
    C_s = []
    u_s = []
    L_s = []
    tb_s = []
    S_s = []
    r_s = []
    F_s = []
    f_s = []
    d_s = []

    if self.radius is not None:
      self.addDistConstraint(np.zeros((self.size_z(), 1)), self.radius)

    self._size_tb = 0

    if self.torque_constraints:
      for tc in self.torque_constraints:
        tc.compute(self)
      L_s, tb_s = zip(*[tc.matrices() for tc in self.torque_constraints])
      self._size_tb = sum(L.shape[0] for L in L_s)
    else:
      L_s, tb_s = [], []

    if self.dist_constraints:
      for dc in self.dist_constraints:
        dc.compute(self)
      S_s, r_s = zip(*[dc.matrices() for dc in self.dist_constraints])
    else:
      S_s, r_s = [], []

    if self.force_constraints:
      for fc in self.force_constraints:
        fc.compute(self)
      F_s, f_s = zip(*[fc.matrices() for fc in self.force_constraints])
    else:
      F_s, f_s = [], []

    if self.cube_constraints:
      for cc in self.cube_constraints:
        cc.compute(self)
      C_s, d_s = zip(*[fc.matrices() for fc in self.cube_constraints])
    else:
      C_s, d_s = [], []

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

    self.F_s = F_s
    self.f_s = f_s

    self.C_s = C_s
    self.d_s = d_s

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
        'l': self.size_tb() + 2*self.size_x() + sum([cc.size for cc in self.cube_constraints]),  # Pure inequality constraints
            # Com cone is now 3d, Size of the cones: x,y,z+1
        'q': [dc.size for dc in self.dist_constraints]+[4]*len(self.contacts)*len(self.gravity_envelope)+[fc.size for fc in self.force_constraints],
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

    #Torque constraint
    if self.L_s:
      g_s.extend(self.L_s)
      h_s.extend(self.tb_s)

    #Linear force constraint
    #       <x> <z>                      <1>
    # <x> [[ I  0] <= <x> lim*mg  2*<x> [ 1]
    # <x> [[ -I 0]
    g_force = np.vstack([np.eye(self.size_x()), -np.eye(self.size_x())])
    g_s.append(np.hstack([g_force, np.zeros((2*self.size_x(), self.size_z()))]))

    h_s.append(self.force_lim*self.mass*9.81*np.ones((2*self.size_x(), 1)))

    #Cube constraints
    if self.C_s:
      g_s.append(np.vstack(self.C_s))
      h_s.append(np.vstack(self.d_s))

    #These are additional dist constraints
    g_s.extend(self.S_s)
    h_s.extend(self.r_s)

    #B = diag{[u_i b_i.T].T}
    blocks = [-np.vstack([u.T, B]) for u, B in zip(u_s, B_s)]*len(self.gravity_envelope)
    block = block_diag(*blocks)

    g_contacts = np.hstack([block, np.zeros((size_cones, self.size_z()))])
    g_s.append(g_contacts)
    h_cones = np.zeros((size_cones, 1))
    h_s.append(h_cones)

    #Force Constraints
    g_s.extend(self.F_s)
    h_s.extend(self.f_s)

    #Concat everything
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
    self.clearAlgo()
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
      self.backend = CDDBackend('scipy')
    elif solver == 'parma':
      self.backend = ParmaBackend('scipy')
    elif solver == 'plain':
      if self.dimension != 2:
        raise ValueError("The plain backend can only be used in 2D")
      self.backend = PlainBackend('scipy')
    else:
      raise ValueError("Only 'cdd' or 'parma' solvers are available")
    self.volume_convex = self.backend.volume_convex
    self.build_polys = partial(self.backend.build_polys, self)
    self.find_direction = partial(self.backend.find_direction, self)
    self.triangulate_polyhedron = self.backend.scipy_triangulate_polyhedron

  def compute(self, mode, maxIter=100, epsilon=1e-4,
              solver='cdd',
              plot_error=False,
              plot_init=False,
              plot_step=False,
              plot_direction=False,
              record_anim=False,
              plot_final=True,
              fname_polys=None):
    """Compute the polygon/polyhedron at a given precision or for a number of
    iterations or at best.

    :param mode: Stopping criterion.
    :param maxIter: Maximum number of iterations.
    :param epsilon: Precision target.
    :param solver: Backend that will be used.
    :param plot_error: Make a running plot of the error during computation.
    :param plot_init: Plot the initial state of the algorithm.
    :param plot_step: Plot the state of the algorithm at each iteration.
    :param plot_direction: Plot the direction found for the next iteration.
    :param record_anim: Record all steps as images in files.
    :param plot_final: Plot the final polygon/polyhedron.
    :param fname_polys: Record successive iterations as text files.
    :type mode: stability.Mode
    :type maxIter: int
    :type precision: double
    :type solver: stability.backends
    """

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
    """Return the inner approximation as a shapely polygon. Only valid in 2D"""
    assert(self.dimension == 2)
    if isinstance(self.backend, CDDBackend):
      gen = np.array(cdd.Polyhedron(self.inner).get_generators())[:, 1:]
    elif isinstance(self.backend, PlainBackend):
      gen = self.inner.vertices
    elif isinstance(self.backend, ParmaBackend):
      gen = self.inner.vrep()[:, 1:]
    else:
      raise NotImplemented("No polygon method defined for this backend")

    p = np.vstack([gen, gen[0, :]])
    if centroid is None:
      x, y = p[:, 0], p[:, 1]
    else:
      x, y = p[:, 0] + centroid.item(0), p[:, 1] + centroid.item(1)
    return geom.Polygon(zip(x, y)).convex_hull

  def polyhedron(self):
    """Return the inner polyhedron as a set of points"""
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

    for dc in self.dist_constraints:
      self.plot_sphere(dc.origin, dc.radius, 'b')

  def show(self):
    plt.show()

  def plot_contacts(self):
    num_points = 50
    radius = 0.1
    angles = np.linspace(0, 2*np.pi, num_points)
    points = np.vstack((radius*np.cos(angles), radius*np.sin(angles), np.zeros((num_points,)))).T

    for c in self.contacts:
      rot = rotate_axis(np.array([[0., 0., 1.]]).T, c.n)
      world_points = radius*c.n+c.r+rot.dot(c.mu*points.T)
      x, y, z = zip(*world_points.T)
      self.ax.plot(x, y, z, 'black')
      self.ax.plot(c.r[0], c.r[1], c.r[2], marker='o', markersize=6, color='black')

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
      self.ax.quiver(X, Y, Z, U, V, W, color='r', length=l, arrow_length_ratio=0.)

  def plot_direction(self, d):
    self.ax.plot([0, d.item(0)], [0, d.item(1)], marker='d')

  def plot_polyhedrons(self):
    self.plot_polyhedron(self.inner, 'red', 0.1)
    self.plot_polyhedron(self.outer, 'blue', 0.1)

  def plot_polyhedron(self, poly, c, a):
    coords, tri = self.triangulate_polyhedron(poly)
    if len(coords) == 3:
      surf = self.ax.plot_trisurf(*coords, triangles=tri, color=c, alpha=a, shade=True)
      surf.set_edgecolor('red')
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
