import numpy as np
from scipy.linalg import block_diag
from cvxopt import matrix, solvers

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa

from constraints import TorqueConstraint, DistConstraint, ForceConstraint, CubeConstraint
from utils import cross_m
from linear_cone import rotate_axis
from printing import Verbosity, Printer

from recursive_projection import RecursiveProjectionProblem


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


class StabilityPolygon(RecursiveProjectionProblem):
  """Algorithm to compute stability polygon according to
  Bretl et al. \"Testing static equilibrium of legged robots\".
  You need to first set some contacts and a robot mass
  Then call compute with desired precision.
  In 2D, computes a static stability polygon without discretizing
  cones. In 3D, computes a robust static stability polyhedron."""

  def __init__(self, robotMass, dimension=3, gravity=-9.81,
               radius=2., force_lim=1., robust_sphere=-1, height=0.):
    """The default constructor to build a polygon/polyhedron.

    :param robotMass: Mass of the robot
    :param dimension: Number of dimensions.
    :param gravity: Intensity of gravity given along upwards z-axis.
    :param radius: Radius of the CoM limitation constraint.
    :param force_lim: Maximum force, expressed as a factor of the robot's weight.
    :param robust_sphere: Robust radius to be used with spherical criterion. Negative disables
    :param height: Height to be used when doing 2D robust
    :type dimension: 2,3
    :type gravity: double
    :type radius: double
    :type force_lim: double
    :type robust_sphere: double
    :type height: double
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
    self.printer = Printer(Verbosity.info)
    self.robust_sphere = robust_sphere

    self.nrSteps = 0

    if dimension == 3:
      shape = [
                  np.array([[-1., 0, 0]]).T,
                  np.array([[1., 0, 0]]).T,
                  np.array([[0, 1., 0]]).T,
                  np.array([[0, -1., 0]]).T
              ]
      self.gravity_envelope = [0.45*s for s in shape]
      self.proj = np.eye(3)
      self.height = 0.

    elif dimension == 2:
      self.gravity_envelope = [
          np.array([[0.0, 0.0, 0.0]]).T
                              ]
      self.proj = np.array([[1, 0, 0], [0, 1, 0]])
      self.height = height
    else:
      raise ValueError("Dimension can only be 2 or 3")

  def nrVars(self):
    return self.size_x() + self.size_z() + self.size_s()

  def size_tb(self):
    if self._size_tb is None:
      return 6*len(self.torque_constraints)
    else:
      return self._size_tb

  def size_x(self):
    return 3*len(self.contacts)*len(self.gravity_envelope)

  def size_z(self):
    return self.dimension

  def size_s(self):
    if self.robust_sphere > 0:
      return len(self.contacts)
    else:
      return 0

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

    if self.robust_sphere > 0:
      self.gravity_envelope = [
          np.array([[0.0, 0.0, 0.0]]).T
          ]

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

    if self.robust_sphere > 0:
      sigmas = np.linalg.svd(self.A1, compute_uv=False)
      self.sigma_bar = 1./min([s for s in sigmas if s > 0])
      thetas = np.linalg.svd(block_diag(*self.B_s).dot(np.linalg.pinv(self.A1)), compute_uv=False)
      self.theta_bar = max(thetas)

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

  def computeT(self, gravity, height=0):
    momentum = cross_m(self.mass*gravity).dot(np.array([[0], [0], [height]]))
    return np.vstack([-self.mass*gravity, momentum])

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
      self.com = vec[self.size_x():self.size_x()+self.size_z()]
      nrForces = len(self.contacts)*len(self.gravity_envelope)
      self.forces = vec[:self.size_x()].reshape((nrForces, 3)).T
      return self.com
    else:
      return None

  def sample(self, p, plot_final=True, plot_step=False):
    """Test if a point is stable by iteratively refining the approximations"""
    nrIter = 0
    p_t = p.reshape((p.size, 1))
    while nrIter < 50:
      if self.backend.inside(self, p):
        if plot_final:
          self.plot()
          self.ax.plot(*p_t, linestyle="none",
                        marker='o', alpha=0.6,
                        markersize=10, markerfacecolor='red')
          self.show()
        return True, nrIter
      elif self.backend.outside(self, p):
        if plot_final:
          self.plot()
          self.ax.plot(*p_t, linestyle="none",
                        marker='o', alpha=0.6,
                        markersize=10, markerfacecolor='red')
          self.show()
        return False, nrIter
      else:
        direction = self.backend.find_point_direction(self, p)
        self.step(direction.T)
        self.build_polys()
        if plot_step:
          self.plot()
          self.ax.plot(*p_t, linestyle="none",
                        marker='o', alpha=0.6,
                        markersize=10, markerfacecolor='red')
          self.show()
      nrIter += 1
    raise SteppingException("Too many steps")

  def single_test(self, p):
    """Test if a single point is robust / non-robust without refining
       the approximations"""
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
        'q': [dc.size for dc in self.dist_constraints]+[4]*len(self.contacts)*len(self.gravity_envelope)+[self.size_z()+1]*self.size_s()+[fc.size for fc in self.force_constraints],
        's': []  # No sd cones
            }
    size_cones = self.size_x()*4 // 3

    #Optimization vector [ f1 ... fn CoM s1 ... sn]

    #Max a^T z ~ min -a^T z
    #Final cost : min -a^T z + s
    c = np.vstack([np.zeros((self.size_x(), 1)), -a, np.zeros((self.size_s(), 1))])
    A1_diag = block_diag(*([A1]*len(self.gravity_envelope)))
    A2 = np.vstack([self.computeA2(self.gravity+e)
                    for e in self.gravity_envelope])

    T = np.vstack([self.computeT(self.gravity+e, height=self.height)
                   for e in self.gravity_envelope])

    A = np.hstack([A1_diag, A2, np.zeros((A2.shape[0], self.size_s()))])

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
    g_s.append(np.hstack([g_force, np.zeros((2*self.size_x(), self.size_z()+self.size_s()))]))

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

    g_contacts = np.hstack([block, np.zeros((size_cones, self.size_z()+self.size_s()))])
    h_cones = np.zeros((size_cones, 1))

    g_slacks, h_slacks = None, None
    if self.robust_sphere > 0:
      g_slacks = np.zeros((4*self.size_s(), self.nrVars()))
      h_slacks = np.zeros((4*self.size_s(), 1))
      for i, contact in enumerate(self.contacts):
        robusticity = -self.robust_sphere*self.mass*(self.theta_bar+self.sigma_bar*contact.mu)
        g_contacts[4*i, self.size_x()+self.size_z()+i] = -robusticity
        h_cones[4*i, 0] = robusticity

        #Add slack variables cones
        # || z || < s
        g_slacks[4*i, self.size_x()+self.size_z()+i] = -1
        g_slacks[4*i+1:4*i+1+self.size_z(), self.size_x():self.size_x()+self.size_z()] = -np.eye(self.size_z())

    g_s.append(g_contacts)
    h_s.append(h_cones)

    if g_slacks is not None:
      g_s.append(g_slacks)
      h_s.append(h_slacks)

    #Force Constraints
    g_s.extend(self.F_s)
    h_s.extend(self.f_s)

    #Concat everything
    g = np.vstack(g_s)
    h = np.vstack(h_s)

    sol = solvers.conelp(matrix(c), G=matrix(g), h=matrix(h),
                         A=matrix(A), b=matrix(T), dims=dims)
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

    T = np.vstack([self.computeT(self.gravity+e, height=self.height)
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

  def iterBound(self, nr_edges_init, error_init, prec):
    c = float(343)/float(243)
    return np.ceil(nr_edges_init*(np.sqrt(c*error_init/prec) - 1))

  def reset_fig(self):
    RecursiveProjectionProblem.reset_fig(self)
    if self.radius is not None:
      tup = [-1.1*self.radius, 1.1*self.radius]
      self.ax.set_xlim(tup)
      self.ax.set_ylim(tup)
      self.ax.set_zlim(tup)

  def plot(self):
    #Custom plot here because we want to plot forces too
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

  def plot_sphere(self, origin, radius, color):
    if self.size_z() == 3:
      u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
      r = radius
      x = r*np.cos(u)*np.sin(v) + origin[0]
      y = r*np.sin(u)*np.sin(v) + origin[1]
      z = r*np.cos(v) + origin[2]
      self.ax.plot_wireframe(x, y, z, color=color, alpha=0.1)
    else:
      u = np.mgrid[0:2*np.pi:40j]
      r = radius
      x = r*np.cos(u) + origin[0]
      y = r*np.sin(u) + origin[1]
      self.ax.plot(x, y, color=color, alpha=0.1)
