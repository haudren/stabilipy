#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2015-2017 CNRS-AIST JRL

# This file is part of StabiliPy.

# StabiliPy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# StabiliPy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with StabiliPy.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import absolute_import
from builtins import map
from builtins import str
from builtins import object
from .backends import CDDBackend, PlainBackend, QhullBackend, ParmaBackend
from functools import partial
from enum import Enum, unique

from .printing import Verbosity, Printer
import numpy as np
import sys
import matplotlib.pyplot as plt

from .utils import normalize


@unique
class Mode(Enum):

  """Enum that shows the three modes of functionment."""

  precision = 1
  iteration = 2
  best = 3


class SteppingException(Exception):
  def __init__(self, m):
    self.message = m

  def __str__(self):
    return self.message


class RecursiveProjectionProblem(object):

  """Base class that encapsulates the recursive projection algorithm.
     To use it, you need to specify your problem by implementing the
     `solve` method. Then, this class will actually perform the projection."""

  def __init__(self, dimension, verbosity=Verbosity.info):
    """Construct a projection problem.

       :param dimension: Dimension of the space on which you project.
       :param verbosity: Verbosity of the output. Default to `info`.
       :type dimension: int
       :type verbosity: Verbosity"""

    self.dimension = dimension
    self.printer = Printer(verbosity)

  def solve(self, d):
    """This method should return an extremal point in the given direction d,
       or None in case of error. You must reimplement this function to
       perform a computation or use one of the pre-implemented instances

       :param d: Search direction
       :type d: np.array((dim, 1))
       """

    raise NotImplementedError("You must re-implement this function")

  def step(self, d):
    extremal = self.solve(d)
    if extremal is not None:
      self.directions.append(d.T)
      self.points.append(extremal)
      self.offsets.append(d.T.dot(extremal))
    else:
      m = ["Failed to step in direction {}".format(d.T),
           "Terminated in {} state".format(self.sol['status'])]
      raise SteppingException('\n'.join(m))

    if self.dimension >= 3:
      self.invalidate_vreps()

  def select_solver(self, solver):
    if solver == 'cdd':
      self.backend = CDDBackend('scipy')
    elif solver == 'parma':
      self.backend = ParmaBackend('scipy')
    elif solver == 'plain':
      if self.dimension != 2:
        raise ValueError("The plain backend can only be used in 2D")
      self.backend = PlainBackend('scipy')
    elif solver == 'qhull':
      self.backend = QhullBackend('scipy')
    else:
      raise ValueError("Only 'cdd' or 'parma' solvers are available")
    self.volume_convex = self.backend.volume_convex
    self.build_polys = partial(self.backend.build_polys, self)
    self.find_direction = partial(self.backend.find_direction, self)
    self.triangulate_polyhedron = self.backend.scipy_triangulate_polyhedron
    self.convexify_polyhedron = self.backend.scipy_convexify_polyhedron

    if self.dimension == 3:
      self.invalidate_vreps = partial(self.backend.invalidate_vreps, self)

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
    :type mode: stabilipy.Mode
    :type maxIter: int
    :type precision: double
    :type solver: stabilipy.backends
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
    self.nrSteps = nrSteps

    if plot_error:
      self.fig_error = plt.figure()
      self.ax_error = self.fig_error.add_subplot(111)
      self.line_error, = self.ax_error.plot([nrSteps], [error], 'r-')

    if(mode is Mode.precision):
      iterBound = self.iterBound(len(self.points), error, epsilon)
      self.printer("Reaching {} should take {} iterations".format(epsilon, iterBound),
                   Verbosity.info)
      stop_condition = lambda: error > epsilon
    elif(mode is Mode.iteration):
      self.printer("This will take {} iterations".format(maxIter), Verbosity.info)
      stop_condition = lambda: nrSteps < maxIter
    elif(mode is Mode.best):
      self.printer("Will try to reach {} under {} iterations".format(epsilon, maxIter),
                   Verbosity.info)
      stop_condition = lambda: nrSteps < maxIter and error > epsilon
    else:
      raise ValueError("Unknown mode, please use a value supplied in enum")

    while(stop_condition()):
      try:
        self.next_edge(plot_step, plot_direction, record_anim,
                       fname_polys, nrSteps)
      except SteppingException as e:
        self.printer("Failure detected... Aborting", Verbosity.error)
        self.printer(e.message, Verbosity.error)
        failure = True
        break
      error = self.volume_convex(self.outer) - self.volume_convex(self.inner)
      self.printer("{} - {} = {} ".format(self.volume_convex(self.outer),
                                          self.volume_convex(self.inner),
                                          error),
                   Verbosity.debug)
      self.printer(error, Verbosity.info)
      sys.stdout.flush()

      nrSteps += 1
      self.nrSteps = nrSteps
      if plot_error:
        self.line_error.set_xdata(np.append(self.line_error.get_xdata(),
                                            nrSteps))
        self.line_error.set_ydata(np.append(self.line_error.get_ydata(),
                                            error))
        self.ax_error.relim()
        self.ax_error.autoscale_view()
        self.fig_error.canvas.draw()
        plt.pause(0.01)

    self.printer("NrIter : {} | Remainder : {}".format(nrSteps, error), Verbosity.info)

    if plot_final and not failure:
      self.plot()
      self.show()

  def make_problem(self):
    """This method is called upon launching the computation. Use it to
       build the complete problem from user-defined quantities. Does nothing
       by default."""
    return

  def init_algo(self):
    self.clearAlgo()

    #Search in "random" directions
    if self.dimension == 3:
      directions = list(map(normalize, [np.array([[1, 0, 0]]).T,
                                   np.array([[-1, 0, 0]]).T,
                                   np.array([[0, 1, 0]]).T,
                                   np.array([[0, -1, 0]]).T,
                                   np.array([[0, 0, 1]]).T,
                                   np.array([[0, 0, -1]]).T]))
    elif self.dimension == 2:
      directions = list(map(normalize, [np.array([[1, 0]]).T,
                                   np.array([[-1, 0]]).T,
                                   np.array([[0, 1]]).T,
                                   np.array([[0, -1]]).T]))

    rdirs = []

    for d in directions:
      try:
        self.step(d)
      except SteppingException as e:
        rdir = np.random.random((self.size_z(), 1))
        self.printer(str(e)+" Will try in a random direction {}".format(rdir.T),
                     Verbosity.error)
        rdirs.append(rdir)

    for d in rdirs:
      self.step(d)

    assert len(self.points) >= self.dimension+1, "Not enough points to form a simplex"

  def clearAlgo(self):
    """Resets internal state"""
    self.volume_dic = {}
    self.vrep_dic = {}
    self.hrep_dic = {}
    self.hull_dic = {}

    self.directions = []
    self.points = []
    self.offsets = []

    self.inner = None
    self.outer = None

  def next_edge(self, plot=False, plot_direction=False,
                record=False, fname_poly=None, number=0):
    d = normalize(self.find_direction(plot_direction).reshape((self.dimension, 1)))
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

  def plot(self):
    """Plot the current solution and polyhedrons"""
    if self.dimension <= 3:
      self.reset_fig()
      self.plot_solution()
      self.plot_polyhedrons()

  def show(self):
    plt.show()

  def plot_solution(self):
    last_pos = self.points[-1]

    self.ax.plot(*last_pos, linestyle="none",
                 marker='o', alpha=0.6,
                 markersize=10, markerfacecolor='black')

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
      #self.ax.plot(*coords, linestyle='+', color=c, alpha=1.0)
      self.ax.triplot(*coords, triangles=tri, color=c, alpha=a)

  def reset_fig(self):
    self.fig = plt.figure()
    self.ax = self.fig.add_subplot('111', aspect='equal', projection='3d')

    tup = [-2, 2]

    self.ax.set_xlim(tup)
    self.ax.set_ylim(tup)
    self.ax.set_zlim(tup)

    self.ax.elev = 30.

  def iterBound(self, nr_edges_init, error_init, prec):
    return None

  def polyhedron(self):
    """Return the inner polyhedron as a set of vertices"""
    p = self.convexify_polyhedron(self.inner)
    return p

  def outer_polyhedron(self):
    """Return the outer polyhedron as a set of vertices"""
    return self.convexify_polyhedron(self.outer)

  def save_polyhedron(self, fname):
    """Save the inner polyhedron as a set of vertices
       :param fname: Filename to which the polyhedron is saved
       :type fname: string"""
    np.savetxt(fname, self.polyhedron())

  def save_outer(self, fname):
    """Save the outer polyhedron as a set of vertices
       :param fname: Filename to which the polyhedron is saved
       :type fname: string"""
    np.savetxt(fname, self.convexify_polyhedron(self.outer))
