import stabilipy as stab


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

class SphereProjection(stab.RecursiveProjectionProblem):
  """Try to approximate a sphere of radius r"""

  def __init__(self, radius):
    """:param radius: Radius of the sphere
       :type radius: double"""
    stab.RecursiveProjectionProblem.__init__(self, dimension=3)
    self.radius = radius

  def solve(self, d):
    """We are computing a sphere so the extremal point in direction d is just r*d"""
    return self.radius*d

if __name__ == '__main__':
  sphere = SphereProjection(1.0)
  sphere.compute(solver='cdd', mode=stab.Mode.iteration, maxIter=50)
