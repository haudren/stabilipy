import stabilipy as stab


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
