import stabilipy as stab
import numpy as np

if __name__ == '__main__':

  A = np.vstack((np.eye(6), -np.eye(6)))
  b = np.ones(12,)

  linear_proj = stab.LinearProjection(3, A, b, None, None)
  linear_proj.compute(stab.Mode.precision, solver='cdd', epsilon=1e-3)
