import numpy as np

def contacts(n):
  p = [np.array([[np.cos(i*2*np.pi/n)],
                 [np.sin(i*2*np.pi/n)],
                 [0]]) for i in range(n)]

  n = [np.array([[0],
                 [0],
                 [1]])]*n
  return p, n
