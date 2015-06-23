import stability as stab
import numpy as np

execfile('box_contacts.py')

def main():
  global pos, normals

  #bar = sum(pos)/float(len(pos))
  #pos = [p-bar for p in pos]

  mu = 1.0
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  poly = stab.StabilityPolygon(60)
  poly.contacts = contacts

  poly.reset_fig()
  poly.plot_contacts()
  poly.show()

  #poly.make_problem()
  #poly.check_sizes()

  poly.compute(1e-2, True, False, False, True)
  return poly
