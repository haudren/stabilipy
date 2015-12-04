import stability as stab
import numpy as np

execfile('stairs_contacts.py')

def main():
  global pos, normals

  #bar = sum(pos)/float(len(pos))
  #pos = [p-bar for p in pos]

  mu = 1.0
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  poly = stab.StabilityPolygon(60, dimension=2)
  poly.contacts = contacts

  poly.reset_fig()
  poly.plot_contacts()
  poly.show()

  #poly.make_problem()
  #poly.check_sizes()

  point = np.array([[0.36510907, 0.31419711, 0.73607441]]).T

  poly.addTorqueConstraint(contacts[-4:-2],
                           point,
                           10*np.ones((3, 1)))

  poly.addTorqueConstraint(contacts[-2:],
                           point,
                           10*np.ones((3, 1)))

  bar1 = sum([c.r for c in poly.contacts[0:4]])/4
  bar2 = sum([c.r for c in poly.contacts[4:8]])/4
  poly.addDistConstraint(bar1, 1.5)
  poly.addDistConstraint(bar2, 1.5)

  poly.make_problem()

  poly.compute(stab.Mode.best, maxIter=100, epsilon=2e-3,
               plot_init=True, plot_final=True)

  return poly

if __name__ == '__main__':
  main()
