#!/usr/bin/env python
# -*- coding: utf-8 -*-

import stability as stab
import numpy as np
import cdd
from plain_polygon import DoublePolygon
import shapely.geometry as geom
import copy

def main():
  normals = map(stab.normalize,
                [np.array([[0.5, 0.0, 1]]).T,
                 np.array([[0.5, 0.0, 1]]).T,
                 np.array([[0.0, 0.0, 1]]).T])

  pos = [np.array([[0.1, -0.2, 0]]).T,
         np.array([[-0.1, -0.2, 0]]).T,
         np.array([[0.2, 0, 1]]).T]

  mu = 0.5
  contacts = [stab.Contact(mu, p, n) for p, n in zip(pos, normals)]

  poly = stab.StabilityPolygon(60, dimension=2)
  poly.contacts = contacts

  poly.select_solver('cdd')
  poly.make_problem()
  poly.init_algo()
  poly.build_polys()

  dbpoly = DoublePolygon(poly.inner, poly.outer)

  print dbpoly.outer.edges

  print "Corresp:"
  for i, j in enumerate(dbpoly.correspondences):
    print "{} --> ({}, {})".format(i, *dbpoly.outer.edges[j][1:])
  for e in dbpoly.inner.edges:
    print "({}, {})".format(*e[1:])

  print "Vertex"
  print "-dbpoly:"
  print dbpoly.inner.vertices

  print "-dineq"
  dineq = [l/abs(l[0]) for l in dbpoly.inner.inequalities]
  print np.vstack(dineq)
  ineq = cdd.Polyhedron(poly.inner).get_inequalities()
  ineq.canonicalize()
  print "-ineq"
  print ineq

  poly.compute(stab.Mode.iteration, maxIter=2, plot_final=False)

  print "Adding..."
  print poly.points[-2]
  #dbpoly.addPointToInner(poly.points[-2], 3)

  print "Vertices"
  print "-dbpoly"
  print dbpoly.inner.vertices
  print "-inner"
  print poly.inner

  print "Ineq"

  dineq = [l/abs(l[0]) for l in dbpoly.inner.inequalities]
  print "-dineq"
  print np.vstack(dineq)

  vertices = cdd.Matrix(np.hstack([np.ones((dbpoly.inner.vertices.shape[0], 1)),
                                   dbpoly.inner.vertices]))
  vertices.rep_type = cdd.RepType.GENERATOR
  ineq = cdd.Polyhedron(vertices).get_inequalities()
  ineq.canonicalize()
  print "-ineq from vertex"
  print ineq

  spoly = geom.Polygon([(l[1], l[2]) for l in vertices])
  print spoly.exterior.is_ccw
  print spoly.convex_hull.exterior.coords.xy
  spoly.exterior.coords = list(spoly.exterior.coords)[::-1]
  print spoly.exterior.is_ccw

  i1 = dbpoly.inner.edges[3][1]
  i2 = dbpoly.inner.edges[3][2]
  print i1, i2

  e1 = dbpoly.outer.edges[dbpoly.correspondences[i1]]
  e2 = dbpoly.outer.edges[dbpoly.correspondences[i2]]
  print e1, e2
  old_edges = copy.deepcopy(dbpoly.outer.edges)
  print dbpoly.outer.vertices
  dbpoly.addInequalityToOuter(poly.offsets[-2], poly.directions[-2], i1, i2)

  print poly.outer
  print dbpoly.outer.inequalities

  print cdd.Polyhedron(poly.outer).get_generators()
  print dbpoly.outer.vertices
  print old_edges
  print dbpoly.outer.edges

  poly.plot()
  x, y = zip(*dbpoly.inner.vertices)
  poly.ax.plot(x, y)
  poly.show()


if __name__ == '__main__':
  main()
