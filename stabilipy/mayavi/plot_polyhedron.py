#!/usr/bin/env python2

from __future__ import print_function
from builtins import zip
import numpy as np
from mayavi import mlab
import os
import fnmatch

import argparse

from scipy.spatial import ConvexHull

def load_data(fname):
  return np.loadtxt(fname)

def plot_poly(mat, alpha):
  cols = [np.array([c]) for c in zip(*mat)]
  mlab.points3d(*cols, scale_factor=0.01)

  tri = ConvexHull(mat)

  surf = mlab.triangular_mesh(cols[0], cols[1], cols[2], tri.simplices,
                              opacity=alpha, colormap='copper')
  return surf

def compute_surface(mat):
  cols = [np.array([c]) for c in zip(*mat)]
  vtk_source = mlab.pipeline.scalar_scatter(*cols, figure=False)
  delaunay = mlab.pipeline.delaunay3d(vtk_source)
  return delaunay

def poly_from_file(directory, fname):
  return load_data(os.path.join(directory, fname))

def poly_generator(directory, inners, outers):
  for i, o in zip(inners, outers):
    A = poly_from_file(directory, i)
    B = poly_from_file(directory, o)
    yield A, B

@mlab.animate(delay=500, ui=True)
def animation_poly(directory, inner, outer, view):
  mlab.view(**view)
  for A, B in poly_generator(directory, inner, outer):
    mlab.clf()
    plot_poly(A, 1.)
    plot_poly(B, 0.2)
    yield

def record_frames(directory, inner, outer, out_dir, view):
  nr_iter = 1
  max_iter = 300

  polys = list(poly_generator(directory, inner, outer))
  nr_poly = len(polys)

  for A, B in polys:
    mlab.clf()
    mlab.view(**view)
    plot_poly(A, 1.)
    plot_poly(B, 0.2)
    print("{} / {}".format(nr_iter, nr_poly))
    fname = os.path.join(out_dir, 'polyhedrons_{0:04d}.png'.format(nr_iter))
    mlab.savefig(fname)
    nr_iter += 1
    if nr_iter > max_iter:
      break

def main(base_name, out_dir, animate):
  directory, filename = os.path.split(base_name)
  files = os.listdir(directory)

  inner_name = filename+'_inner_*'
  outer_name = filename+'_outer_*'

  inner = sorted(fnmatch.filter(files, inner_name))
  outer = sorted(fnmatch.filter(files, outer_name))

  assert len(inner) == len(outer),\
      "Not the same number of inner and outer polygons"

  view = {'azimuth': 160,
          'elevation': 20,
          'distance': 0.5}

  if animate:
    mlab.figure(size=(800, 700))
    print("Animating")
    a = animation_poly(directory, inner, outer, view)  # noqa
    mlab.show()
    print("Done animating")
  else:
    mlab.options.offscreen = True
    mlab.figure(size=(800, 700))
    record_frames(directory, inner, outer, out_dir, view)

if __name__ == '__main__':
  parser = argparse.ArgumentParser('plot-polyhedrons', description='Use this program to generate successive convex hulls from the output of 3D convex hulls')  # noqa
  parser.add_argument('base_name', help='This is the base name of your polyhedrons, typically stairs_poly or example poly. This will be expanded to base_name_{inner,outer}_{0:04d}')  # noqa
  parser.add_argument('out_dir', help='This is the folder where successive polyhedrons will be saved')  # noqa
  parser.add_argument('--save', action='store_false', help='Set it to render all frames offscreen')  # noqa

  args = parser.parse_args()

  main(args.base_name, args.out_dir, args.save)
