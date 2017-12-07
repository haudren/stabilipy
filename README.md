StabiliPy
=========

[![Build Status](https://travis-ci.org/haudren/stabilipy.svg?branch=master)](https://travis-ci.org/haudren/stabilipy) [![Doc Status](https://readthedocs.org/projects/stabilipy/badge/)](https://stabilipy.readthedocs.io)

StabiliPy is a Python package that provides a simple interface to compute
static equilibrium polygons and robust static equilibrium polyhedrons.

This problem is a projection problem, and uses the recursive projection
method to compute those regions. You can thus use this package to compute
arbitrary convex projections by recursive projection.

The documentation is available at:
https://stabilipy.readthedocs.io

This package is distributed under the GNU Lesser General Public License version 3
or above (LGPLv3+).

Installation
============

In order to install all possible backends use:
> sudo apt-get install liblapack-dev libatlas-dev libblas-dev libgmp-dev libppl-dev

Followed by:
> pip install cython && pip install -r requirements.txt

Then:
> pip install .

Should do the trick !

Notes:
- You need LAPACK, Atlas and BLAS to install CVXOPT
- You need Cython to install [pycddlib](https://github.com/mcmtroffaes/pycddlib)
- You need Cython, gmp, and the Parma Polyhedra Library to install [pyparma](https://github.com/haudren/pyparma)
