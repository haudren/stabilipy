#!/usr/bin/env python2
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

from __future__ import print_function
from builtins import object
from enum import IntEnum, unique

@unique
class Verbosity(IntEnum):

  """Enum representing verbosity levels"""
  none = 0
  error = 1
  info = 2
  debug = 3

class Printer(object):

  """Docstring for Printer. """

  def __init__(self, verbosity):
    """Initialize with verbosity level """
    self.verbosity = verbosity

  def __call__(self, text, verbosity):
    """Print only if verbosity level is higher
       than current verbosity"""
    if verbosity <= self.verbosity:
      print(text)
