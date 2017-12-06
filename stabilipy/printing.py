#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
