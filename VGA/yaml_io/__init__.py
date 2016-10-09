"""
=======================
YAML-formatted data I/O
=======================

This module handles input from YAML formatted data.  Because this code is in
bad need of cleanup and because in the long run, I (Stephen Edie) hope to
replace the input mechanism with something cleaner and faster (not involving
YAML), the documentation for this module is limited at this time.


---------
Reference
---------

.. autodoc puts stuff extracted from docstrings here.
"""


from .yaml_io import *

__all__ = yaml_io.__all__[:]

del yaml_io
