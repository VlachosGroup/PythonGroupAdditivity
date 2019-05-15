"""
=======================================
Physical quantities and unit conversion
=======================================

This module provides routines for representing and manipulating physical
quantities with units.  A physical quantity is represented internally by its
numerical value in SI units combined with a :class:`FundamentalUnits` that
specifies the exponents of meters, kilograms, seconds, ampere, kelvin,
gram-moles, and candela represented in the quantity.

All standard prefixes (*e.g.* k for kilo, u for micro, etc.) are supported as
are most derived units in the metric and standard systems.  Notable absences
include celsius and fahrenheit as these do not represent temperatures in the
absolute sense.  Proper handling of these units in arithmetic is complex and
beyond the intended scope of this library.


-------
Summary
-------

.. autosummary::

    eval_quantity
    eval_qty
    Quantity
    FundamentalUnits
    with_units
    in_units
    has_units
    from_SI_to
    to_SI_from


--------
Examples
--------

>>> from chemtk.units import with_units, in_units, has_units
>>> g = with_units(9.81, 'm/s^2')
>>> m = with_units(10.0, 'lb')
>>> print m*g
44.497411497 m*kg/s^(-2)
>>> print in_units(m, 'g')
4535.9237
>>> print (m*g).fmt_in_units('N')
44.4974 N
>>> print (m*g).fmt_in_units('lbf')
10.0034 lbf
>>> print m.has_units('ug')
True
>>> print m.has_units('J s')
False
>>>


---------
Reference
---------

.. autodoc puts stuff extracted from docstrings here.
"""

from .qty import *
from .helpers import *
from . import builtin

__all__ = qty.__all__[:] + helpers.__all__[:]

del qty, helpers, builtin
