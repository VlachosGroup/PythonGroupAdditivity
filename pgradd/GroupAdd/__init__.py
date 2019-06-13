# 2016. Vlachos Group Geun Ho Gu. University of Delaware.
"""
=============================
Group-additivity computations
=============================

This module is used for using Group additivity method to compute energies

--------
Examples
--------

>>> from pgradd.GroupAdd.Library import GroupLibrary
>>> import pgradd.ThermoChem
>>> lib = GroupLibrary.Load('benson')
>>> groups = lib.GetDescriptors('CC')
>>> print groups
>>> thermochem = lib.Estimate(groups,'thermochem')
>>> print thermochem.get_HoRT(298.15)
defaultdict(<type 'int'>, {'C(C)(H)3': 2})
-34.4280812417

"""

from . Library import GroupLibrary
from . Scheme import GroupAdditivityScheme

__all__ = ['GroupLibrary', 'GroupAdditivityScheme']
