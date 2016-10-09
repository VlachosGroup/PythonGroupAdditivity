"""
===============================================================
Thermochemical property correlations (:mod:`chemtk.thermochem`)
===============================================================

This module provides a variety of different correlations used to represent
thermochemical properties (in non-dimensional form), i.e. heat capacity
|eq_ND_Cp_T|, enthalpy of formation |eq_ND_H_T|, and absolute |eq_ND_S_T| all
as functions of temperature |eq_T|.

This module also provides a group-additivity *property set* for estimation of
thermochemical properties.  Upon importing this module, the ``thermochem``
property set is registered for use within
:class:`chemtk.groupadd.GroupLibrary` instances.  Group contributions may also
be computed by fitting to thermochemical data via
:class:`ThermochemGroupFittingSet()`.objects.


-------
Summary
-------

.. autosummary::

    ThermochemBase
    ThermochemConstCp
    ThermochemNASA
    ThermochemNIST
    ThermochemPiecewise
    ThermochemRawData
    ThermochemIncomplete
    ThermochemGroup
    ThermochemGroupAdditive
    ThermochemGroupFittingSet



--------
Examples
--------

>>> from chemtk.thermochem import ThermochemNIST
>>> from chemtk.consts import GAS_CONSTANT as R
>>> from chemtk.units import with_units
>>> T_ref = with_units(298.15, 'K')
>>> # Correlation for formaldehyde direct from NIST webbook.
... 
>>> tc = ThermochemNIST([
...     5.193767, 93.23249, -44.85457, 7.882279, 0.551175, -119.3591, 202.4663,
...     -115.8972], range=(298.0, 1200.0))
>>> print (R*T_ref*tc.eval_ND_H(298.15)).fmt_in_units('kcal/mol')
-27.6998 kcal/mol
>>> print (R*tc.eval_ND_S(298.15)).fmt_in_units('cal/mol/K')
52.3313 cal/mol/K
>>> print (R*tc.eval_ND_Cp(400.0)).fmt_in_units('cal/mol/K')
9.38321 cal/mol/K
>>> print (R*tc.eval_ND_Cp(600.0)).fmt_in_units('cal/mol/K')
11.5247 cal/mol/K
>>> 


---------
Reference
---------

.. autodoc puts stuff extracted from docstrings here.
"""

# TODO: Test ThermochemQuartic code and add to correlations.
from .base import *
from .const_cp import *
from .group_data import *

__all__ = sum((getattr(globals()[module], '__all__') for module in [
    'base', 'const_cp', 'group_data']), [])

del base, const_cp, group_data

