"""
====================================
Thermochemical property correlations
====================================

This module provides a variety of different correlations used to represent
thermochemical properties (in non-dimensional form), i.e. heat capacity
|eq_ND_Cp_T|, enthalpy of formation |eq_ND_H_T|, and absolute |eq_ND_S_T| all
as functions of temperature |eq_T|.

This module provides a group-additivity *property set* for estimation of
thermochemical properties.  Upon importing this module, the ``thermochem``
property set is registered for use within
:class:`chemtk.groupadd.GroupLibrary` instances.  Group contributions may also
be computed by fitting to thermochemical data via
:class:`ThermochemGroupFittingSet()`.objects.

---------
Reference
---------

.. autodoc puts stuff extracted from docstrings here.
"""

# TODO: Test ThermochemQuartic code and add to correlations.
from .base import *
from .raw_data import *
from .incomplete import *
from .group_data import *

__all__ = sum((getattr(globals()[module], '__all__') for module in [
    'base', 'raw_data', 'incomplete', 'group_data']), [])

del base, raw_data, incomplete, group_data
