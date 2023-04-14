import abc

import numpy as np
from pmutt import constants as c

from ..Units import eval_qty
from ..Error import OutsideCorrelationError
from .. import yaml_io


class ThermochemBase(object):
    metaclass = abc.ABCMeta
    """Manage thermochemical properties correlation (abstract base class).

    This class cannot be instantiated directly; rather, it serves as a
    common base for the rest of the correlations.
    """
    # Abstract base class

    def init_params(self, range=None, units=None):
        if range is not None:
            self.range = (
                eval_qty(range[0], 'K').in_units('K'),
                eval_qty(range[1], 'K').in_units('K'))
        else:
            self.range = None

    def __init__(self, range=None, units=None):
        """
        Initialize generic thermochemical property correlation.

        Parameters
        ----------
        range : tuple(float, float), optional
            ``(lb, ub) = range`` where lb and ub are respectively the lower and
             uppers bounds of temperatures [K] for which the correlation is
             valid.
        """
        if range is not None:
            range = tuple(range)
            assert len(range) == 2
            assert range[1] >= range[0]
        self.range = range

    def check_range(self, T):
        """
        Verify that one or more temperatures is within the valid range for this
        correlation.

        Parameters
        ----------
        T : float or float array
            Specify temperature(s) to test.

        Returns
        -------
        None

        Raises
        ------
        chemtk.error.OutsideCorrelationError
            If `T` (or any temperature in array `T`) is outside the valid range
            of this correlation.
        """
        if self.range is None:
            return
        if np.any(T < self.range[0]) or np.any(T > self.range[1]):
            raise OutsideCorrelationError(
                "In %s().check_range: one or more of the following "
                "temperatures is outside the valid interval [%r,%r]: %r"
                % (type(self).__name__, self.range[0], self.range[1], T))

    def get_range(self):
        """
        Return range of validity for this thermochemical properties
        correlation.

        Returns
        -------
        lower_bound : float
            Lower bound of range.
        upper_bound : float
            Upper bound of range.
        """
        return self.range

    def set_range(self, range=None):
        """
        Set range of validity for this thermochemical properties correlation.

        Parameters
        ----------
        range : tuple(float, float), optional
            ``(lb, ub) = range`` where lb and ub are respectively the lower and
             uppers bounds of temperatures [K] for which the correlation is
             valid.
        """
        self.range = range

    @abc.abstractmethod
    def get_CpoR(self, T):
        """Return non-dimensional standard state heat capacity |eq_ND_Cp_T|."""

    @abc.abstractmethod
    def get_SoR(self, T, S_elements=None):
        """Return non-dimensional standard state entropy |eq_ND_S_T|."""

    @abc.abstractmethod
    def get_HoRT(self, T):
        """Return non-dimensional standard heat of formation |eq_ND_H_T|."""

    def get_GoRT(self, T, S_elements=None):
        """
        Return non-dimensional standard Gibbs energy of formation |eq_ND_G_T|.
        """
        return self.get_HoRT(T) - self.get_SoR(T, S_elements=S_elements)

    def get_H(self, T, units):
        """
        Return dimensional enthalpy at T with units specified by 'units'
        """
        return self.get_HoRT(T)*T*c.R('{}/K'.format(units))

    def get_G(self, T, units, S_elements=None):
        """
        Return dimensional Gibbs free energy at T with units
        specified by 'units'
        """
        return self.get_GoRT(T, S_elements=S_elements)*T*c.R('{}/K'.format(units))

    def get_S(self, T, units, S_elements=None):
        """
        Return dimensional entropy with units specified by 'units'
        """
        return self.get_SoR(T, S_elements=S_elements)*c.R(units)

    def get_Cp(self, T, units):
        """
        Return dimensional heat capacity with units specified by 'units'
        """
        return self.get_CpoR(T)*c.R(units)

    _yaml_schema = """
range:
  type: tuple
  item_types: [{type: qty, kind: temperature},
               {type: qty, kind: temperature}]
  optional: True
  desc: range of valid temperatures
"""


yaml_io.register_class('ThermochemBase',
                       yaml_io.parse(ThermochemBase._yaml_schema),
                       ThermochemBase)


__all__ = ['ThermochemBase']
