import abc

import numpy as np

from ..Units import eval_qty
from ..Error import OutsideCorrelationError



class ThermochemBase(object, metaclass=abc.ABCMeta):
    """Manage thermochemical properties correlation (abstract base class).

    This class cannot be instantiated directly; rather, it serves as a
    common base for the rest of the correlations.
    """
    # Abstract base class

    def init_params(self, range=None):
        if range is not None:
            self.range = (
                eval_qty(range[0], 'K').in_units('K'),
                eval_qty(range[1], 'K').in_units('K'))
        else:
            self.range = None

    def __init__(self, range=None):
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
                %(type(self).__name__, self.range[0], self.range[1], T))

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
    def eval_ND_Cp(self, T):
        """Return non-dimensional standard state heat capacity |eq_ND_Cp_T|."""

    @abc.abstractmethod
    def eval_ND_S(self, T):
        """Return non-dimensional standard state entropy |eq_ND_S_T|."""

    @abc.abstractmethod
    def eval_ND_H(self, T):
        """Return non-dimensional standard heat of formation |eq_ND_H_T|."""

    def eval_ND_G(self, T):
        """Return non-dimensional standard Gibbs energy of formation
        |eq_ND_G_T|."""
        return self.eval_ND_H(T) - self.eval_ND_S(T)

    _yaml_schema = """
range:
  type: tuple
  item_types: [{type: qty, kind: temperature},
               {type: qty, kind: temperature}]
  optional: True
  desc: range of valid temperatures
"""


from .. import yaml_io
yaml_io.register_class('ThermochemBase',
                       yaml_io.parse(ThermochemBase._yaml_schema),
                       ThermochemBase)


__all__ = ['ThermochemBase']
