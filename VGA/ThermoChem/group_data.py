import numpy as np

from .. consts import GAS_CONSTANT as R
from .. ThermoChem import ThermochemBase

from ..GroupAdd.Library import GroupLibrary


__all__ = ['ThermochemGroupAdditive', 'ThermochemGroupFittingSet']


class ThermochemGroupAdditive(ThermochemBase):
    """
    Implement thermochemical property correlation using group-additivity.

    Properties are evaluated by taking a linear combination of contributing
    properties of the constituent groups.

    Normally this class is not instantiated directly; rather, instances are
    created and returned by :meth:`GroupLibrary.estimate()`.
    """
    def __init__(self, lib, groups):
        """
        Initialize thermochemical property correlation using group-additivity.

        Parameters
        ----------
        lib : :class:`GroupLibrary`
            Specify library from which to access group data.
        groups : mapping
            Map from :class:`Group` to int specifying counts of each group in
            the chemical structure.
        """
        self.correlations = []
        common_min = None
        common_max = None

        for group in groups:
            count = groups[group]
            correlation = lib[group]['thermochem']
            self.correlations.append((correlation, count))
            data_range = correlation.get_range()
            if data_range is not None:
                if common_min is None:
                    (common_min, common_max) = data_range
                else:
                    common_min = max(common_min, data_range[0])
                    common_max = min(common_max, data_range[1])

        if common_min is not None:
            ThermochemBase.__init__(self, range=(common_min, common_max))
        else:
            ThermochemBase.__init__(self, range=None)

    def eval_ND_Cp(self, T):
        return sum((count*correlation.eval_ND_Cp(T)
            for (correlation,count) in self.correlations))
    eval_ND_Cp.__doc__ = ThermochemBase.eval_ND_Cp.__doc__

    def eval_ND_H(self, T):
        return sum((count*correlation.eval_ND_H(T)
            for (correlation,count) in self.correlations))
    eval_ND_H.__doc__ = ThermochemBase.eval_ND_H.__doc__

    def eval_ND_S(self, T):
        return sum((count*correlation.eval_ND_S(T)
            for (correlation,count) in self.correlations))
    eval_ND_S.__doc__ = ThermochemBase.eval_ND_S.__doc__



GroupLibrary.register_property_set_type(
    'thermochem', ThermochemGroupAdditive)
