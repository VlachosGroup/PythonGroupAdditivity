from .. ThermoChem import ThermochemBase, ThermochemIncomplete
from .. import yaml_io
import numpy as np
from ..GroupAdd.Library import GroupLibrary

__all__ = ['ThermochemGroup', 'ThermochemGroupAdditive']


class ThermochemGroup(ThermochemIncomplete):
    """
    Implement thermochemical property correlation for
    one group's contributions.

    See base class :class:`ThermochemIncomplete` for additional documentation.
    """


yaml_io.register_class('ThermochemGroup',
                       yaml_io.parse(ThermochemGroup._yaml_schema),
                       ThermochemGroup)


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
        if lib.uq_contents:
            self.RMSE = lib.uq_contents['RMSE'].thermochem
            xp = np.zeros((len(lib.uq_contents['descriptors']), 1))
            for group in groups:
                count = groups[group]
                i = lib.uq_contents['descriptors'].index(group)
                xp[i] = count
            self.Xp_invXX_Xp = np.dot(np.dot(np.transpose(xp),
                                             lib.uq_contents['mat']), xp)
            self.dof = lib.uq_contents['dof']
        if common_min is not None:
            ThermochemBase.__init__(self, range=(common_min, common_max))
        else:
            ThermochemBase.__init__(self, range=None)

    def get_CpoR(self, T):
        return sum((count*correlation.get_CpoR(T)
                    for (correlation, count) in self.correlations))
    get_CpoR.__doc__ = ThermochemBase.get_CpoR.__doc__

    def get_HoRT(self, T):
        return sum((count*correlation.get_HoRT(T)
                    for (correlation, count) in self.correlations))
    get_HoRT.__doc__ = ThermochemBase.get_HoRT.__doc__

    def get_SoR(self, T):
        return sum((count*correlation.get_SoR(T)
                    for (correlation, count) in self.correlations))
    get_SoR.__doc__ = ThermochemBase.get_SoR.__doc__

    def get_CpoR_SE(self, T):
        return float(np.sqrt(np.square(self.RMSE.get_CpoR(T)) *
                             self.Xp_invXX_Xp))

    def get_HoRT_SE(self, T):
        return float(np.sqrt(np.square(self.RMSE.get_HoRT(T)) *
                             self.Xp_invXX_Xp))

    def get_SoR_SE(self, T):
        return float(np.sqrt(np.square(self.RMSE.get_SoR(T)) *
                             self.Xp_invXX_Xp))


GroupLibrary.register_property_set_type(
    'thermochem', 'ThermochemGroup', ThermochemGroupAdditive)
