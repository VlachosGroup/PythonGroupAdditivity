from .. ThermoChem import ThermochemBase, ThermochemIncomplete
from .. import yaml_io
import numpy as np
from ..GroupAdd.Library import GroupLibrary
from scipy import stats

__all__ = ['ThermochemGroup', 'ThermochemGroupAdditive']


class ThermochemGroup(ThermochemIncomplete):
    """
    Implement thermochemical property correlation for one group's contributions.

    See base class :class:`ThermochemIncomplete` for additional documentation.
    """
yaml_io.register_class('ThermochemGroup',
    yaml_io.parse(ThermochemGroup._yaml_schema), ThermochemGroup)


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
            xp = np.zeros((len(lib.uq_contents['descriptors']),1))
            for group in groups:
                count = groups[group]
                i = lib.uq_contents['descriptors'].index(group)
                xp[i] = count
            self.Xp_invXX_Xp = np.dot(np.dot(np.transpose(xp),lib.uq_contents['mat']),xp)
            self.dof = lib.uq_contents['dof']
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
    
    def eval_ND_Cp_CI(self, T, pt): #pt is percentile
        se = float(np.sqrt(np.square(self.RMSE.eval_ND_Cp(T))*self.Xp_invXX_Xp))
        return stats.t.ppf(1-(100.0-pt)/200.0,self.dof)*se
            
    def eval_ND_H_CI(self, T, pt):
        se = float(np.sqrt(np.square(self.RMSE.eval_ND_H(T))*self.Xp_invXX_Xp))
        return stats.t.ppf(1-(100.0-pt)/200.0,self.dof)*se
    
    def eval_ND_S_CI(self, T, pt):
        se = float(np.sqrt(np.square(self.RMSE.eval_ND_S(T))*self.Xp_invXX_Xp))
        return stats.t.ppf(1-(100.0-pt)/200.0,self.dof)*se



GroupLibrary.register_property_set_type(
    'thermochem', 'ThermochemGroup', ThermochemGroupAdditive)
