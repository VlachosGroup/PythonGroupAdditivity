import numpy as np

from .base import ThermochemBase

from ..consts import GAS_CONSTANT as R
from ..units import *


class ThermochemConstCp(ThermochemBase):
    """
    Implement constant heat capacity thermochemical property correlation.

    See :ref:`correlations documentation <correlations>` for a precise
    definition.
    """
    def __init__(self, ND_H_ref, ND_S_ref, ND_Cp, T_ref=298.15, range=None):
        """
        Initialize constant heat capacity thermochemical property correlation.

        Parameters
        ----------
        ND_H_ref : float
            Non-dimensional standard heat of formation |eq_ND_H_ref|
        ND_S_ref : float
            Non-dimensional standard state reference entropy |eq_ND_S_ref|
        ND_Cp : float
            Non-dimensional standard state heat capacity |eq_ND_Cp_const|
        T_ref : float, optional
            Reference temperature |eq_Tref| for `ND_H_ref` and `ND_S_ref`
            (default: 298.15K).
        range : tuple(float, float), optional
            ``(lb, ub) = range`` where lb and ub are respectively the lower and
             uppers bounds of temperatures [K] for which the correlation is
             valid.
        """
        ThermochemBase.__init__(self, range)
        self.ND_Cp = ND_Cp
        self.ND_H_ref = ND_H_ref
        self.ND_S_ref = ND_S_ref
        self.T_ref = T_ref

    def eval_ND_Cp(self, T):
        """Return non-dimensional standard state heat capacity |eq_ND_Cp_T|."""
        self.check_range(T)
        if not np.isscalar(T):
            return self.ND_Cp*np.ones(T.shape)
        else:
            return self.ND_Cp

    def eval_ND_S(self, T):
        """Return non-dimensional standard state entropy |eq_ND_S_T|."""
        self.check_range(T)
        return self.ND_S_ref + self.ND_Cp*np.log(T/self.T_ref)

    def eval_ND_H(self, T):
        """Return non-dimensional standard heat of formation |eq_ND_H_T|."""
        self.check_range(T)
        return (self.ND_H_ref*self.T_ref + self.ND_Cp*(T - self.T_ref))/T

    @classmethod
    def yaml_construct(cls, params, context):
        if 'ND_Cp' in params:
            ND_Cp = params['ND_Cp']
        else:
            ND_Cp = params['Cp']/R
        T_ref = params['T_ref']
        if 'ND_H_ref' in params:
            ND_H_ref = params['ND_H_ref']
        else:
            ND_H_ref = params['H_ref']/(R*T_ref)
        if 'ND_S_ref' in params:
            ND_S_ref = params['ND_S_ref']
        else:
            ND_S_ref = params['S_ref']/R

        range = params.get('range')
        if range is not None:
            range = range[0].in_units('K'), range[1].in_units('K')
        return cls(ND_H_ref, ND_S_ref, ND_Cp, T_ref.in_units('K'), range)

    _yaml_schema = """
range:
  type: tuple
  item_types: [{type: qty, kind: temperature},
               {type: qty, kind: temperature}]
  optional: True
  desc: range of valid temperatures

T_ref:
  type: qty
  kind: temperature
  default: 298.15 K
  desc: reference temperature for reference enthalpy and entropy

Cp:
  type: qty
  kind: molar heat capacity
  alts: ND_Cp
  desc: constant heat capacity

H_ref:
  type: qty
  kind: molar enthalpy
  alts: ND_H_ref
  desc: reference heat of formation

S_ref:
  type: qty
  kind: molar entropy
  alts: ND_S_ref
  desc: reference entropy

ND_Cp:
  type: float
  alts: Cp
  desc: Cp/R where R is the gas constant

ND_H_ref:
  type: float
  alts: H_ref
  desc: H_ref/(R*T_ref) where R is the gas constant

ND_S_ref:
  type: float
  alts: S_ref
  desc: S_ref/R where R is the gas constant
"""

__all__ = ['ThermochemConstCp']
