import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.integrate import quad as integrate
from .. import yaml_io
from .base import ThermochemBase
from ..Units import eval_qty

# from ..Consts import GAS_CONSTANT as R
import pmutt as pmutt
R = pmutt.constants.R(units='J/mol/K')
# R is now sourced from pmutt.constants instead of Consts


class ThermochemRawData(ThermochemBase):
    """
    Implement a thermochemical property correlation from raw data.

    Evaluated quantities are interpolated using a B-spline as discussed in
    :ref:`correlations documentation <correlations>`.
    """
    def __init__(self, ND_H_ref, ND_S_ref, Ts, ND_Cps,
                 T_ref=pmutt.constants.T0(units='K'), range=None):
        """
        Initialize a thermochemical property correlation from raw data.

        Parameters
        ----------
        ND_H_ref : float
            Non-dimensional standard heat of formation |eq_ND_H_ref|
        ND_S_ref : float
            Non-dimensional standard state reference entropy |eq_ND_S_ref|
        Ts : float array
            Temperatures at which `ND_Cps` are evaluated.
        ND_Cps : float
            Non-dimensional standard state heat capacities |eq_ND_Cp_T|
            evaluated at each temperature in `Ts`.
        T_ref : float, optional
            Reference temperature |eq_Tref| for `ND_H_ref` and `ND_S_ref`
            (default: room temperature according to pmutt, likely 298.15K).
        range : tuple(float, float), optional
            ``(lb, ub) = range`` where lb and ub are respectively the lower and
             uppers bounds of temperatures [K] for which the correlation is
             valid.    If specified, this range must contain T_ref and all data
             points in ND_Cp.
        """
        (self.Ts, self.ND_Cps) = list(zip(*sorted(
            zip(Ts, ND_Cps), key=lambda T_ND_Cps: T_ND_Cps[0])))
        self.min_T = Ts[0]
        self.max_T = Ts[-1]

        if range is None:
            range = (self.min_T, self.max_T)
        else:
            if self.min_T < range[0] or self.max_T > range[1]:
                raise ValueError(
                    'Heat capacity data points %g or %g lie outside of range'
                    ' [%g,%g].' % (self.min_T, self.max_T, range[0], range[1]))
        if T_ref < range[0] or T_ref > range[1]:
            raise ValueError(
                'T_ref=%g is outside the valid correlation range [%g,%g].'
                % (T_ref, range[0], range[1]))

        ThermochemBase.__init__(self, range)

        self.min_ND_Cp = self.ND_Cps[0]
        self.max_ND_Cp = self.ND_Cps[-1]

        self.ND_H_ref = ND_H_ref
        self.ND_S_ref = ND_S_ref
        self.T_ref = T_ref

        N = len(self.Ts)
        if N == 1:
            self.spline = ConstantSpline(self.ND_Cps[0])
        else:
            self.spline = InterpolatedUnivariateSpline(
                self.Ts, self.ND_Cps, k=(3 if N > 3 else N - 1))

    def get_CpoR(self, T):
        """Return non-dimensional standard state heat capacity |eq_ND_Cp_T|."""
        self.check_range(T)
        if not np.isscalar(T):
            return self._get_CpoR_ar(T)

        if T < self.min_T:
            return self.min_ND_Cp
        if T > self.max_T:
            return self.max_ND_Cp

        # Work-around for SciPy bug (?):
        # return self.spline(T)
        return float(self.spline(T))

    def _get_CpoR_ar(self, T):
        ND_Cp = np.empty(T.shape)
        T_below = T < self.min_T
        T_above = T > self.max_T
        T_middle = np.logical_not(np.logical_or(T_below, T_above))

        ND_Cp[T_below] = self.min_ND_Cp
        ND_Cp[T_above] = self.max_ND_Cp
        ND_Cp[T_middle] = self.spline(T[T_middle])

        return ND_Cp

    def get_SoR(self, T):
        """Return non-dimensional standard state entropy |eq_ND_S_T|."""
        self.check_range(T)
        T_a = self.T_ref
        T_b = T
        min_T = self.min_T
        max_T = self.max_T

        ND_S = self.ND_S_ref

        if T_a <= min_T:
            if T_b <= min_T:
                return ND_S + self.min_ND_Cp*np.log(T_b/T_a)
            ND_S += self.min_ND_Cp*np.log(min_T/T_a)
            T_a = min_T
        elif T_b <= min_T:
            ND_S += self.min_ND_Cp*np.log(T_b/min_T)
            T_b = min_T

        if T_a >= max_T:
            if T_b >= max_T:
                return ND_S + self.max_ND_Cp*np.log(T_b/T_a)
            ND_S += self.max_ND_Cp*np.log(max_T/T_a)
            T_a = max_T
        elif T_b >= max_T:
            ND_S += self.max_ND_Cp*np.log(T_b/max_T)
            T_b = max_T

        # The easiest, albeit not necessarily the best thing to do here is to
        # use numerical integration, so that's what we do.
        return ND_S + integrate(lambda t: self.spline(t)/t, T_a, T_b)[0]

    def get_HoRT(self, T):
        """Return non-dimensional standard heat of formation |eq_ND_H_T|."""
        self.check_range(T)
        T_a = self.T_ref
        T_b = T
        min_T = self.min_T
        max_T = self.max_T

        # This value represents the accumulated H/R (has temperature units).
        rH = self.ND_H_ref*T_a

        if T_a <= min_T:
            if T_b <= min_T:
                return (rH + self.min_ND_Cp*(T_b - T_a))/T
            rH += self.min_ND_Cp*(min_T - T_a)
            T_a = min_T
        elif T_b <= min_T:
            rH += self.min_ND_Cp*(T_b - min_T)
            T_b = min_T

        if T_a >= max_T:
            if T_b >= max_T:
                return rH + self.max_ND_Cp*(T_b - T_a)/T
            rH += self.max_ND_Cp*(max_T - T_a)
            T_a = max_T
        elif T_b >= max_T:
            rH += self.max_ND_Cp*(T_b - max_T)
            T_b = max_T

        return (rH + self.spline.integral(T_a, T_b))/T

    @classmethod
    def yaml_construct(cls, params, context):
        if 'T_ref' in params:
            T_ref = params['T_ref']
        else:
            #T_ref = eval_qty('298.15 K')  #fixed from eval_qty(298.15, 'K')
            T_ref = pmutt.constants.T0(units='K')
            #replaced getting room temp (298K) from eval_qty to pmutt.constants
        if 'ND_H_ref' in params:
            ND_H_ref = params['ND_H_ref']
        else:
            ND_H_ref = params['H_ref']/(R*T_ref)
        if 'ND_S_ref' in params:
            ND_S_ref = params['ND_S_ref']
        else:
            ND_S_ref = params['S_ref']/R

        if 'ND_Cp_data' in params:
            T_data, ND_Cp_data = list(zip(*params['ND_Cp_data']))
            Ts = np.array([T.in_units('K') for T in T_data])
            ND_Cps = np.array(ND_Cp_data)
        else:
            T_data, Cp_data = list(zip(*params['Cp_data']))
            Ts = np.array([T.in_units('K') for T in T_data])
            ND_Cps = np.array(
                [Cp for Cp in Cp_data])/R

        range = params.get('range')
        if range is not None:
            range = range[0].in_units('K'), range[1].in_units('K')
        else:
            range = Ts.min(), Ts.max()

        return cls(ND_H_ref, ND_S_ref, Ts, ND_Cps, T_ref.in_units('K'), range)

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

ND_Cp_data:
  type: list
  item_type:
     type: tuple
     item_types:
       - type: qty
         kind: temperature
       - type: float
  alts: Cp_data
  desc: set of (T, Cp(T)/R) data point pairs, where R is the gas constant

ND_H_ref:
  type: float
  alts: H_ref
  desc: H_ref/(R*T_ref) where R is the gas constant

ND_S_ref:
  type: float
  alts: ND_S_ref
  desc: S_ref/R where R is the gas constant

Cp_data:
  type: list
  item_type:
     type: tuple
     item_types:
       - type: qty
         kind: temperature
       - type: qty
         kind: molar heat capacity
  alts: ND_Cp_data
  desc: set of (T, Cp(T)) data point pairs

H_ref:
  type: qty
  kind: molar enthalpy
  alts: ND_H_ref
  desc: reference enthalpy

S_ref:
  type: qty
  kind: molar entropy
  alts: ND_S_ref
  desc: reference entropy
"""


yaml_io.register_class('ThermochemRawData',
                       yaml_io.parse(ThermochemRawData._yaml_schema),
                       ThermochemRawData)


class ConstantSpline(object):
    # This class emulates the interface to UnivariateSpline in the case of a
    # single data point (no interpolation).
    def __init__(self, ND_Cp):
        self.ND_Cp = ND_Cp

    def __call__(self, Ts):
        return self.ND_Cp*np.ones_like(Ts)

    def integral(self, T_a, T_b):
        return self.ND_Cp*(T_b - T_a)


__all__ = ['ThermochemRawData']
