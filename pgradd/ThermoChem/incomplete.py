from warnings import warn

import numpy as np

from .. Error import (
    IncompleteDataError, IncompleteDataWarning, ReadOnlyDataError,
    OutsideCorrelationError)
#from .. Consts import GAS_CONSTANT as R
import pmutt.constants as c
from .. Units import with_units
from .. Units import eval_qty
from .. ThermoChem import ThermochemBase, ThermochemRawData
from .. import yaml_io
R = eval_qty(str(c.R(units='J/mol/K')) + ' J/(mol K)')

__all__ = ['ThermochemIncomplete']


class ThermochemIncomplete(ThermochemBase):
    """
    Implement thermochemical property correlation for data missing
    heat capacity.

    Attempts to evaluate |eq_ND_H_T| or |eq_ND_S_T| will fail by raising
    :exc:`chemtk.error.IncompleteDataError` if `ND_H_ref` or `ND_S_ref`
    are respectively missing.

    If `ND_Cp_data` are missing, then attempts to evaluate |eq_ND_Cp_T|
    will also fail.  Attempts to evaluate |eq_ND_H_T| or |eq_ND_S_T| at
    temperatures other than |eq_Tref| will assume |eq_ND_Cp_T| is zero and
    will issue an :exc:`chemtk.error.IncompleteDataWarning`.

    If `ND_Cp_data` are present, quantities will be evaluated
    (where otherwise possible ) using a
    :class:`chemtk.thermochem.ThermochemRawData` correlation.
    """
    def __init__(self, ND_H_ref=None, ND_S_ref=None, ND_Cp_data={},
                 T_ref=298.15, range=None):
        """
        Initialize thermochemical property correlation for incomplete data.

        Parameters
        ----------
        ND_H_ref : float, optional
            Non-dimensional standard heat of formation |eq_ND_H_ref|
        ND_S_ref : float, optional
            Non-dimensional standard state reference entropy |eq_ND_S_ref|
        ND_Cp_data : mapping, optional
            Non-dimensional standard state heat capacity data |eq_ND_Cp_T|,
            provided as a mapping from temperature [K] to non-dimensional
            heat capacity.
        T_ref : float, optional
            Reference temperature |eq_Tref| for `ND_H_ref` and `ND_S_ref`
            (default: 298.15K).
        range : tuple(float, float), optional
            ``(lb, ub) = range`` where lb and ub are respectively the lower and
             uppers bounds of temperatures [K] for which the correlation is
             valid.
        """
        ThermochemBase.__init__(self, range=range)

        self.ND_H_ref = ND_H_ref
        self.ND_S_ref = ND_S_ref
        self.ND_Cp_data = ND_Cp_data.copy()
        self.T_ref = T_ref
        self._setup_correlation()

    def _expand_ND_Cp_data(self, ND_Cp_data):
        if not ND_Cp_data:
            return [], []
        else:
            return list(zip(*sorted(
                list(ND_Cp_data.items()), key=lambda item: item[0])))

    def _setup_correlation(self):
        if hasattr(self, '_correlation'):
            del self._correlation

        if self.ND_Cp_data:
            ND_H_ref = self.ND_H_ref if self.ND_H_ref is not None else 0.0
            ND_S_ref = self.ND_S_ref if self.ND_S_ref is not None else 0.0
            (Ts, ND_Cps) = self._expand_ND_Cp_data(self.ND_Cp_data)
            self._correlation = ThermochemRawData(
                ND_H_ref, ND_S_ref, Ts, ND_Cps, self.T_ref, self.get_range())

    def has_ND_Cp(self, T=None):
        """Return True if correlation has |eq_ND_Cp_T| data (possibly at `T`).

        Parameters
        ----------
        T : float, optional
            If `T` is specified, then check for a data point explicitly at `T`.
            Otherwise, check for the existence of any data points.
        """
        if T is None:
            return bool(self.ND_Cp_data)
        else:
            return T in self.ND_Cp_data

    def has_ND_H(self):
        """Return True if correlation has |eq_ND_H_T| data."""
        return bool(self.ND_H_ref)

    def has_ND_S(self):
        """Return True if correlation has |eq_ND_S_T| data."""
        return bool(self.ND_S_ref)

    def del_ND_Cp(self, T=None):
        """Delete |eq_ND_Cp_T| data (possibly at `T` only).

        Parameters
        ----------
        T : float, optional
            If `T` is specified, then delete data point at `T`.  Otherwise
            delete all |eq_ND_Cp_T| data.
        """
        if T is None:
            self.ND_Cp_data = None
        else:
            del self.ND_Cp_data[T]
        self._setup_correlation()

    def del_ND_H_ref(self):
        """Delete |eq_ND_H_ref| data."""
        self.ND_H_ref = None

    def del_ND_S_ref(self):
        """Delete |eq_ND_S_ref| data."""
        self.ND_S_ref = None

    def get_CpoR(self, T):
        if not self.ND_Cp_data:
            raise IncompleteDataError(
                "Cannot evaluate ND_Cp: no heat capacity data is available")
        else:
            try:
                return self._correlation.get_CpoR(T)
            except OutsideCorrelationError:
                raise IncompleteDataError(
                    "Cannot evaluate ND_Cp: no heat capacity data for T=%g"
                    % T)

    get_CpoR.__doc__ = ThermochemBase.get_CpoR.__doc__

    def get_HoRT(self, T):
        if self.ND_H_ref is None:
            raise IncompleteDataError(
                "Cannot evaluate ND_H: no enthalpy data is available")
        elif not self.ND_Cp_data:
            if T != self.T_ref:
                warn(
                    "Evaluation of ND_H_ref with (T=%g <=> T_ref=%g) will not"
                    " be corrected because heat capacity data is not"
                    " available." % (T, self.T_ref),
                    IncompleteDataWarning)
            return self.ND_H_ref
        else:
            try:
                return self._correlation.get_HoRT(T)
            except OutsideCorrelationError:
                raise IncompleteDataError(
                    "Cannot evaluate ND_H: no heat capacity data for T=%g"
                    % T)
    get_HoRT.__doc__ = ThermochemBase.get_HoRT.__doc__

    def get_SoR(self, T):
        if self.ND_S_ref is None:
            raise IncompleteDataError(
                "Cannot evaluate ND_S: no entropy data is available")
        elif not self.ND_Cp_data:
            if T != self.T_ref:
                warn(
                    "Evaluation of ND_S_ref with (T=%g <=> T_ref=%g) will not"
                    " be corrected because heat capacity data is not"
                    " available." % (T, self.T_ref),
                    IncompleteDataWarning)
            return self.ND_S_ref
        else:
            try:
                return self._correlation.get_SoR(T)
            except OutsideCorrelationError:
                raise IncompleteDataError(
                    "Cannot evaluate ND_S: no heat capacity data for T=%g"
                    % T)
    get_SoR.__doc__ = ThermochemBase.get_SoR.__doc__

    def copy(self):
        """Return a copy of this correlation.

        Returns
        -------
        thermochem : :class:`ThermochemIncomplete`
            Copy of this correlation.
        """
        return type(self)(
            self.ND_H_ref, self.ND_S_ref, self.ND_Cp_data, self.T_ref,
            self.range)

    def update(self, correlation, overwrite=False):
        """Add data contained in `correlation` into this correlation object.

        Parameters
        ----------
        correlation : :class:`ThermochemIncomplete` instance
            Specify correlation to import data from.
        overwrite : bool, optional
            True if existing data should be overwritten by data in
            `correlation`.

        Raises
        ------
        chemtk.error.ReadOnlyDataError
            If `correlation` contains data that differs from data already in
            this correlation and `overwrite` is False.
        """
        if not isinstance(correlation, type(self)):
            raise TypeError("In ThermochemIncomplete.update():",
                            " argument not instance of ThermochemIncomplete.")

        # Copy existing data first.
        data_range = self.get_range()
        T_ref = self.T_ref
        ND_H_ref = self.ND_H_ref
        ND_S_ref = self.ND_S_ref
        ND_Cp_data = self.ND_Cp_data.copy()

        # Compute updated range (as the union of the two).
        other_data_range = correlation.get_range()
        if other_data_range is not None:
            if data_range is None:
                data_range = other_data_range
            else:
                data_range = (
                    min(data_range[0], other_data_range[0]),
                    max(data_range[1], other_data_range[1]))

        # Compute updated ND_Cp_data.
        # (This is done first so the new data is available below.)
        if correlation.has_ND_Cp():
            other_ND_Cp_data = correlation.ND_Cp_data
            for T in other_ND_Cp_data:
                if(not overwrite
                   and T in self.ND_Cp_data
                   and other_ND_Cp_data[T] != ND_Cp_data[T]):
                    raise ReadOnlyDataError(
                        "In ThermochemIncomplete.update():"
                        " property ND_Cp(T=%g) already exists and differs from"
                        " new value." % T)
                ND_Cp_data[T] = other_ND_Cp_data[T]
            (Ts, ND_Cps) = self._expand_ND_Cp_data(ND_Cp_data)

        else:
            (Ts, ND_Cps) = self._expand_ND_Cp_data(ND_Cp_data)

        # Compute updated ND_H_ref and ND_S_ref.
        if correlation.has_ND_H() or correlation.has_ND_S():
            other_T_ref = correlation.T_ref
            other_ND_H_ref = correlation.ND_H_ref
            other_ND_S_ref = correlation.ND_S_ref

            # Create temporary incomplete correlation using:
            #    - Aggregated Cp_data and range
            #    - ND_H_ref, ND_S_ref and T_ref from other correlation
            # This is done to translate ND_H_ref and ND_S_ref between
            # potentially differing references temperatures.
            test_correlation = type(self)(
                other_ND_H_ref, other_ND_S_ref, ND_Cp_data, other_T_ref,
                data_range)

            # Update ND_H_ref.
            if correlation.has_ND_H():
                new_ND_H_ref = test_correlation.get_HoRT(T_ref)
                if(not overwrite
                   and ND_H_ref is not None
                   and new_ND_H_ref != ND_H_ref):
                    raise ReadOnlyDataError(
                        "In ThermochemIncomplete.update():"
                        " property ND_H_ref already exists and differs from"
                        " new value.")
                ND_H_ref = new_ND_H_ref

            # Update ND_S_ref.
            if correlation.has_ND_S():
                new_ND_S_ref = test_correlation.get_SoR(T_ref)
                if(not overwrite
                   and ND_S_ref is not None
                   and new_ND_S_ref != ND_S_ref):
                    raise ReadOnlyDataError(
                        "In ThermochemIncomplete.update():"
                        " property ND_S_ref already exists and differs from"
                        " new value.")
                ND_S_ref = new_ND_S_ref

        # Now store new data and update internal correlation.
        self.set_range(data_range)
        self.T_ref = T_ref
        self.ND_H_ref = ND_H_ref
        self.ND_S_ref = ND_S_ref
        self.ND_Cp_data = ND_Cp_data
        self._setup_correlation()

    @classmethod
    def yaml_construct(cls, params, context):
        T_ref = params['T_ref']
        if params.get('ND_H_ref') is not None:
            ND_H_ref = params['ND_H_ref']
        elif params.get('H_ref') is not None:
            ND_H_ref = params['H_ref']/(R*T_ref)
        else:
            ND_H_ref = None

        if params.get('ND_S_ref') is not None:
            ND_S_ref = params['ND_S_ref']
        elif params.get('S_ref') is not None:
            ND_S_ref = params['S_ref']/R
        else:
            ND_S_ref = None

        if params.get('ND_Cp_data'):
            (T_data, ND_Cp_data) = list(zip(*params['ND_Cp_data']))
            Ts = np.array([T.in_units('K') for T in T_data])
            ND_Cps = np.array(ND_Cp_data)
            ND_Cp_data = dict(list(zip(Ts, ND_Cps)))
        elif params.get('Cp_data'):
            (T_data, Cp_data) = list(zip(*params['Cp_data']))
            Ts = np.array([T.in_units('K') for T in T_data])
            ND_Cps = np.array([Cp/R for Cp in Cp_data])
            ND_Cp_data = dict(list(zip(Ts, ND_Cps)))
        else:
            ND_Cp_data = {}

        range = params.get('range')
        if range is not None:
            range = range[0].in_units('K'), range[1].in_units('K')
        T_ref = T_ref.in_units('K')

        return cls(ND_H_ref, ND_S_ref, ND_Cp_data, T_ref, range)

    def yaml_format(self, units={}):
        """Output YAML format definition of this object.

        Other Parameters
        ----------------
        units : dict
            Map containing one or more of the following keys:

                :``molar enthalpy``: units to use for enthalpy values
                :``molar entropy``: units to use for entropy values
                :``molar heat capacity``: units to use for heat capacity values
                :``temperature``: units to use for temperature values

            If any of ``molar enthalpy``, ``molar entropy``, or
            ``molar heat capacity`` is unspecified or None, then write these as
            dimensionless values.  If ``temperature`` is unspecified, then
            write temperatures in [K].
        """
        lines = []
        T_ref = with_units(self.T_ref, 'K')
        T_units = units.get('temperature', 'K')
        lines.append('T_ref: %s' % T_ref.fmt_in_units(T_units))

        if self.has_ND_H():
            H_units = units.get('molar enthalpy')
            if H_units:
                lines.append(
                    'H_ref: %s' % (R*T_ref *
                                   self.ND_H_ref).fmt_in_units(H_units))
            else:
                lines.append('ND_H_ref: %r' % self.ND_H_ref)

        if self.has_ND_S():
            S_units = units.get('molar entropy')
            if S_units:
                lines.append(
                    'S_ref: %s' % (R*self.ND_S_ref).fmt_in_units(S_units))
            else:
                lines.append('ND_S_ref: %r' % self.ND_S_ref)

        if self.has_ND_Cp():
            Cp_units = units.get('molar heat capacity')
            if Cp_units:
                lines.append('Cp_data:')
                for T in sorted(self.ND_Cp_data):
                    lines.append('    - [%s, %s]' % (
                        with_units(T, 'K').fmt_in_units(T_units),
                        (R*self.ND_Cp_data[T]).fmt_in_units(Cp_units)))
            else:
                lines.append('ND_Cp_data:')
                for T in sorted(self.ND_Cp_data):
                    lines.append('    - [%s, %r]' % (
                        with_units(T, 'K').fmt_in_units(T_units),
                        self.ND_Cp_data[T]))

        range = self.get_range()
        if range is not None:
            lines.append('range: [%s, %s]' % (
                with_units(range[0], 'K').fmt_in_units(T_units),
                with_units(range[1], 'K').fmt_in_units(T_units)))
        return '\n'.join(lines)

    _yaml_schema = """
range:
  type: tuple
  item_types: [{type: qty, kind: temperature},
               {type: qty, kind: temperature}]
  optional: true
  desc: range of valid temperatures

T_ref:
  type: qty
  kind: temperature
  default: 298.15 K
  desc: reference temperature for reference enthalpy and entropy

ND_Cp_data:
  optional: true
  type: list
  item_type:
     type: tuple
     item_types: [{type: qty, kind: temperature}, {type: float}]
  desc: set of (T, Cp(T)/R) data point pairs, where R is the gas constant

ND_H_ref:
  optional: true
  type: float
  desc: H_ref/(R*T_ref) where R is the gas constant

ND_S_ref:
  optional: true
  type: float
  desc: S_ref/R where R is the gas constant

Cp_data:
  optional: true
  type: list
  item_type:
     type: tuple
     item_types:
       - type: qty
         kind: temperature
       - type: qty
         kind: molar heat capacity
  desc: set of (T, Cp(T)) data point pairs

H_ref:
  optional: true
  type: qty
  kind: molar enthalpy
  desc: reference enthalpy

S_ref:
  optional: true
  type: qty
  kind: molar entropy
  desc: reference entropy
"""


yaml_io.register_class('ThermochemIncomplete',
                       yaml_io.parse(ThermochemIncomplete._yaml_schema),
                       ThermochemIncomplete)
