import numpy as np

from .. error import UnitsError

from .utils import is_zero
from .parser import *


__all__ = [
    'eval_quantity', 'eval_qty', 'Quantity', 'ArrayQuantity', 
    'FundamentalUnits']


def eval_quantity(expr):
    """
    Evaluate string expression of physical quantity or units.

    If `expr` is not a string type, then `expr` is returned unchanged.

    Parameters
    ----------
    expr : str
        Describe physical quantity

    Returns
    -------
    qty : :class:`Quantity`
        Representation of physical quantity
    """
    if isinstance(expr, basestring):
        return eval_expr(expr)
    else:
        return expr
eval_qty = eval_quantity  # Abbreviation for above (gets used a lot).


class GenericQuantity(object):
    @staticmethod
    def _unpack_qty(other):
        if isinstance(other, Quantity):
            value = other.value
            units = other.units
        elif isinstance(other, ArrayQuantity):
            value = other.view(np.ndarray)
            units = other._units
        elif isinstance(other, (list, tuple)):
            value = np.array(other)
            units = FundamentalUnits.null()
        else:
            value = other
            units = FundamentalUnits.null()
        return (value, units)

    @classmethod
    def _build(cls, value, units):
        if not units:
            return value

        if isinstance(value, (float, int, long)):
            return Quantity(value, units)
        elif isinstance(value, np.ndarray):
            new_value = value.view(ArrayQuantity)
            new_value._units = units
            return new_value
        else:
            raise AssertionError(
                'Unknown value type to GenericQuantity._build(): %r'%value)

    """
    def __new__(cls, value, units):
        q = cls._numeric_class.__new__(cls, value)
        q._units = units
        return q
    """

    def in_units(self, units):
        """
        Convert to specified `units` and return numerical value.

        Parameters
        ----------
        units : str or instance of :class:`Quantity`
            Describe the units to convert to before formatting

        Returns
        -------
        formatted : int or float
            The numerical value of the converted quantity

        Raises
        ------
        UnitsError
            If `units` are incompatible with the units of this quantity.
        """
        value = self/eval_qty(units)
        if isinstance(value, GenericQuantity):
            raise UnitsError("Units of '%s' are incompatible with '%s'"
                %(self, units))
        return value

    def has_units(self, units):
        """
        Test for compatibility with given `units`.

        Parameters
        ----------
        units : str or instance of :class:`Quantity`
            Describe the units this quantity must be convertible to

        Returns
        -------
        is_compatible : bool
            `True` if and only if conversion to `units` is possible.
        """
        (self_value, self_units) = self._unpack_qty(self)
        if isinstance(units, FundamentalUnits):
            return self_units == units
        else:
            return self_units == self._unpack_qty(eval_qty(units))[1]

    def __eq__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        if not is_zero(other) and (not other_units or
                not self.has_units(other_units)):
            return False
        return self_value == other_value

    def __ne__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        if not is_zero(other) and (not other_units or
                not self.has_units(other_units)):
            return True
        return self_value != other_value

    def __lt__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        if (not is_zero(other_value) and
                (not other_units or not self.has_units(other_units))):
            raise UnitsError(
                'Incompatible units %s vs %s in comparison'
                %(self_units, other_units))
        return self_value < other_value

    def __lt__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        if (not is_zero(other_value) and
                (not other_units or not self.has_units(other_units))):
            raise UnitsError(
                'Incompatible units %s vs %s in comparison'
                %(self_units, other_units))
        return self_value > other_value

    def __ge__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        if (not is_zero(other_value) and
                (not other_units or not self.has_units(other_units))):
            raise UnitsError(
                'Incompatible units %s vs %s in comparison'
                %(self_units, other_units))
        return self_value >= other_value

    def __le__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        if (not is_zero(other_value) and
                (not other_units or not self.has_units(other_units))):
            raise UnitsError(
                'Incompatible units %s vs %s in comparison'
                %(self_units, other_units))
        return self_value <= other_value

    def __add__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        if (not is_zero(other_value) and
                (not other_units or not self.has_units(other_units))):
            raise UnitsError(
                'Incompatible units %s vs %s in addition'
                %(self_units, other_units))
        return self._build(self_value + other_value, self_units)

    def __radd__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        if (not is_zero(other_value) and
                (not other_units or not self.has_units(other_units))):
            raise UnitsError(
                'Incompatible units %s vs %s in addition'
                %(self_units, other_units))
        return self._build(other_value + self_value, self_units)

    def __sub__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        if (not is_zero(other_value) and
                (not other_units or not self.has_units(other_units))):
            raise UnitsError(
                'Incompatible units %s vs %s in subtraction'
                %(self_units, other_units))
        return self._build(self_value - other_value, self_units)

    def __rsub__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        if (not is_zero(other_value) and
                (not other_units or not self.has_units(other_units))):
            raise UnitsError(
                'Incompatible units %s vs %s in subtraction'
                %(self_units, other_units))
        return self._build(other_value - self_value, self_units)

    def __mul__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        return self._build(self_value*other_value, self_units*other_units)

    def __rmul__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        return self._build(other_value*self_value, other_units*self_units)

    def __div__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        return self._build(self_value/other_value, self_units/other_units)
    __truediv__ = __div__

    def __rdiv__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        return self._build(other_value/self_value, other_units/self_units)
    __rtruediv__ = __rdiv__

    def __pow__(self, other):
        (self_value, self_units) = self._unpack_qty(self)
        (other_value, other_units) = self._unpack_qty(other)
        if other_units:
            raise TypeError('Invalid operation: exponentiation by quantity.')
        return self._build(
            self_value**other_value, self_units**other_value)

    def __rpow__(self, other):
        raise TypeError('Invalid operation: exponentiation by quantity.')

    def __neg__(self):
        (self_value, self_units) = self._unpack_qty(self)
        return self._build(-self_value, self_units)

    def __abs__(self):
        (self_value, self_units) = self._unpack_qty(self)
        return self._build(abs(self_value), self_units)

    def __str__(self):
        (value, units) = self._unpack_qty(self)
        return str(value) + ' ' + str(units)

    def __repr__(self):
        (value, units) = self._unpack_qty(self)
        return '%s(%r, %r)'%(type(self).__name__, value, units)


class Quantity(GenericQuantity):
    """
    Represent a physical quantity.

    .. note::
        This class should not be instantiated directly by users.  Rather, users
        should construct quantities and units using :func:`eval_qty()` or
        another helper function.

    Physical quantities usually possess units of some kind, and these are
    stored internally as an instance of :class:`FundamentalUnits`.
    Effectively, physical quantities are always stored in SI units to simplify
    arithmetic.  Various method are available to format or numerically convert
    them to any available compatible units.

    To facilitate meaningful computation with unitful quantities, this
    representation overloads several arithmetic and comparison operators:

        * Comparison operators (``==``, ``!=``, ``<``, ``<=``, ``>=``, and
          ``>``) are supported between quantities with compatible units.
        * Addition (``+``) and subtraction (``-``) are supported between
          quantities with compatible units.
        * Multiplication (``*``) and division (``/``) are supported between
          two quantities and between a scalar and a quantity.
        * Exponentation (``**``) is supported for a quantity as the base and a
          scalar as the exponent.
        * Conversion to type `bool`: The converted value is `True`, if and only
          if the quantity is non-zero.
        * Conversion to type `int` or type `float`: The converted value is
          the value of the quantity in equilvalent SI units.

    If two quantities involved have incompatible units and such compatiblity
    is required, the operation raises :exc:`chemtk.error.UnitsError`.

    If an arithmetic operation leads to a unitless result, then a number will
    be returned instead of a :class:`Quantity`.
    """
    __array_priority__ = 1.0

    def __init__(self, value, units):
        self.value = value
        if isinstance(units, basestring):
            units = eval_qty(units).units
        self.units = units

    def __hash__(self):
        return hash(self.value)

    def __int__(self):
        return int(self.value)

    def __float__(self):
        return float(self.value)

    def fmt_in_units(self, units):
        """
        Format in specified `units`.

        Convert this quantity to `units` and return a string representation of
        the converted quantity.

        Parameters
        ----------
        units : str or instance of :class:`Quantity`
            Describe the units to convert to before formatting

        Returns
        -------
        formatted : str
            String representation of quantity after converting to `units`

        Raises
        ------
        UnitsError
            If `units` are incompatible with the units of this quantity.
        """
        return '%g %s'%(self.in_units(units), units)


class ArrayQuantity(GenericQuantity, np.ndarray):
    def __new__(cls, data, *args, **kwargs):
        units = kwargs.pop('units', None)
        if units:
            units = eval_qty(units).units
        if(hasattr(data, '__iter__')
                and all((isinstance(datum, GenericQuantity) or datum == 0)
                    for datum in data)):
            # Bundled together quantities.
            values = []
            use_units = None
            for datum in data:
                (value, datum_units) = cls._unpack_qty(datum)
                if not use_units and datum_units:
                    use_units = datum_units
                if(use_units and datum_units and use_units != datum_units
                        and value != 0 and not units):
                    raise UnitsError(
                        'Inconsistent units in contents provided to '
                        'ArrayQuantity initializer')
                values.append(value)

            data = values
            if units:
                if use_units != units:
                    raise UnitsError(
                        'Units in contents do not match specified units in '
                        'ArrayQuantity initializer')
            else:
                units = use_units

        arrqty = np.array(data, *args, **kwargs)
        if not issubclass(arrqty.dtype.type, np.number):
            raise TypeError('Units may only be assigned to numerical values.')

        arrqty = arrqty.view(cls)
        if units:
            arrqty._units = units
        elif hasattr(arrqty, '_units'):
            arrqty._units = arr._units
        else:
            raise TypeError(
                'Units must be specified when initializing ArrayQuantity.')
        return arrqty

    def __array_finalize__(self, arrqty):
        if hasattr(arrqty, '_units'):
            self._units = arrqty._units

    def __array_wrap__(self, outarr, context=None):
        if not issubclass(outarr.dtype.type, np.number):
            return np.asarray(outarr)
        return outarr

    def __getitem__(self, idx):
        result = np.ndarray.__getitem__(self, idx)
        if not isinstance(result, np.ndarray):
            return Quantity(result, self._units)
        return result

    def fmt_in_units(self, units):
        """
        Format in specified `units`.

        Convert this quantity to `units` and return a string representation of
        the converted quantity.

        Parameters
        ----------
        units : str or instance of :class:`Quantity`
            Describe the units to convert to before formatting

        Returns
        -------
        formatted : str
            String representation of quantity after converting to `units`

        Raises
        ------
        UnitsError
            If `units` are incompatible with the units of this quantity.
        """
        return '%s %s'%(self.in_units(units), units)


class FundamentalUnits(object):
    """
    Represent the fundamental units of a physical quantity.

    .. note::
        This class should not be instantiated directly by users.  Rather, users
        should construct quantities and units using :func:`eval_qty()` or
        another helper function.

    Units are represented internally as an array of numerical exponents, one
    each for the standard SI units:

        * meter (m)
        * kilogram (kg)
        * second (s)
        * ampere (A)
        * kevlin (K)
        * gram mole (mol)
        * candela (cd)

    The representation overloads multiplication (``*``), division (``/``),
    and equality (``==``) and inequality (``!=``) comparisons for use between
    instances of :class:`FundamentalUnits`.  Additionally, the exponentation
    (``**``) operator is suported between a :class:`FundamentalUnits` instance
    and a numerical exponent.  The truth value (*i.e.* ``bool(units)``) of an
    instance of :class:`FundamentalUnits` is `True` if and only if at least one
    exponent is non-zero.

    This abstraction allows for floating point numbers on the exponents,
    however, exponent values within a *threshold* (currently 10^-7) of an
    integer value are rounded to avoid comparison failures arising due to the
    inexactness of floating point arithmetic.
    """
    THRESHOLD_INTEGER = 1e-7
    _primitive_units = ['m', 'kg', 's', 'A', 'K', 'mol', 'cd']
    _primitive_units_index = dict((unit, i)
        for (i, unit) in enumerate(_primitive_units))

    @classmethod
    def new(cls, name):
        exps = np.zeros(len(cls._primitive_units))
        are_floats = np.zeros(len(cls._primitive_units), dtype='bool')
        exps[cls._primitive_units_index[name]] = 1
        return cls(exps, are_floats)

    @classmethod
    def null(cls):
        exps = np.zeros(len(cls._primitive_units))
        are_floats = np.zeros(len(cls._primitive_units), dtype='bool')
        return cls(exps, are_floats)

    @classmethod
    def _build(cls, exps):
        exps_rounded = exps.round()
        are_floats = np.abs(exps - exps_rounded) > cls.THRESHOLD_INTEGER
        return cls(np.select(
            [are_floats, np.logical_not(are_floats)],
            [exps, exps_rounded]), are_floats)

    def __init__(self, exps, are_floats):
        self.exps = exps
        self.are_floats = are_floats

    def __mul__(self, other):
        return self._build(self.exps + other.exps)

    def __div__(self, other):
        return self._build(self.exps - other.exps)

    def __pow__(self, other):
        return self._build(other*self.exps)

    def __eq__(self, other):
        return (self.exps == other.exps).all()

    def __ne__(self, other):
        return not (self == other)

    def __nonzero__(self):
        return bool(self.exps.any())

    def __str__(self):
        up = []
        dn = []
        exps = self.exps
        for i,unit in enumerate(self._primitive_units):
            if exps[i] == 1:
                up.append(unit)
            elif exps[i] > 0:
                if self.are_floats[i]:
                    up.append('%s^%s'%(unit, exps[i]))
                else:
                    up.append('%s^%s'%(unit, int(exps[i])))
            elif exps[i] == -1:
                dn.append(unit)
            elif exps[i] < 0:
                if self.are_floats[i]:
                    dn.append('%s^(%s)'%(unit, exps[i]))
                else:
                    dn.append('%s^(%s)'%(unit, int(exps[i])))
        if not up:
            s = '1'
        else:
            s = '*'.join(up)
        if not dn:
            return s
        elif len(dn) == 1:
            return s + '/' + dn[0]
        else:
            return s + '/(' + '*'.join(dn) + ')'

    def __repr__(self):
        return "<%s: '%s'>"%(type(self).__name__, str(self))
