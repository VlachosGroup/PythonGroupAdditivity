from .utils import is_zero
from .qty import eval_qty


__all__ = ['with_units', 'in_units', 'has_units', 'to_SI_from', 'from_SI_to']


def with_units(number, units):
    """
    Construct a physical quantity from the given `number` and `units`.

    If `number` is `None`, then return `None`, regardless of `units`.

    .. note::
        If `number` is not `None`, then the following two expressions are
        equivalent:

            1. ``chemtk.units.with_units(number, units)``
            2. ``number*chemtk.units.eval_qty(units)``

    Parameters
    ----------
    number : int, float, or NoneType
        Specify the numerical part of the new quantity
    units : str or :class:`Quantity`
        Describe the units of the new quantity

    Returns
    -------
    qty : :class:`Quantity` or NoneType
        The new quantity if `number` is not None; otherwise `None`
    """
    if number is None:
        return None
    if is_zero(number):
        return number
    return number*eval_qty(units)


def in_units(qty, units):
    """
    Convert quantity to specified `units` and return numerical value.

    If `qty` is `None`, then return `None`, regardless of `units`.

    .. note::
        If `qty` is not `None`, then the following two expressions are
        equivalent:

            1. ``chemtk.units.in_units(qty, units)``
            2. ``qty.in_units(units)``

    Parameters
    ----------
    qty : :class:`Quantity` or NoneType
        Specify the quantity to convert
    units : str or :class:`Quantity`
        Describe the units to convert `qty` to

    Returns
    -------
    number : int, float, or NoneType
        The numerical value of the converted quantity, if `qty` is not None;
        otherwise `None`
    """
    if qty is None:
        return None
    return qty.in_units(units)


def has_units(qty, units):
    """
    Test quantity (`qty`) for compatibility with given `units`.

    .. note::
        The following two expressions are equivalent:

            1. ``chemtk.units.has_units(qty, units)``
            2. ``qty.has_units(units)``

    Parameters
    ----------
    qty : :class:`Quantity` or NoneType
        Specify the quantity to test
    units : str or :class:`Quantity`
        Describe the units `qty` must be convertible to

    Returns
    -------
    is_compatible : bool
        True if and only if the units of `qty` can be converted to `units`.
    """
    return qty.has_units(units)


def to_SI_from(value, units):
    """
    Convert `value` having specified `units` to the equivalent SI units.

    Parameters
    ----------
    value : int or float
        Specify value of quantity to convert to SI units
    units : str or :class:`Quantity`
        Describe the units associated with `value`

    Returns
    -------
    new_value : int or float
        Numerical value of SI units quantity
    """
    return value*eval_qty(units).value


def from_SI_to(value, units):
    """
    Convert `value` in SI units to specified compatible `units`.

    Parameters
    ----------
    value : int or float
        Specify value of SI units quantity to convert
    units : str or :class:`Quantity`
        Describe the units to convert SI quantity to

    Returns
    -------
    new_value : int or float
        Numerical value of converted quantity.
    """
    return value/eval_qty(units).value
