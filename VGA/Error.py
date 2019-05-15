# 2016. Vlachos Group Geun Ho Gu. University of Delaware.
"""
==========
Exceptions
==========

All the exceptions are listed here.

"""
import re
drawmolonerror = 0
if drawmolonerror:
    from .DrawMol import moltosvg

__all__ = []


class SemiEmpiricalMethodError(Exception):
    """ errors related to semi-empirical methods """


class GroupSyntaxError(Exception):
    """Exception raised when the group parser encounters invalid input."""


class GroupSchemeError(Exception):
    """Exception raised when a group scheme-related error occurs."""


class GroupMissingDataError(Exception):
    """
    Exception raised when a group library lacks the requested data for
    some group.

    Attributes
    ----------
    groups : list of :class:`chemtk.groupadd.Group`
        Groups for which data are missing.
    property_set_name : str
        Name of property set for which data are missing.
    """
    def __init__(self, groups, property_set_name):
        self.groups = groups
        self.property_set_name = property_set_name

    def __str__(self):
        if len(self.groups) == 1:
            return ('Library does not define property set %r for group %s'
                    % (str(self.property_set_name), str(self.groups[0])))
        else:
            return ('Library does not define property set %r for groups { %s }'
                    % (str(self.property_set_name),
                       ', '.join(repr(str(group)) for group in self.groups)))


class PatternMatchError(Exception):
    """
    Exception raised when no known pattern matches part of a chemical
    structure.

    The atom and (possibly) the bond at which the matching failed are stored
    as attributes :attr:`atom` and :attr:`bond`.

    For center pattern matches, :attr:`atom` specifies the atom at which the
    matching failed and :attr:`bond` is None.

    For peripheral pattern matches, :attr:`atom` is the atom adjacent to
    the bond (:attr:`bond`), that is being considered in the peripheral match.

    Attributes
    ----------
    atom : :class:`chemtk.structure.Atom`
        The afflicted atom
    bond : :class:`chemtk.structure.Bond`
        The afflicted bond (possibly None)
    """
    def __init__(self, mesg, atom):
        self.mesg = mesg
        self.atom = atom
        self.mol = atom.GetOwningMol()
        self.visualize()

    def visualize(self, *args, **kwargs):
        if drawmolonerror:
            moltosvg(self.mol, highlight=[self.atom.GetIdx()], kekulize=False)

    def __str__(self):
        return '%s at atom number %s' % (self.mesg, self.atom.GetIdx())


__all__ += ['SemiEmpiricalMethodError', 'GroupSyntaxError', 'GroupSchemeError',
            'GroupMissingDataError', 'PatternMatchError']


# RING-related errors.
class RINGError(Exception):
    """Base exception for Parser errors.

    All other Parser-related exceptions inherit from this.
    """


class RINGSyntaxError(RINGError):
    """
    Exception raised due to invalid (but syntactically correct) input.
    """
    _parse = re.compile('[\n]')

    def __init__(self, tok, lineno, colno, stream):
        self.toks = set([tok])
        self.lineno = lineno
        self.colno = colno
        self.stream = stream

    def update(self, other):
        if (self.lineno < other.lineno or
                (self.lineno == other.lineno and self.colno < other.colno)):
            self.lineno = other.lineno
            self.colno = other.colno
            self.toks = other.toks.copy()
        elif (self.lineno == other.lineno and self.colno == other.colno):
            self.toks |= other.toks

    def __str__(self):
        s = 'Expected ' + ' | '.join(
            sorted(str(tok) for tok in self.toks if tok))\
            + ' at line %d column %d:\n' % (self.lineno, self.colno)
        s += self._parse.split(self.stream)[self.lineno-1] + '\n'
        s += ' '*(self.colno-1) + '^\n'
        return s

    def __repr__(self):
        return '%s(%r, %r, %r, %r)' % (
            type(self).__name__, self.toks, self.lineno, self.colno,
            self.smiles)


class RINGReaderError(RINGError):
    """
    Exception raised when input does not conform to RING syntax.
    """
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message

    def __repr__(self):
        return '%s(%r, %r)' % (type(self).__name__, self.message)


class MolQueryError(Exception):
    """
    Exception raised when input does not conform to RING syntax.
    """
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message

    def __repr__(self):
        return '%s(%r, %r)' % (type(self).__name__, self.message)


class ReactionQueryError(Exception):
    """
    Exception raised when input does not conform to RING syntax.
    """
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message

    def __repr__(self):
        return '%s(%r, %r)' % (type(self).__name__, self.message)


__all__ += ['RINGError', 'RINGSyntaxError', 'RINGReaderError',
            'MolQueryError', 'ReactionQueryError']


# Units errors.
class UnitsParseError(Exception):
    """Error from providing bad input to the units parser."""


class UnitsError(Exception):
    """Error from operation involving incompatible physical units."""


__all__ += ['UnitsParseError', 'UnitsError']


# Generic errors.
class OutsideCorrelationError(Exception):
    """Error from attempt to evaluate correlation outside its valid range."""


class ReadOnlyDataError(Exception):
    """Error raised by attempt to modify read-only data."""


class IncompleteDataError(Exception):
    """Error raised when a computation requires more data than is available."""


class IncompleteDataWarning(Warning):
    """
    Warning issued when a computation proceeds using less data than is optimal.

    This is raised when the absence of certain non-required data may lead to
    pontetially severe assumptions in later computations.
    """


__all__ += ['OutsideCorrelationError', 'ReadOnlyDataError',
            'IncompleteDataError', 'IncompleteDataWarning']
