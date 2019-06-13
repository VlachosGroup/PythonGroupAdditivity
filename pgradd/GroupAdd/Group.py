import re
from collections import defaultdict

from .. Error import GroupSyntaxError


# Note: This class is not currently documented explicitly and is only used to
# inherit Group() from.
class Descriptor(object):
    """
    Represent a chemical descriptor upon which properties may
    linearly depend on.
    """
    def __init__(self, scheme, name):
        """Initialize a chemical descriptor with given string name.

        Parameters
        ----------
        scheme : instance of :class:`GroupScheme`
            Specify group-additivity scheme for which the group belongs to.
        name : str
            Specify name of descriptor
        """
        self.scheme = scheme
        self.name = name

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        if isinstance(other, type(self)):
            return self.name == other.name
        else:
            return self.name == other

    def __neq__(self, other):
        return not (self == other)

    def __str__(self):
        return self.name

    def __repr__(self):
        return '<%s: %s>' % (type(self).__name__, str(self))


class Group(Descriptor):
    """Represent a chemical functional group.

    This object is hashable (meaning that it may be used as a dictionary key)
    and can be compared for equality against others of its kind, regardless of
    the order in which each the peripheral subgroups were specified at
    initialization.

    In addition to the default constructor (:meth:`Group.__init__()`), one can
    instantiate :class:`Group` objects using the alternative constructor
    :meth:`Group.parse()` which accepts a string name for the group that is
    parsed to determine the contained `csg` and `psgs`.
    """
    _parser_re = re.compile('[()]')

    @classmethod
    def parse(cls, scheme, text):
        """Initialize a chemical functional group given a string name.

        The syntax consists of the center-subgroup name followed by one or more
        peripheral subgroup names enclosed in parantheses, each of which is
        optionally followed by a numerical repeat count.

        Parameters
        ----------
        text : str
            Name of the group

        Examples
        --------
        Using the (builtin) ``'benson'`` scheme::

            >>> from chemtk.groupadd import Group, GroupScheme
            >>> scheme = GroupScheme.load('benson')
            >>> Group.parse(scheme, 'C(C)(H)3')
            Group('C', ['C', 'H', 'H', 'H'])
            >>> Group.parse(scheme, 'C(C[d])(C)2(H)')
            Group('C', ['C[d]', 'C', 'C', 'H'])
            >>> Group.parse(scheme, 'O(C)(H)')
            Group('O', ['C', 'H'])
        """
        def append_to_psgs(count, psg):
            for i in range(count):
                psgs.append(psg)
        parts = cls._parser_re.split(text)
        csg = parts[0]
        next_psg = None  # next nearest neighbor to be appended to nns
        psgs = []  # list of nearest neighbors
        for part in parts[1:]:
            if not part:
                continue
            if part.isdigit():
                if next_psg is None:
                    raise GroupSyntaxError('Expected base atom for Benson',
                                           'group got number instead.')
                append_to_psgs(int(part), next_psg)
                next_psg = None
            else:
                if next_psg is not None:
                    append_to_psgs(1, next_psg)
                    next_psg = None
                next_psg = part
        if next_psg is not None:
            append_to_psgs(1, next_psg)
        return cls(scheme, csg, psgs)

    def __init__(self, scheme, csg, psgs):
        """Initialize a chemical functional group.

        Parameters
        ----------
        scheme : instance of :class:`GroupScheme`
            Specify group-additivity scheme for which the group belongs to.
        csg : str
            Name of center-subgroup
        psgs : iterable of strs
            Names of peripheral-subgroups
        """
        self.csg = csg
        self.psgs = psgs
        canon_name = self._canonical_name()
        Descriptor.__init__(self, scheme, canon_name)

    def _canonical_name(self):
        """Return canonical name of group.

        Parameters
        ----------
        scheme : :class:`GroupScheme`
            Specify group-additivity scheme for which the group belongs to.

        Returns
        -------
        name : str
            Canonical name of group
        """
        canon_name = self.csg
        psg_counts = defaultdict(int)
        for psg in self.psgs:
            psg_counts[psg] += 1
        for name in sorted(psg_counts):
            if psg_counts[name] == 1:
                canon_name += '(' + name + ')'
            else:
                canon_name += '(' + name + ')' + '%d' % psg_counts[name]
        return canon_name

    def __repr__(self):
        return '%s(%r, %r)' % (type(self).__name__, self.csg, self.psgs)
