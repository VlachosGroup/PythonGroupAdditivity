from .. Error import UnitsParseError


__all__ = ['units_db']


class UnitsDB(object):
    prefixes = {
        'Y': 1e24,
        'Z': 1e21,
        'E': 1e18,
        'P': 1e15,
        'T': 1e12,
        'G': 1e9,
        'M': 1e6,
        'k': 1e3,
        'h': 1e2,
        'da': 1e1,
        'd': 1e-1,
        'c': 1e-2,
        'm': 1e-3,
        'u': 1e-6,
        'n': 1e-9,
        'p': 1e-12,
        'f': 1e-15,
        'a': 1e-18,
        'z': 1e-21,
        'y': 1e-24,
    }

    def __init__(self):
        self.db = {}

    def lookup(self, name):
        # Try without prefix first
        if name in self.db:
            return self.db[name]
        # Try with single letter prefix:
        if name[1:] in self.db and name[:1] in self.prefixes:
            return self.prefixes[name[:1]]*self.db[name[1:]]
        # Try with double letter prefix (just 'da'):
        if name[2:] in self.db and name[:2] in self.prefixes:
            return self.prefixes[name[1:]]*self.db[name[2:]]
        # Raise error.
        raise UnitsParseError('Unknown units: %r' % name)

    def add(self, name, val):
        self.db[name] = val

# Define global units database.


units_db = UnitsDB()
