import sys


__all__ = [
    'string', 'SchemaDefinitionError', 'InputDataError', 'InputDataWarning']


if sys.version_info[0] >= 3:
    string = str
else:
    string = str
    # string = unicode


class SchemaDefinitionError(Exception):
    pass


class InputDataError(Exception):
    pass


class InputDataWarning(Warning):
    pass
