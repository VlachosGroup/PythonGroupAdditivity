from collections.abc import Mapping, Sequence

from .. import Units

from .common import string, InputDataError
# from .schema import *
# from .lib_interface import *


# All classes named _*_schema define builtin schemas for SchemaRepository
# objects.


class any_loader(object):
    def __init__(self, repo):
        pass

    def __call__(self, name, value, context):
        return value


class qty_loader(object):
    def __init__(self, repo, kind=None):
        self.kind = kind

    def __call__(self, name, value, context):
        if value is None:
            return None
        if 'units' in context:
            kind_units = context['units'].get(self.kind)
        else:
            kind_units = None

        qty = Units.eval_qty(value)
        if not isinstance(qty, Units.Quantity):
            if kind_units is not None:
                return Units.with_units(qty, kind_units)
            else:
                raise InputDataError(
                    'Cannot determine units of quantity: %r' % qty)
        return qty


class bool_loader(object):
    def __init__(self, repo):
        pass

    def __call__(self, name, value, context):
        if value is None:
            return None
        try:
            return bool(value)
        except ValueError:
            raise InputDataError('Invalid int: %r' % value)


class float_loader(object):
    def __init__(self, repo):
        pass

    def __call__(self, name, value, context):
        if value is None:
            return None
        try:
            return float(value)
        except ValueError:
            raise InputDataError('Invalid float: %r' % value)


class int_loader(object):
    def __init__(self, repo):
        pass

    def __call__(self, name, value, context):
        if value is None:
            return None
        try:
            return int(value)
        except ValueError:
            raise InputDataError('Invalid int: %r' % value)


class number_loader(object):
    def __init__(self, repo):
        pass

    def __call__(self, name, value, context):
        if value is None:
            return None
        if isinstance(value, int) or isinstance(value, float):
            return value
        try:
            return int(value)
        except ValueError:
            try:
                return float(value)
            except ValueError:
                raise InputDataError('Invalid int: %r' % value)


class string_loader(object):
    def __init__(self, repo):
        pass

    def __call__(self, name, value, context):
        if value is None:
            return None
        try:
            return string(value)
        except ValueError:
            raise InputDataError('Invalid string: %r' % value)


class enum_loader(object):
    # TODO: more elegant -- support type coercion of values to
    # allowable types.
    def __init__(self, repo, values):
        # if values is None:
        #    raise InputDataError("In schema: enum type description is missing"
        #        "required attribute 'values'")
        # elif not isinstance(values, Sequence):
        if not isinstance(values, Sequence):
            raise InputDataError("In schema: enum type description attribute",
                                 "'values' must be a sequence type.  Instead",
                                 "got %r" % values)
        self.values = set(values)

    def __call__(self, name, value, context):
        if value not in self.values:
            raise InputDataError('Enumerated value %r is not one of %r'
                                 % (value, self.values))
        return value


class tuple_loader(object):
    def __init__(self, repo, item_types):
        if not isinstance(item_types, Sequence):
            raise InputDataError(
                'In tuple schema: item_types must be a sequence')
        self.loaders = [repo.make_loader(item_type)
                        for item_type in item_types]

    def __call__(self, name, values, context):
        if (not isinstance(values, Sequence)
                or len(values) != len(self.loaders)):
            raise InputDataError(
                'Expected sequence of %d values, got %r instead'
                % (len(self.loaders), values))
        items = []
        for (i, loader) in enumerate(self.loaders):
            items.append(loader(None, values[i], context))
        return items


class list_loader(object):
    def __init__(self, repo, item_type=None):
        self.loader = repo.make_loader(item_type)

    def __call__(self, name, input_values, context):
        if not isinstance(input_values, Sequence):
            raise InputDataError(
                'Expected sequence, got %r instead' % input_values)
        return [self.loader(string(i), item, context)
                for (i, item) in enumerate(input_values)]


class array_loader(object):
    def __init__(self, repo, shape=None, dtype='d'):
        if not isinstance(shape, Sequence):
            shape = (shape,)
        self.shape = shape
        self.dtype = dtype
        # self.loader = repo.make_loader(item_type)

    def __call__(self, name, input_values, context):
        import numpy
        if not isinstance(input_values, Sequence):
            raise InputDataError(
                'Expected sequence, got %r instead' % input_values)
        ar = numpy.array(input_values, dtype=self.dtype)
        if self.shape is not None and ar.shape != self.shape:
            raise InputDataError("Shape of array %r does not conform to "
                                 "required shape %r" % (ar.shape, self.shape))
        return ar


class mapping_loader(object):
    def __init__(self, repo, name_type='string', value_type=None):
        self.name_loader = repo.make_loader(name_type)
        self.value_loader = repo.make_loader(value_type)

    def __call__(self, name, dct, context):
        if not isinstance(dct, Mapping):
            raise InputDataError(
                'Expected mapping, got %r instead' % (dct))
        return dict(
            (self.name_loader(name, name, context),
                self.value_loader(name, dct[name], context))
            for name in dct)


class named_object_loader(object):
    def __init__(self, repo, members):
        self.loader = repo.make_object_loader(members)

    def __call__(self, name, input_value, context):
        k, v = list(input_value.items())[0]
        return self.loader(k, v, context)


class named_value_loader(object):
    def __init__(self, repo, value_type, name_type='string'):
        self.name_loader = repo.make_loader(name_type)
        self.value_loader = repo.make_loader(value_type)

    def __call__(self, name, input_value, context):
        name, value = list(input_value.items())[0]
        return (self.name_loader(None, name, context),
                self.value_loader(None, value, context))


_globals = globals().copy()
loaders = dict((name[:-7], _globals[name]) for name in _globals
               if name.endswith('_loader') and
               isinstance(_globals[name], type))
