from .lib_interface import YAMLTaggedValue, parse
from .schema import SchemaRepository
# from .common import *


__all__ = ['Error', 'load', 'make_loader', 'make_object_loader',
           'register_class', 'register_list_type', 'register_type', 'parse']


class Error(Exception):
    """Exception raised upon YAML input/output related error."""


def load(data, context=None, name=None, tag=None, loader=None):
    """Construct object from abstract tree obtained using `parse()`."""
    if context is None:
        context = {}
    context['repo'] = _repository

    if tag is not None:
        if loader is not None:
            raise Error("Specification of keyword argument 'tag' conflicts "
                        "with specification of 'loader' in yaml_io.load()")
        return _repository.load_tagged(name, data, context, tag)
    elif loader is not None:
        if isinstance(data, YAMLTaggedValue):
            raise Error("Object %r with tag %r cannot be loaded with "
                        "loader %r in yaml_io.load().", loader)
        return loader(name, data, context)
    elif not isinstance(data, YAMLTaggedValue):
        raise Error("Cannot load untagged object %r; no tag or loader "
                    "specified to yaml_io.load()." % data)
    else:
        return _repository.load_tagged(name, data.value, context, data.tag)


_repository = SchemaRepository()

make_loader = _repository.make_loader
make_object_loader = _repository.make_object_loader
register_class = _repository.register_class
register_list_type = _repository.register_list_type
register_type = _repository.register_type
