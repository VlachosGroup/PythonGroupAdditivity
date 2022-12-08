from collections.abc import Mapping, Sequence
from warnings import warn

from .common import string, SchemaDefinitionError, InputDataError,\
    InputDataWarning
from .lib_interface import YAMLTaggedValue


__all__ = ['SchemaRepository', 'SchemaDefinitionError']


class AnonymousClass(Mapping):
    def __contains__(self, key):
        return key in self.__dict__

    def __len__(self):
        return len(self.__dict__)

    def __iter__(self):
        return iter(self.__dict__)

    def __getitem__(self, item):
        return getattr(self, item)


def _get_type_and_attrs(schema_def):
    if schema_def is None:
        return ('any', {})
    elif isinstance(schema_def, string):
        return (schema_def, {})
    elif isinstance(schema_def, Mapping):
        attrs = schema_def.copy()
        if 'type' not in attrs:
            raise SchemaDefinitionError("missing 'type' attribute in "
                                        "schema definition %r" % schema_def)
        del attrs['type']
        return (schema_def['type'], attrs)
    else:
        raise SchemaDefinitionError("invalid schema definition %r"
                                    % schema_def)


class ListLoader(object):
    def __init__(self, repo, item_type, object_class=None, tag=None):
        self.object_class = object_class
        self.tag = tag
        self.item_loader = repo.make_loader(item_type)

    def __call__(self, name, data, context):
        # params = {}
        if isinstance(data, YAMLTaggedValue):
            # Check proper subclass ...
            tag = data.tag
            data = data.value
            if self.tag is not None and tag != self.tag:
                repo = context['repo']
                if repo.is_tag_sub_type(tag, self.tag):
                    return repo.load_tagged(name, data, context, tag)
                else:
                    raise InputDataError(
                        "Expected %r object (or sub type), "
                        "got %r object" % (self.tag, tag))
        elif not isinstance(data, Sequence):
            raise InputDataError('Input data is not sequence type: %r' % data)

        # Load parameters with multiple alternatives.
        result = [self.item_loader(name, item, context) for item in data]

        # Return constructed value.
        if self.object_class is not None:
            # Class of object is defined.
            if hasattr(self.object_class, 'yaml_construct'):
                # User provided a constructor.
                try:
                    return self.object_class.yaml_construct(result, context)
                except Exception as exc:
                    from sys import exc_traceback
                    from traceback import format_tb
                    raise InputDataError(
                        '%s.yaml_construct() reported error '
                        '"%s" with traceback:\n%s'
                        % (self.object_class.__name__, exc,
                            ''.join(format_tb(exc_traceback))))
            else:
                return result
        else:
            return result


class ObjectLoader(object):
    def __init__(
            self, repo, members, object_class=None, tag=None,
            has_open_namespace=False):
        if not isinstance(members, Mapping):
            raise InputDataError('Provided members %r is not a mapping type'
                                 % members)

        self.object_class = object_class
        self.tag = tag
        self.has_open_namespace = has_open_namespace
        members = members.copy()

        # Convert parameter names to unicode strings.
        renames = []
        for name in members:
            if not isinstance(name, str):
                old_name = name
                renames.append((old_name, string(name)))
        for (old_name, new_name) in renames:
            entry = members[old_name]
            del members[old_name]
            members[new_name] = entry

        # Validate schema attributes and collect conflicts and alternatives.
        # conflicts = {}  -- TODO
        self.loaders = {}
        self.descs = {}
        self.requireds = set()
        self.optionals = set()
        self.defaults = {}
        self.alt_sets = set()
        alt_sets_index = {}

        for name in members:
            # Determine attributes of parameter.
            schema_def = members[name]
            typ, attrs = _get_type_and_attrs(schema_def)

            # Extract description of parameter (if one is available).
            self.descs[name] = attrs.get('desc')

            # Look for special attributes: optional, default, or alts:
            if attrs.get('optional'):
                if 'default' in attrs or attrs.get('alts'):
                    raise SchemaDefinitionError("use of 'default' or 'alts'",
                                                "attributes conflict with use",
                                                "of 'optional' in %r"
                                                % schema_def)
                self.optionals.add(name)
            elif 'default' in attrs:
                if attrs.get('optional') or attrs.get('alts'):
                    raise SchemaDefinitionError("use of 'optional' or 'alts'",
                                                "attributes conflicts with",
                                                "use of 'optional' in %r"
                                                % schema_def)
                self.defaults[name] = attrs['default']
            elif 'alts' in attrs:
                if attrs.get('optional') or 'default' in attrs:
                    raise SchemaDefinitionError("use of 'optional' or",
                                                "'default' attributes",
                                                "conflicts with use of 'alts'",
                                                "in %r" % schema_def)

                # Get list of alternative parameter names:
                alts = attrs['alts']
                if isinstance(alts, str):
                    alts = [alts]
                if not isinstance(alts, Sequence):
                    raise SchemaDefinitionError("'alts' attribute must be a",
                                                "sequence instead of %r."
                                                % alts)

                # Compute the complete set of alternatives.
                class hashable_set(set):
                    def __hash__(self):
                        return id(self)

                    def __eq__(self, other):
                        return id(self) == id(other)

                alt_set = hashable_set(alts)
                self.alt_sets.add(alt_set)
                alt_set.add(name)
                for alt in alts:
                    if alt in alt_sets_index:
                        if alt_sets_index[alt] in self.alt_sets:
                            self.alt_sets.remove(alt_sets_index[alt])
                        alt_set.update(alt_sets_index[alt])  # A = A U B

                # Make sure all alt_sets_index entries refer to the newly
                # created set.
                for alt in alt_set:
                    alt_sets_index[alt] = alt_set
            else:
                self.requireds.add(name)

            # Create new schema.
            self.loaders[name] = repo.make_loader(schema_def)

    def __call__(self, name, data, context):
        params = {}
        if isinstance(data, YAMLTaggedValue):
            # Check proper subclass ...
            tag = data.tag
            data = data.value
            if self.tag is not None and tag != self.tag:
                repo = context['repo']
                if repo.is_tag_sub_type(tag, self.tag):
                    return repo.load_tagged(name, data, context, tag)
                else:
                    raise InputDataError(
                        "Expected %r object (or sub type), "
                        "got %r object" % (self.tag, tag))
        # elif data is None:
        #    return None
        elif not isinstance(data, Mapping):
            raise InputDataError('Input data for %r is not mapping type: %r'
                                 % (name, data))

        # Check for parameters not defined in schema.
        for name in data:
            if name not in self.loaders and not self.has_open_namespace:
                warn('Invalid parameter %r specified in YAML input'
                     % name, InputDataWarning, 3)
                # raise InputDataError('Invalid parameter %r specified in %r'
                #    %(name, data))

        # Load required parameters.
        for name in self.requireds:
            if name not in data:
                raise InputDataError('Missing required parameter %r in %r'
                                     % (name, data))
            params[name] = self.loaders[name](name, data[name], context)

        # Load optional parameters.
        for name in self.optionals:
            if name in data:
                params[name] = self.loaders[name](name, data[name], context)

        # Load parameters with defaults.
        for name in self.defaults:
            if name in data:
                params[name] = self.loaders[name](name, data[name], context)
            else:
                params[name] = self.loaders[name](
                    name, self.defaults[name], context)

        # Load parameters with multiple alternatives.
        for alt_set in self.alt_sets:
            alt_names = [name for name in alt_set if name in data]
            if len(alt_names) != 1:
                raise InputDataError("Number of alternatives %r specified in",
                                     "%r does not equal 1" % (alt_set, data))
            name = alt_names[0]
            params[name] = self.loaders[name](name, data[name], context)

        # Return constructed value.
        if self.object_class is not None:
            # Class of object is defined.
            if hasattr(self.object_class, 'yaml_construct'):
                # User provided a constructor.
                params['name'] = name
                try:
                    return self.object_class.yaml_construct(params, context)
                except Exception as exc:
                    from sys import exc_traceback
                    from traceback import format_tb
                    raise InputDataError(
                        '%s.yaml_construct() reported error',
                        '"%s" with traceback:\n%s'
                        % (self.object_class.__name__, exc,
                            ''.join(format_tb(exc_traceback))))
            else:
                # Default constructor.
                obj = type(self.object_class.__name__, (object,), {})()
                for name in params:
                    setattr(obj, name, params[name])
                obj.__class__ = self.object_class
                return obj
        else:
            # Class of object is undefined.
            obj = AnonymousClass()
            for name in params:
                setattr(obj, name, params[name])
            return obj


class SchemaRepository(object):
    def __init__(self, builtins=True):
        if builtins:
            # Setup builtin schemata.
            from .builtins import loaders
            self._loaders = loaders
        else:
            # No builtin schemata.
            self._loaders = {}
        self._loaders['object'] = ObjectLoader

    def make_loader(self, schema_def):
        """Return a function that accepts an abstract tree conforming to the
        provided schema definition and constructs a Python object.
        """
        type_name, attrs = _get_type_and_attrs(schema_def)
        loader_class = self._loaders.get(type_name)
        if loader_class is None:
            raise SchemaDefinitionError('undefined data type %r' % type_name)

        def del_key(attrs, key):
            if key in attrs:
                del attrs[key]
        for key in ['optional', 'default', 'alts', 'type', 'desc']:
            del_key(attrs, key)
        try:
            return loader_class(self, **attrs)
        except TypeError:
            import sys
            from traceback import extract_tb
            typ, value, tb = sys.exc_info()
            mod_name, lineno, func_name, code = extract_tb(tb)[-1]
            if(mod_name == __file__
                    and func_name == sys._getframe().f_code.co_name):
                raise SchemaDefinitionError(
                    'Invalid attributes to type %r: %r' % (type_name, attrs))
            else:
                raise

    class ObjectLoaderFactory(object):
        def __init__(self, members, object_class=None, tag=None):
            self.members = members
            self.object_class = object_class
            self.tag = tag

        def __call__(self, repo, **attrs):
            return ObjectLoader(repo, self.members, self.object_class,
                                self.tag, **attrs)

    class ListLoaderFactory(object):
        def __init__(self, item_type, object_class=None, tag=None):
            self.item_type = item_type
            self.object_class = object_class
            self.tag = tag

        def __call__(self, repo, **attrs):
            return ListLoader(repo, self.item_type, self.object_class,
                              self.tag, **attrs)

    def load_tagged(self, name, data, context, tag):
        loader = self._loaders.get(tag[1:])
        if loader is None:
            raise InputDataError('Unknown tag: %r' % tag)
        return loader(self)(name, data, context)

    def is_tag_sub_type(self, tag, base_tag):
        loader = self._loaders.get(tag[1:])
        base_loader = self._loaders.get(base_tag[1:])
        if loader is None:
            raise InputDataError('Invalid tag: %r' % tag)
        if base_loader is None:
            raise InputDataError('Invalid tag: %r' % base_tag)
        object_class = loader.object_class
        base_object_class = base_loader.object_class
        if object_class is not None and base_object_class is not None:
            return issubclass(object_class, base_object_class)

    def make_object_loader(self, members, object_class=None,
                           open_namespace=False):
        """Return a function that accepts an abstract tree conforming to a
        schema for an object whose `members` are defined by a yaml mapping.
        """
        return self.ObjectLoaderFactory(members, object_class)(
            self, has_open_namespace=open_namespace)

    def register_class(self, name, schema_def, class_):
        """Register a Python class to be used to construct objects having
        specified YAML tag (`name`) with provided schema definition.
        """
        self.register_type(name,
                           self.ObjectLoaderFactory(schema_def, class_,
                                                    tag='!' + name))

    def register_list_type(self, name, item_type, class_):
        """Register a Python class to be used to construct objects having
        specified YAML tag (`name`) in which the data is expected to be a list
        containing items of specified `item_type`.
        """
        self.register_type(name,
                           self.ListLoaderFactory(item_type,
                                                  class_, tag='!' + name))

    def register_type(self, name, loader):
        """Register a function to be called (`loader`) to construct new
        objects from YAML data tagged with `name`.
        """
        self._loaders[name] = loader

    def unregister_type(self, name):
        del self._schemata[name]
