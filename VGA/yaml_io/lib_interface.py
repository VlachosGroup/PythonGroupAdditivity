import yaml

from .common import string


__all__ = ['parse', 'YAMLTaggedValue']


class YAMLTaggedValue(object):
    def __init__(self, tag, value):
        self.tag = tag
        self.value = value

    def __str__(self):
        return '<%s>: %s' % (self.tag, self.value)

    def __repr__(self):
        return 'YAMLTaggedValue(%r, %r)' % (self.tag, self.value)


def _construct_str(loader, node):
    return string(loader.construct_scalar(node))


def _construct_seq(loader, node):
    result = loader.construct_sequence(node)
    return result


def _construct_map(loader, node):
    result = dict([(tuple(key), value) if isinstance(key, list)
                  else (key, value)
                  for (key, value) in loader.construct_pairs(node, deep=True)])
    return result


def _construct_tagged(loader, tag, node):
    if isinstance(node, yaml.MappingNode):
        value = _construct_map(loader, node)
    elif isinstance(node, yaml.SequenceNode):
        value = _construct_seq(loader, node)
    elif isinstance(node, yaml.ScalarNode):
        value = loader.construct_scalar(node)
    else:
        raise TypeError('Unknown node type: %s' % type(node))
    return YAMLTaggedValue(tag, value)


try:
    # Use the class that's atop libyaml for speed...
    class Loader(yaml.CSafeLoader):
        pass
except AttributeError:
    # ...and fall back on the Python implementation otherwise:
    class Loader(yaml.SafeLoader):
        pass

# Note: other possible tags:
#   'tag:yaml.org,2002:int', 'tag:yaml.org,2002:float',
#   'tag:yaml.org,2002:bool', tag:yaml.org,2002:null',
#   and 'tag:yaml.org,2002:seq'
Loader.add_constructor('tag:yaml.org,2002:str', _construct_str)
Loader.add_constructor('tag:yaml.org,2002:map', _construct_map)

# Captures all things with an explicit tag.
Loader.add_multi_constructor('', _construct_tagged)

# Might look into these later?
# VGToolsLoader.add_multi_constructor('!', VGToolsLoader.constructor)
# VGToolsLoader.add_constructor('!!seq', VGToolsLoader.seq_constructor)


# We want YAML 1.2 support but the Python YAML library only supports YAML 1.1.
# YAML 1.1 has lots of implicit resolvers that try to automatically convert
# too many things.  We address this by replacing the implicit_resolvers with
# only those that are explicitly in the YAML 1.2 core set).
YAML12_core_implicit_tags = [
    'tag:yaml.org,2002:' + tag for tag in
    ['str', 'seq', 'map', 'null', 'bool', 'int', 'float']]
for first in Loader.yaml_implicit_resolvers:
    Loader.yaml_implicit_resolvers[first] = [
        (tag, regexp)
        for (tag, regexp) in Loader.yaml_implicit_resolvers[first]
        if tag in YAML12_core_implicit_tags]


def parse(yaml_input):
    """Parse data in YAML format to an abstract tree."""
    return yaml.load(yaml_input, Loader=Loader)
