import os
from warnings import warn
from collections import Mapping
from .. import yaml_io
import numpy as np

from .. Error import GroupMissingDataError
from . Group import Group, Descriptor
from . Scheme import GroupAdditivityScheme
from . DataDir import get_data_dir


class GroupLibrary(Mapping):
    """Represent library of contributing properties organized by group.

    The set of properties that may be represented in a :class:`GroupLibrary` is
    extensible.  Contributing properties are represented by *property sets*
    organized by group.  See the manual for a list of available property sets.

    .. note::
        Because the *property set* system is extensible, the module within
        which a particular *property set* is defined and registered must be
        imported before loading a group library that contains data for that
        type of *property set*.

    To estimate properties, call :meth:`GroupLibrary.estimate()` with the
    *property set* `name` and  set of groups contained in the chemical
    structure of interest.  The properties estimate will be returned as an
    object whose type depends on the particular *property set*.

    To determine which groups are present in a particular chemical structure,
    use :meth:`GroupLibrary.match_groups()`.

    Data in multiple group libraries can be combined so long as the groups
    they contain are defined within compatible schemes.  See
    :meth:`GroupLibrary.update()`.

    """
    _property_set_estimator_types = {}
    _property_set_group_yaml_types = {}

    @classmethod
    def register_property_set_type(cls, name, group_yaml_type, estimator_type):
        """(class method) Register new property set type.

        Parameters
        ----------
        name : str
            Name of new property set type.
        group_yaml_type : str
            Name of property set type in the YAML type namespace.
        estimator_type : class
            The provided class is instantiated when an estimate is to be made
            for a particular set of groups.  The constructor should accept the
            following parameters:

                library : :class:`GroupLibrary`
                    library which is to be used to estimate these properties.
                groups : mapping
                    Map from :class:`Group` to int or float specifying counts
                    of each group in the chemical structure.
        """
        if name in cls._property_set_group_yaml_types:
            raise KeyError('Property set %r already registered.' % name)
        cls._property_set_group_yaml_types[name] = group_yaml_type
        cls._property_set_estimator_types[name] = estimator_type

    def __init__(self, scheme, contents={}, uq_contents={}, path=None):
        """Initialize library of contributing properties organized by group.

        Parameters
        ----------
        scheme : :class:`GroupScheme`
            Specify group-additivity scheme to use.
        contents : mapping or list
            Define initial contents of the library either as mapping or list of
            (`key`, `value`) pairs.  See the last paragraph of the class
            documentation for information on the format.

        Other Parameters
        ----------------
        path : str
            File-system path library was loaded from or should be saved to by
            default.
        """
        self.scheme = scheme
        self.path = path
        if isinstance(contents, Mapping):
            contents = list(contents.items())
        self.contents = dict((group, property_sets)
                             for (group, property_sets) in contents)
        self.uq_contents = uq_contents

    def GetDescriptors(self, mol):
        """Determine groups appearing in chemical structure `chem`.

        Parameters
        ----------
        mol : :class:`rdkit.mol`
            Specify chemical structure to match groups for.
        manual_descriptors : mapping, optional
            Specify value(s)/degree(s) of influence of additional descriptors
            to include.

        Returns
        -------
        groups : mapping
            Map from :class:`Group` to int or float identifying groups and
            their number of occurence in the structure.
        """
        return self.scheme.GetDescriptors(mol)

    def Estimate(self, groups, property_set_name):
        """Estimate set of properties for chemical.

        Parameters
        ----------
        groups : mapping (dictionary)
            Map from :class:`Group` to int or float specifying counts of each
            group in the chemical structure.
        property_set_name : str
            Name of property set to estimate.

        Returns
        -------
        estimated_properties : (varies)
            The estimated properties, an object whose type depends on the
            particular property set.
        """
        if property_set_name not in self._property_set_estimator_types:
            raise KeyError('Invalid property_set name: %r' % property_set_name)
        # Verify groups present.
        missing_groups = [group for group in groups
                          if property_set_name not in self[group]]
        if missing_groups:
            raise GroupMissingDataError(missing_groups, property_set_name)
        estimator_type = self._property_set_estimator_types[property_set_name]
        return estimator_type(self, groups)

    def __contains__(self, group):
        """Test if this library contains contributing properties for `group`.

        Parameters
        ----------
        group : :class:`Group`
            Group whose membership is being tested.

        Returns
        -------
        result : bool
            True if this library has properties for `group`.
        """
        return group in self.contents

    def __iter__(self):
        """Return iterator over all groups with property data in this library.
        """
        return iter(self.contents)

    def __len__(self):
        """Return number of groups with properties in this library."""
        return len(self.contents)

    def __getitem__(self, group):
        """Return contributing properties sets for `group`.

        If no properties exist for `group`, then return ``{}`` instead of
        raising an exception.

        Parameters
        ----------
        group : :class:`Group`
            Identify group whose property sets are to be retrieved.

        Returns
        -------
        property_sets : dict
            Sets of contributing properties for `group`.
        """
        return self.contents.get(group, {})

    @classmethod
    def Load(cls, path):
        """(class method) Load group-additivity library from file-system
        `path` or builtin.

        Parameters
        ----------
        path : str
            Specify either the path to a file containing the data or a symbolic
            name of a builtin library to load (*e.g.* ``gas_benson`` to load
            gas phase Benson groups.)

        Returns
        -------
        lib : :class:`GroupLibrary`
            Group library containing the loaded data.
        """
        if os.sep not in path and '.' not in path and not os.path.exists(path):
            # [JTF] where's our data directory?
            base_path = os.path.join(get_data_dir(), path)
            # We want to load the library.yaml in that directory:
            path = os.path.join(base_path, 'library.yaml')
        else:
            # The base path is the directory containing whatever file/directory
            # is referenced by path:
            base_path = os.path.dirname(path)

        # Load the scheme.yaml from the selected data directory:
        scheme = GroupAdditivityScheme.Load(os.path.join(base_path,
                                                         'scheme.yaml'))

        # Use that scheme to load the rest of the library:
        return cls._do_load(path, base_path, scheme)

    @classmethod
    def _Load(cls, path, scheme):
        if os.sep not in path and '.' not in path and not os.path.exists(path):
            # [JTF] where's our data directory?
            base_path = os.path.join(get_data_dir(), path)
            # We want to load the library.yaml in that directory:
            path = os.path.join(base_path, 'library.yaml')
        else:
            # The base path is the directory containing whatever file/directory
            # is referenced by path:
            base_path = os.path.dirname(path)

        # Use the scheme passed to us to load the rest of the library:
        return cls._do_load(path, base_path, scheme)

    @classmethod
    def _do_load(cls, path, base_path, scheme):
        # Read data from file.
        context = {'base_path': base_path}
        with open(path) as f:
            lib_data = yaml_io.load(
                yaml_io.parse(f.read()), context, loader=cls._yaml_loader)

        context['units'] = lib_data.units

        group_properties = lib_data.groups
        other_descriptor_properties = lib_data.other_descriptors
        UQ = lib_data.UQ

        if cls._property_set_group_yaml_types:
            # Prepare property_sets loader.
            property_sets_loader = yaml_io.make_object_loader(yaml_io.parse(
                '\n'.join(('%r:\n    type: %r\n    optional: true'
                          % (str(name),
                             str(cls._property_set_group_yaml_types[name])))
                          for name in cls._property_set_group_yaml_types)))
            # Read all properties.
            lib_contents = {}
            for name in group_properties:
                group = Group.parse(scheme, name)
                if group in lib_contents:
                    raise KeyError('Multiple definitions of group %s' % group)
                property_sets = yaml_io.load(
                    group_properties[name], context,
                    loader=property_sets_loader)
                lib_contents[group] = property_sets
            for name in other_descriptor_properties:
                descriptor = Descriptor(scheme, name)
                if descriptor in lib_contents:
                    raise KeyError('Multiple definitions of descriptor %s' %
                                   descriptor)
                property_sets = yaml_io.load(
                    other_descriptor_properties[name], context,
                    loader=property_sets_loader)
                lib_contents[descriptor] = property_sets

            # Read UQ data
            uq_contents = {}
            if UQ:
                uq_contents['RMSE'] = yaml_io.load(
                    UQ['RMSE'], context,
                    loader=property_sets_loader)
                uq_contents['descriptors'] = UQ['InvCovMat']['groups']
                uq_contents['mat'] = np.array(UQ['InvCovMat']['mat'])
                uq_contents['dof'] = UQ['DOF']
        else:
            # No property sets defined.
            warn('GroupLibrary.load(): No property sets defined.')
            lib_contents = {}
            uq_contents = {}

        new_lib = cls(scheme, lib_contents, uq_contents, path=path)
        # Update with included content.
        for include_path in lib_data.include:
            new_lib.Update(cls._Load(os.path.join(base_path,
                                                  include_path), scheme))
        return new_lib

    def Update(self, lib, overwrite=False):
        """Add complete contents of `lib` into this library.

        Parameters
        ----------
        lib : :class:`GroupLibrary`
            Library to import from.
        overwrite : bool
            If True, then existing data may be overwritten by data from `lib`.
        """
        for (group, other_property_sets) in list(lib.items()):
            if group not in self.contents:
                self.contents[group] = {}
            property_sets = self.contents[group]

            for name in other_property_sets:
                if name not in property_sets:
                    property_sets[name] = other_property_sets[name].copy()
                else:
                    property_sets[name].update(
                        other_property_sets[name], overwrite)
        # UQ stuff can only be loaded once
        if self.uq_contents and lib.uq_contents:
            raise ValueError('More than one uncertainty quantification',
                             'information provided')
        if not self.uq_contents:
            self.uq_contents = lib.uq_contents
    _yaml_loader = yaml_io.make_object_loader(yaml_io.parse("""
units:
    type: mapping
    default: {}

include:
    type: list
    item_type: string
    default: []

groups:
    type: mapping
    default: {}

other_descriptors:
    type: mapping
    default: {}

UQ:
    type: mapping
    default: {}
"""))
