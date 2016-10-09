import os
from warnings import warn
from collections import Mapping
import yaml

from .. Error import GroupMissingDataError
from . Group import Group
from . Scheme import GroupAdditivityScheme

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

    In addition to the above methods, the :class:`GroupLibrary()` type may be
    treated as a mapping from :class:`Group()` to a secondary mapping.  The
    secondary mapping is from *property set* `name` to an object that defines
    the contributions of those properties for the group.  The type of this
    object depends on the particular *property set* requested.
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
            raise KeyError('Property set %r already registered.'%name)
        cls._property_set_group_yaml_types[name] = group_yaml_type
        cls._property_set_estimator_types[name] = estimator_type

    def __init__(self, scheme, contents={}, path=None):
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
            contents = contents.items()
        self.contents = dict((group, property_sets)
            for (group, property_sets) in contents)

    def match_groups(self, chem, manual_descriptors={}):
        """Determine groups appearing in chemical structure `chem`.

        Parameters
        ----------
        chem : :class:`chemtk.structure.Structure`
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
        return self.scheme.match_groups(chem, manual_descriptors)

    def estimate(self, groups, property_set_name):
        """Estimate set of properties for chemical.

        Parameters
        ----------
        groups : mapping
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
            raise KeyError('Invalid property_set name: %r'%property_set_name)
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
        """(class method) Load group-additivity library from file-system `path` or builtin.

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
        if '/' not in path and '.' not in path and not os.path.exists(path):
            # Look for folder with name, benson.
            base_path = os.path.join(os.path.split(os.path.split(os.path.split(__file__)[0])[0])[0],'data')
            base_path = os.path.join(base_path,path)
            path = os.path.join(base_path,'library.yaml')
        else:
            base_path = os.path.split(path)[0]

        return cls._DoLoad(path, base_path)

    @classmethod
    def _DoLoad(cls, path, base_path):
        # Read data from file.
        context = {'base_path': base_path}
        print path
        with open(path) as f:
            lib_data = yaml.load(f)

        context['units'] = lib_data.get('units')
        scheme = lib_data.get('scheme')
        group_properties = lib_data.get('contents')

        if cls._property_set_group_yaml_types:
            # Prepare property_sets loader.
            property_sets_loader = yaml.load(
                '\n'.join(('%r:\n    type: %r\n    optional: true'
                        %(str(name),
                            str(cls._property_set_group_yaml_types[name])))
                    for name in cls._property_set_group_yaml_types))

            # Read all properties.
            lib_contents = {}
            for name in group_properties:
                group = Group.parse(scheme, name)
                if group in lib_contents:
                    raise KeyError('Multiple definitions of group %s'%group)
                property_sets = yaml.load(
                    group_properties[name], context,
	                loader=property_sets_loader)
                lib_contents[group] = property_sets
        else:
            # No property sets defined.
            lib_contents = {}

        new_lib = cls(scheme, lib_contents, path=path)
        # Update with included content.
        if lib_data.has_key('include'):
            for include_path in lib_data.get('include'):
                new_lib.Update(cls.Load(os.path.join(base_path, include_path)))
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
        self.scheme = GroupAdditivityScheme(include=[self.scheme, lib.scheme])

        for (group, other_property_sets) in lib.items():
            if group not in self.contents:
                self.contents[group] = {}
            property_sets = self.contents[group]

            for name in other_property_sets:
                if name not in property_sets:
                    property_sets[name] = other_property_sets[name].copy()
                else:
                    property_sets[name].update(
                        other_property_sets[name], overwrite)