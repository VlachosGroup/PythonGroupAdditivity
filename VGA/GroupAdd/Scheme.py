from collections import defaultdict
from rdkit import Chem
from .. RINGParser import Read
import yaml
import os
from .. Error import PatternMatchError
from . Group import Group
from . DataDir import get_data_dir

# Scheme contains the group pattern information, and can assign groups based on the smiles.

class Scheme(object):
    pass

class GroupAdditivityScheme(Scheme):
    def __init__(self, patterns=[], pretreatment_rules=[], remaps ={}, other_descriptors=[],smiles_based_descriptors=[],smarts_based_descriptors=[], include =[]):
        """Load group-additivity scheme from file-system `path` or builtin.

        Parameters
        ----------
        path : str
            Specify either the path to a file containing the data or a symbolic
            name of a builtin scheme to load (*e.g.* ``benson`` to load the
            standard Benson groups scheme.)

        Returns
        -------
        scheme : :class:`GroupScheme`
            Group scheme containing the loaded data.

        Raises
        ------
        VGA.error.GroupSchemeError
            If `signature` is specified and does not match the signature of
            the loaded data.
        """
        self.patterns = patterns
        self.pretreatment_rules = pretreatment_rules
        self.remaps = remaps
        self.other_descriptors = other_descriptors
        self.smiles_based_descriptors = smiles_based_descriptors
        self.smarts_based_descriptors = smarts_based_descriptors
        for scheme in include:
            if isinstance(scheme,GroupAdditivityScheme):
                scheme_object = scheme
            else:
                scheme_object = self.Load(scheme['include'][0])
            self.patterns += scheme_object.patterns
            self.pretreatment_rules += scheme_object.pretreatment_rules
            self.remaps.update(scheme_object.remaps)
            self.other_descriptors += scheme_object.other_descriptors
            self.smiles_based_descriptors = scheme_object.smiles_based_descriptors
            self.smarts_based_descriptors = scheme_object.smarts_based_descriptors

    @classmethod
    def Load(cls, path):
        if os.sep not in path and '.' not in path and not os.path.exists(path):
            # [JTF] where's our data directory?
            base_path = os.path.join(get_data_dir(), path)
            # Load the scheme.yaml file in that directory:
            actual_path = os.path.join(base_path,'scheme.yaml')
        else:
            actual_path = path
        abs_path = os.path.abspath(actual_path)
        scheme_data = yaml.load(open(abs_path,'r'))
        # change to molquery
        for i in range(0,len(scheme_data['patterns'])):
            scheme_data['patterns'][i]['connectivity'] = \
                Read(scheme_data['patterns'][i]['connectivity'])

        patterns = scheme_data['patterns']
        pretreatment_rules=[]
        remaps ={}
        other_descriptors=[]
        smiles_based_descriptors=[]
        smarts_based_descriptors=[]
        if 'pretreatment_rules' in scheme_data:
            pretreatment_rules = scheme_data['pretreatment_rules']
        if 'remaps' in scheme_data:
            remaps = scheme_data['remaps']
        if 'other_descriptors' in scheme_data:
            for i in range(0,len(scheme_data['other_descriptors'])):
                scheme_data['other_descriptors'][i]['connectivity'] = \
                    Read(scheme_data['other_descriptors'][i]['connectivity'])
            other_descriptors = scheme_data['other_descriptors']
        if 'smiles_based_descriptors' in scheme_data:
            for i in range(0,len(scheme_data['smiles_based_descriptors'])):
                scheme_data['smiles_based_descriptors'][i]['smarts'] = \
                    Chem.MolFromSmarts(scheme_data['smiles_based_descriptors'][i]['smarts'])
            smiles_based_descriptors = scheme_data['smiles_based_descriptors']
        if 'smarts_based_descriptors' in scheme_data:
            for i in range(0,len(scheme_data['smarts_based_descriptors'])):
                scheme_data['smarts_based_descriptors'][i]['smarts'] = \
                    Chem.MolFromSmarts(scheme_data['smarts_based_descriptors'][i]['smarts'])
            smarts_based_descriptors = scheme_data['smarts_based_descriptors']
        return cls(patterns, pretreatment_rules, remaps, other_descriptors,smiles_based_descriptors,smarts_based_descriptors)

    def GetDescriptors(self, mol, debug=0):
        if isinstance(mol,Chem.Mol):

            mol = Chem.AddHs(mol)
            Chem.Kekulize(mol)
        elif isinstance(mol,str):
            clean_mol = Chem.MolFromSmiles(mol)
            # sanitize false is necessary to read adsorbed hydrogen species. (e.g. [Pt][H])
            mol.replace('{','')
            mol.replace('}','')
            #mol = Chem.MolFromSmiles(mol,sanitize=False) # Stereochemistry get's erased. :/
            mol = Chem.MolFromSmiles(mol)
            if mol is None:
                print('Smiles could not be loaded.')
                raise Exception
            sanitize_except_aromatization(mol)
            mol = Chem.AddHs(mol)
            Chem.Kekulize(mol)
        # pretreatment
        #TODO
        #reaction_query = Read(self.pretreatment_rules[0][0])
        #mol = reaction_query.ModifyMolInPlace(mol).
        # Change unspecified bond to ZERO, or weak bonds.
        for bond in mol.GetBonds():
            if bond.GetBondType().__str__()=='UNSPECIFIED':
                bond.SetBondType(Chem.BondType.ZERO)
        # aromatize C6 ring for benson
        _aromatization_Benson(mol)
        # assign groups
        self._AssignCenterPattern(mol,debug)
        groups = self._AssignGroup(mol)
        descriptors = self._AssignDescriptor(mol,clean_mol)
        all_descriptors = groups.copy()
        all_descriptors.update(descriptors)
        return all_descriptors

    def _AssignCenterPattern(self, mol,debug=0):

        for pattern in self.patterns:

            matches = pattern['connectivity'].GetQueryMatches(mol)
            # only the first atom is the center atom
            matches = set([match[0] for match in matches])
            for match in matches:
                atom = mol.GetAtomWithIdx(match)
                if atom.HasProp('Group_Center_Name'):
                    s = 'Group center overwritten:' + \
                        atom.GetProp('Group_Center_Name') + ' to ' + \
                        pattern['center_name']
                    s += '(' + pattern['connectivity'].__str__() + ')'
                    raise PatternMatchError(s,atom)
                else:
                    atom.SetProp('Group_Center_Name',pattern['center_name'])
                    atom.SetProp('Group_Periph_Name',pattern['periph_name'])
            if debug:
                s = '\nPattern: ' + pattern['center_name'] + ',' + pattern['periph_name']
                s += '\nMatches: ' + str(matches)
                print(s)
        # exception spitted if any atoms have no groups assigned
        for atom in mol.GetAtoms():
            if not atom.HasProp('Group_Center_Name'):
                s = '\nSub group not assigned:\n' + Chem.MolToSmiles(mol) +\
                    ' Atom symbol: ' + atom.GetSymbol() + \
                    ' #radicalE: ' + atom.GetNumRadicalElectrons().__str__() +\
                    '\nConnectivity:'
                for bond in atom.GetBonds():
                    if bond.GetBeginAtomIdx() == atom.GetIdx():
                        s += '\n' + bond.GetEndAtom().GetSymbol() + ' ' + \
                            bond.GetBondType().__str__()
                    else:
                        s += '\n' + bond.GetBeginAtom().GetSymbol() + ' ' + \
                            bond.GetBondType().__str__()
                s += '\n Occured '
                raise PatternMatchError(s,atom)

    def _AssignGroup(self,mol):
        groups = defaultdict(int)
        for atom in mol.GetAtoms():
            csg = atom.GetProp('Group_Center_Name')
            if csg != 'none':
                psgs = list()
                for neighboratom in atom.GetNeighbors():
                    psg = neighboratom.GetProp('Group_Periph_Name')
                    if psg != 'none':
                        psgs.append(psg)
                group = Group(self,csg,psgs)
                groups[group.name] +=1
                atom.SetProp('Group_name',group.name)
            else:
                atom.SetProp('Group_name','none')
        # remapping
        if hasattr(self,'remaps'):
            for atom in mol.GetAtoms():
                if atom.GetProp('Group_name') in self.remaps:
                    atom.SetProp('Group_name',self.remaps[atom.GetProp('Group_name')][0][1])
            for group in list(groups.keys()):
                if group in self.remaps:
                    n = groups.pop(group)
                    for remap in self.remaps[group]:
                        nn = n*remap[0]
                        groups[remap[1]] += nn
        return groups

    def _AssignDescriptor(self,mol,clean_mol):
        descriptors = defaultdict(int)
        for descriptor in self.other_descriptors:
            matches = descriptor['connectivity'].GetQueryMatches(mol)
            matches = set([tuple(set(match)) for match in matches])
            if matches:
                descriptors[descriptor['name']] += len(matches)
        for descriptor in self.smiles_based_descriptors:
            matches = clean_mol.GetSubstructMatches(descriptor['smiles'],useChirality=descriptor['useChirality'])
            matches = set([tuple(set(match)) for match in matches])
            if matches:
                descriptors[descriptor['name']] += len(matches)
        for descriptor in self.smarts_based_descriptors:
            matches = mol.GetSubstructMatches(descriptor['smarts'],useChirality=descriptor['useChirality'])
            matches = set([tuple(set(match)) for match in matches])
            if matches:
                descriptors[descriptor['name']] += len(matches)
        # remaps
        if hasattr(self,'remaps'):
            for descriptor in list(descriptors.keys()):
                if descriptor in self.remaps:
                    n = descriptors.pop(descriptor)
                    for remap in self.remaps[descriptor]:
                        nn = n*remap[0]
                        descriptors[remap[1]] += nn
        return descriptors


def sanitize_except_aromatization(mol):
    Chem.SanitizeMol(mol,sanitizeOps= Chem.rdmolops.SanitizeFlags.SANITIZE_ADJUSTHS)
    Chem.SanitizeMol(mol,sanitizeOps= Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUP)
    Chem.SanitizeMol(mol,sanitizeOps= Chem.rdmolops.SanitizeFlags.SANITIZE_CLEANUPCHIRALITY)
    Chem.SanitizeMol(mol,sanitizeOps= Chem.rdmolops.SanitizeFlags.SANITIZE_FINDRADICALS)
    Chem.SanitizeMol(mol,sanitizeOps= Chem.rdmolops.SanitizeFlags.SANITIZE_KEKULIZE)
    Chem.SanitizeMol(mol,sanitizeOps= Chem.rdmolops.SanitizeFlags.SANITIZE_PROPERTIES)
    Chem.SanitizeMol(mol,sanitizeOps= Chem.rdmolops.SanitizeFlags.SANITIZE_SETCONJUGATION)
    Chem.SanitizeMol(mol,sanitizeOps= Chem.rdmolops.SanitizeFlags.SANITIZE_SETHYBRIDIZATION)
    Chem.SanitizeMol(mol,sanitizeOps= Chem.rdmolops.SanitizeFlags.SANITIZE_SYMMRINGS)

def _aromatization_Benson(mol):
    # Benson's aromatic bond criteria:
    # C6 ring
    # No oxygen atoms
    # No triple bond (radical is okay)
    # Evenly spaced 3 double bond
    rings = Chem.GetSymmSSSR(mol)
    for i in range(0,len(rings)):
        # ring size check
        if len(rings[i]) != 6:
            continue
        # atom check
        if mol.GetAtomWithIdx(rings[i][0]).GetSymbol() != 'C':
            continue
        if mol.GetAtomWithIdx(rings[i][1]).GetSymbol() != 'C':
            continue
        if mol.GetAtomWithIdx(rings[i][2]).GetSymbol() != 'C':
            continue
        if mol.GetAtomWithIdx(rings[i][3]).GetSymbol() != 'C':
            continue
        if mol.GetAtomWithIdx(rings[i][4]).GetSymbol() != 'C':
            continue
        if mol.GetAtomWithIdx(rings[i][5]).GetSymbol() != 'C':
            continue

        # bond check
        if mol.GetBondBetweenAtoms(rings[i][0],rings[i][1]).GetBondType().__str__() == 'SINGLE':
            if mol.GetBondBetweenAtoms(rings[i][1],rings[i][2]).GetBondType().__str__() == 'DOUBLE':
                if mol.GetBondBetweenAtoms(rings[i][2],rings[i][3]).GetBondType().__str__() == 'SINGLE':
                    if mol.GetBondBetweenAtoms(rings[i][3],rings[i][4]).GetBondType().__str__() == 'DOUBLE':
                        if mol.GetBondBetweenAtoms(rings[i][4],rings[i][5]).GetBondType().__str__() == 'SINGLE':
                            if mol.GetBondBetweenAtoms(rings[i][5],rings[i][0]).GetBondType().__str__() == 'DOUBLE':
                                pass
                            else:
                                continue
                        else:
                            continue
                    else:
                        continue
                else:
                    continue
            else:
                continue
        elif mol.GetBondBetweenAtoms(rings[i][0],rings[i][1]).GetBondType().__str__() == 'DOUBLE':
            if mol.GetBondBetweenAtoms(rings[i][1],rings[i][2]).GetBondType().__str__() == 'SINGLE':
                if mol.GetBondBetweenAtoms(rings[i][2],rings[i][3]).GetBondType().__str__() == 'DOUBLE':
                    if mol.GetBondBetweenAtoms(rings[i][3],rings[i][4]).GetBondType().__str__() == 'SINGLE':
                        if mol.GetBondBetweenAtoms(rings[i][4],rings[i][5]).GetBondType().__str__() == 'DOUBLE':
                            if mol.GetBondBetweenAtoms(rings[i][5],rings[i][0]).GetBondType().__str__() == 'SINGLE':
                                pass
                            else:
                                continue
                        else:
                            continue
                    else:
                        continue
                else:
                    continue
            else:
                continue
        else:
            continue
        # if it survived all these continue, then this molecule is aromatic
        mol.GetAtomWithIdx(rings[i][0]).SetIsAromatic(True)
        mol.GetAtomWithIdx(rings[i][1]).SetIsAromatic(True)
        mol.GetAtomWithIdx(rings[i][2]).SetIsAromatic(True)
        mol.GetAtomWithIdx(rings[i][3]).SetIsAromatic(True)
        mol.GetAtomWithIdx(rings[i][4]).SetIsAromatic(True)
        mol.GetAtomWithIdx(rings[i][5]).SetIsAromatic(True)
        mol.GetBondBetweenAtoms(rings[i][0],rings[i][1]).SetIsAromatic(True)
        mol.GetBondBetweenAtoms(rings[i][1],rings[i][2]).SetIsAromatic(True)
        mol.GetBondBetweenAtoms(rings[i][2],rings[i][3]).SetIsAromatic(True)
        mol.GetBondBetweenAtoms(rings[i][3],rings[i][4]).SetIsAromatic(True)
        mol.GetBondBetweenAtoms(rings[i][4],rings[i][5]).SetIsAromatic(True)
        mol.GetBondBetweenAtoms(rings[i][5],rings[i][0]).SetIsAromatic(True)
        mol.GetBondBetweenAtoms(rings[i][0],rings[i][1]).SetIsConjugated(True)
        mol.GetBondBetweenAtoms(rings[i][1],rings[i][2]).SetIsConjugated(True)
        mol.GetBondBetweenAtoms(rings[i][2],rings[i][3]).SetIsConjugated(True)
        mol.GetBondBetweenAtoms(rings[i][3],rings[i][4]).SetIsConjugated(True)
        mol.GetBondBetweenAtoms(rings[i][4],rings[i][5]).SetIsConjugated(True)
        mol.GetBondBetweenAtoms(rings[i][5],rings[i][0]).SetIsConjugated(True)
        mol.GetBondBetweenAtoms(rings[i][0],rings[i][1]).SetBondType(Chem.rdchem.BondType.AROMATIC)
        mol.GetBondBetweenAtoms(rings[i][1],rings[i][2]).SetBondType(Chem.rdchem.BondType.AROMATIC)
        mol.GetBondBetweenAtoms(rings[i][2],rings[i][3]).SetBondType(Chem.rdchem.BondType.AROMATIC)
        mol.GetBondBetweenAtoms(rings[i][3],rings[i][4]).SetBondType(Chem.rdchem.BondType.AROMATIC)
        mol.GetBondBetweenAtoms(rings[i][4],rings[i][5]).SetBondType(Chem.rdchem.BondType.AROMATIC)
        mol.GetBondBetweenAtoms(rings[i][5],rings[i][0]).SetBondType(Chem.rdchem.BondType.AROMATIC)
