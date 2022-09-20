from ..Error import MolQueryError
from rdkit import Chem
from collections import defaultdict
import operator as op
"""
Written by Geun Ho Gu 8/4/2016
MolQuery class is equipped with enhanced query comparison over
RDkit query comparison.
This module uses a graph description language called R.I.N.G.
Rangarajan, S., Kaminski, T., Van Wyk, E., Bhan A., and Daoutidis, P.,
"Language-oriented rule-based reaction network generation and analysis:
    Algorithms of RING", Computers and Chemical Engineering 64 (2014) 124,
    10.1016/j.compchemeng.2014.02.007
Rangarajan, S., Bhan, A., and Daoutidis, P.,
"Language-oriented rule-based reaction network generation and analysis:
    Descrpition of RING", Computers and Chemical Engineering 45 (2012) 114,
    10.1016/j.compchemeng.2012.06.008
"""
ops = {'>': op.gt,
       '<': op.lt,
       '>=': op.ge,
       '<=': op.le,
       '=': op.eq}
digit = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '+', '-', 'E', 'e']


class ConstraintNumber(object):
    def __init__(self, s):
        if type(s) == list and len(s) == 2:  # given a list ['>=',2]
            self.operator = s[0]
            self.n = s[1]
        elif type(s) == list and len(s) == 1:
            self.operator = '='
            self.n = s[0]
        elif type(s) == str and s[0].isdigit():  # given a string number [3]
            self.operator = '='
            self.n = int(s)
        elif type(s) == int:
            self.operator = '='
            self.n = s
        else:  # given a string '>= 3'
            i = 1
            while not s[i] in digit and s[i] != ' ':
                i += 1
            self.operator = s[:i]
            self.n = int(s[i:])

    def __call__(self, inp):
        return ops[self.operator](inp, self.n)

    def __str__(self):
        return '<ConstraintNumber(%s%s)>' % (self.operator, self.n)

    def __repr__(self):
        return 'ConstraintNumber(%s%s)' % (self.operator, self.n)


class BondConstraint(object):
    def __init__(self, bondquery):
        assert isinstance(bondquery, BondQuery)
        self.bondquery = bondquery

    def __call__(self, idx1, idx2, mol):
        try:
            bond = mol.GetBondBetweenAtoms(idx1, idx2)
        except Exception:
            raise MolQueryError('There is no bond between atom'
                                + idx1+' and atom'+idx2)
        return self.bondquery(bond)

    def __str__(self):
        return '<BondConstraint(%s)>' % (self.bondquery)

    def __repr__(self):
        return 'BondConstraint(%s)' % (self.bondquery)

# BondQuery is not combined with BondConstraint because
# Bondquery is used elsewhere


class BondQuery(object):
    def __init__(self, RINGbondname):
        if RINGbondname not in ['single', 'double', 'triple',
                                'quadruple', 'ring', 'nonring',
                                'aromatic', 'any', 'strong', 'partial']:
            raise NotImplementedError("Unsupported bond type '"
                                      + RINGbondname + "'")
        self.RINGbondname = RINGbondname

    def __call__(self, rdkitbond):
        assert isinstance(rdkitbond, Chem.Bond)
        if self.RINGbondname == 'single':
            if rdkitbond.GetBondType() == Chem.BondType.SINGLE:
                return True
            else:
                raise MolQueryError('BondQuery: False.')
        elif self.RINGbondname == 'double':
            if rdkitbond.GetBondType() == Chem.BondType.DOUBLE:
                return True
            else:
                raise MolQueryError('BondQuery: False.')
        elif self.RINGbondname == 'triple':
            if rdkitbond.GetBondType() == Chem.BondType.TRIPLE:
                return True
            else:
                raise MolQueryError('BondQuery: False.')
        elif self.RINGbondname == 'quadruple':
            if rdkitbond.GetBondType() == Chem.BondType.QUADRUPLE:
                return True
            else:
                raise MolQueryError('BondQuery: False.')
        elif self.RINGbondname == 'ring':
            if rdkitbond.IsInRing():
                return True
            else:
                raise MolQueryError('BondQuery: False.')
        elif self.RINGbondname == 'nonring':
            if not rdkitbond.IsInRing():
                return True
            else:
                raise MolQueryError('BondQuery: False.')
        elif self.RINGbondname == 'aromatic':
            if rdkitbond.GetBondType() == Chem.BondType.AROMATIC:
                return True
            else:
                raise MolQueryError('BondQuery: False.')
        elif self.RINGbondname == 'any':
            return True
        elif self.RINGbondname == 'strong':
            if rdkitbond.GetBondType() in [Chem.BondType.DOUBLE,
                                           Chem.BondType.TRIPLE,
                                           Chem.BondType.QUADRUPLE,
                                           Chem.BondType.AROMATIC]:
                return True
            else:
                raise MolQueryError('BondQuery: False.')
        elif self.RINGbondname == 'partial':
            if rdkitbond.GetBondType() in [Chem.BondType.DATIVE,
                                           Chem.BondType.OTHER,
                                           Chem.BondType.ZERO]:
                return True
            else:
                raise MolQueryError('BondQuery: False.')
        else:
            raise NotImplementedError("Unsupported bond type '"
                                      + rdkitbond.GetBondType().name+"'")

    def __str__(self):
        return '<BondQuery(%s)>' % (self.RINGbondname)

    def __repr__(self):
        return 'BondQuery(%s)' % (self.RINGbondname)


class DoubleBondStereoConstraint(object):
    def __init__(self, negate, stereobondtype):
        self.negate = negate
        self.stereobondtype = stereobondtype

    def __call__(self, idx1, idx2, idx3, idx4, mol):
        double_bond = mol.GetBondBetweenAtoms(idx3, idx4)
        stereo = double_bond.GetStereo()
        # Chem.rdchem.BondStereo.STEREOZ = cis
        # Chem.rdchem.BondStereo.STEREOE = trans
        # Chem.rdchem.BondStereo.STEREONONE = Not specified
        # non-type check
        if stereo == Chem.rdchem.BondStereo.STEREONONE:
            if self.stereobondtype == stereo and not self.negate:
                return True
            elif self.stereobondtype != stereo and self.negate:
                return True
            else:
                raise MolQueryError('DoubleBondStereoConstraint: False.')
        # Cis trans. Here I use a slightly clever algorithm.
        # I look at the intersectio between the one specified by rdkit.
        nmatch = len(set(double_bond.GetStereoAtoms()).
                     intersection([idx1, idx2]))
        if nmatch == 2 or nmatch == 0:

            if self.stereobondtype == stereo and not self.negate:
                return True
            elif self.stereobondtype != stereo and self.negate:
                return True
            else:
                raise MolQueryError('DoubleBondStereoConstraint: False.')
        elif nmatch == 1:
            # if there is only one match, stereo is opposite.
            if stereo == Chem.rdchem.BondStereo.STEREOZ:
                stereo = Chem.rdchem.BondStereo.STEREOE
            elif stereo == Chem.rdchem.BondStereo.STEREOE:
                stereo = Chem.rdchem.BondStereo.STEREOZ
            if self.stereobondtype == stereo and not self.negate:
                return True
            elif self.stereobondtype != stereo and self.negate:
                return True
            else:
                raise MolQueryError('DoubleBondStereoConstraint: False.')

    def __str__(self):
        return '<DoubleBondStereoConstraint(%s%s)>' % (self.negate,
                                                       self.stereobondtype)

    def __repr__(self):
        return 'DoubleBondStereoConstraint(%s%s)' % (self.negate,
                                                     self.stereobondtype)


class AtomConstraint(object):
    pass


class AtomRadical(AtomConstraint):
    def __init__(self, negate, CN):
        self.negate = negate
        self.CN = CN

    def __call__(self, atom):
        if self.negate and self.CN(atom.GetNumRadicalElectrons()):
            raise MolQueryError('AtomRadical: False.')
        elif not self.negate and not self.CN(atom.GetNumRadicalElectrons()):
            raise MolQueryError('AtomRadical: False.')

    def __repr__(self):
        return 'AtomRadical()'


class AtomIsInRing(AtomConstraint):
    def __init__(self, negate):
        self.negate = negate

    def __call__(self, atom):
        if self.negate and atom.IsInRing():
            raise MolQueryError('AtomIsInRing: False.')
        elif not self.negate and not atom.IsInRing():
            raise MolQueryError('AtomIsInRing: False.')

    def __repr__(self):
        return 'AtomIsInRing()'


class AtomIsAromatic(AtomConstraint):
    def __init__(self, negate):
        self.negate = negate

    def __call__(self, atom):
        if self.negate and atom.GetIsAromatic():
            raise MolQueryError('AtomIsAromatic: False.')
        elif not self.negate and not atom.GetIsAromatic():
            raise MolQueryError('AtomIsAromatic: False.')

    def __repr__(self):
        return 'AtomAromatic()'


class AtomIsAllylic(AtomConstraint):
    def __init__(self, negate):
        self.negate = negate

    def __call__(self, atom):
        yes = False
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                yes = True
                break
        if self.negate and yes:
            raise MolQueryError('AtomAllylic: False.')
        elif not self.negate and not yes:
            raise MolQueryError('AtomAllylic: False.')

    def __repr__(self):
        return 'AtomAllylic()'


class AtomRing(AtomConstraint):
    def __init__(self, negate, ring_sizeCN):
        self.negate = negate
        self.ring_sizeCN = ring_sizeCN

    def __call__(self, atom):
        if self.negate:
            rings = atom.GetOwningMol().GetRingInfo().AtomRings()
            for ring in rings:
                if atom.GetIdx() in ring:
                    if self.ring_sizeCN(len(ring)):
                        raise MolQueryError('AtomRing: False.')
        else:
            if not atom.IsInRing():
                raise MolQueryError('AtomRing: False. There is no ring.')
            else:
                rings = atom.GetOwningMol().GetRingInfo().AtomRings()
                for ring in rings:
                    if atom.GetIdx() in ring:
                        if self.ring_sizeCN(len(ring)):
                            return
                raise MolQueryError('AtomRing: False.')

    def __repr__(self):
        return 'AtomRing()'


class AtomNRing(AtomConstraint):
    def __init__(self, negate, NringCN):
        self.negate = negate
        self.NringCN = NringCN

    def __call__(self, atom):
        # get number of rings
        rings = atom.GetOwningMol().GetRingInfo().AtomRings()
        n = 0
        for ring in rings:
            for ring_atom_idx in ring:
                if atom.GetIdx() == ring_atom_idx:
                    n += 1
        if self.negate and self.NringCN(n):
            raise MolQueryError('AtomNRing: False.')
        elif not self.NringCN(n):
            raise MolQueryError('AtomNRing: False.')

    def __repr__(self):
        return 'AtomNRing()'


class AtomConnectivityAtom(AtomConstraint):
    def __init__(self, negate, CN, connected, bondquery, constraints):
        assert isinstance(negate, bool)
        assert isinstance(CN, ConstraintNumber)
        self.negate = negate
        self.ConstraintNumber = CN
        self.connected = connected
        if isinstance(connected, str):
            NotImplementedError("GroupName in AtomConstraint Not supported")
        self.bondquery = bondquery
        self.constraints = constraints

    def __call__(self, atom):
        count = 0
        for bond in atom.GetBonds():
            # atom check
            connected_atom = bond.GetOtherAtom(atom)
            if not self.connected.Match(connected_atom):
                continue
            try:
                for constraint in self.constraints:
                    constraint(connected_atom)
            except Exception:
                continue
            # bond check
            try:
                if self.bondquery(bond):
                    count += 1
            except Exception:
                continue
        if self.negate:
            if self.ConstraintNumber(count):
                raise MolQueryError('AtomConnectivity: False. (Negate)',
                                    'Matching bond: %s. Constraint:%s.'
                                    % (count, self.ConstraintNumber))
        else:
            if not self.ConstraintNumber(count):
                raise MolQueryError('AtomConnectivity: False. Matching'
                                    'bond: %s. Constraint:%s.'
                                    % (count, self.ConstraintNumber))

    def __repr__(self):
        return 'AtomConnectivityAtom()'


class AtomConnectivityGroup(AtomConstraint):
    def __init__(self, negate, CN, connected, bondquery):
        assert isinstance(negate, bool)
        assert isinstance(CN, ConstraintNumber)
        self.negate = negate
        self.ConstraintNumber = CN
        self.connected = connected
        if isinstance(connected, str):
            NotImplementedError("GroupName in AtomConstraint Not supported")
        self.bondquery = bondquery

    def __call__(self, atom):
        count = 0
        # Do mol match first.
        mol = atom.GetOwningMol()
        for match_index in mol.GetSubstructMatches(self.connected):
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), match_index[0])
            try:
                if bond and self.bondquery(bond):
                    count += 1
            except Exception:
                continue
        if self.negate:
            if self.ConstraintNumber(count):
                raise MolQueryError('AtomConnectivity: False. (Negate)',
                                    'Matching bond: %s. Constraint:%s.'
                                    % (count, self.ConstraintNumber))
        else:
            if not self.ConstraintNumber(count):
                raise MolQueryError('AtomConnectivity: False. Matching bond:',
                                    '%s. Constraint:%s.'
                                    % (count, self.ConstraintNumber))

    def __repr__(self):
        return 'AtomConnectivityAtom()'


class MolConstraint(object):
    pass


class MolCharge(MolConstraint):
    def __init__(self, CN):
        assert isinstance(CN, ConstraintNumber)
        self.ConstraintNumber = CN

    def __call__(self, mol):
        total_charge = 0
        for atom in mol.GetAtoms():
            total_charge += atom.GetFormalCharge()
        if self.ConstraintNumber(total_charge):
            return self.ConstraintNumber(total_charge)
        else:
            raise MolQueryError('MolCharge: False. Total Charge: %s.',
                                'Constraint:%s.'
                                % (total_charge, self.ConstraintNumber))

    def __repr__(self):
        return 'MolCharge(%r)' % self.n


class MolAromatic(MolConstraint):
    def __call__(self, mol):
        if not mol.GetAromaticAtoms():
            raise MolQueryError("MolAromatic: mol not aromatic")

    def __repr__(self):
        return 'MolAromatic()'


class MolOlefinic(MolConstraint):
    def __call__(self, mol):
        if not mol.HasSubstructMatch(Chem.MolFromSmiles('C=C')):
            raise MolQueryError("MolOlefinic: No double bond")

    def __repr__(self):
        return 'MolOlefinic()'


class MolParaffinic(MolConstraint):
    def __call__(self, mol):
        if mol.HasSubstructMatch(Chem.MolFromSmiles('C=C')):
            raise MolQueryError("MolParaffinic: There is double bond")

    def __repr__(self):
        return 'MolParaffinic()'


class MolCyclic(MolConstraint):
    def __call__(self, mol):
        if not mol.GetRingInfo().NumRings():
            raise MolQueryError("MolCyclic: Mol isn't cyclic")

    def __repr__(self):
        return 'MolCyclic()'


class MolLinear(MolConstraint):
    def __call__(self, mol):
        if mol.GetRingInfo().NumRings():
            raise MolQueryError("MolLinear: Mol is cyclic")

    def __repr__(self):
        return 'MolCyclic()'


class MolQuery(object):
    """Goes through constraint and check for match"""
    def __init__(self):
        self.mol = Chem.RWMol(Chem.Mol())
        self.atom_names = list()
        self.mol_constraints = list()
        self.atom_constraints = defaultdict(list)
        """
        atom_constraints = { %atomidx% : [constraint1,constrain2],...}
        """
        self.bond_constraints = list()
        """
        bond_constraints = [[atomqueryidx1,atomqueryidx2,constraint],...]
        """
        self.double_bond_stereo_constraints = list()

    def __repr__(self):
        s = "%s(" % (self.__class__.__name__)
        if hasattr(self, 'name'):
            s += "name='%s', " % (self.name)
        s += "smiles='%s, " % Chem.MolToSmiles(self.mol)
        nconstraint = len(self.mol_constraints)
        for constraint in self.atom_constraints:
            nconstraint += len(self.atom_constraints[constraint])
        s += "NumberConstraints=('%s')" % (nconstraint)
        s += ")"
        return s

    def AppendMolConstraint(self, constraint):
        assert isinstance(constraint, MolConstraint)
        self.mol_constraints.append(constraint)

    def AppendAtomConstraint(self, constraint, idx_or_atomname):
        assert isinstance(constraint, AtomConstraint)
        if isinstance(idx_or_atomname, str):
            try:
                idx = self.atom_names.index(idx_or_atomname)
            except Exception:
                s = "Undeclared atom name: '" + idx_or_atomname + "'"
                raise MolQueryError(s)
        elif isinstance(idx_or_atomname, int):
            idx = idx_or_atomname
        else:
            MolQueryError('Unrecognized atom name instance')
        # If empty, initialize

        self.atom_constraints[idx].append(constraint)

    def AppendBondConstraint(self,
                             idx_or_atomname1,
                             idx_or_atomname2,
                             bondconstraint):
        assert isinstance(bondconstraint, BondConstraint)
        if isinstance(idx_or_atomname1, str):
            try:
                idx1 = self.atom_names.index(idx_or_atomname1)
            except Exception:
                s = "Undeclared atom name: '" + idx_or_atomname1 + "'"
                raise MolQueryError(s)
        elif isinstance(idx_or_atomname1, int):
            idx1 = idx_or_atomname1
        else:
            MolQueryError('Unrecognized atom name instance')
        if isinstance(idx_or_atomname2, str):
            try:
                idx2 = self.atom_names.index(idx_or_atomname2)
            except Exception:
                s = "Undeclared atom name: '" + idx_or_atomname2 + "'"
                raise MolQueryError(s)
        elif isinstance(idx_or_atomname2, int):
            idx2 = idx_or_atomname2
        else:
            MolQueryError('Unrecognized atom name instance')
        # If empty, initialize

        self.bond_constraints.append([idx1, idx2, bondconstraint])

    def AppendDoubleBondStereoConstraint(self,
                                         idx_or_atomname1,
                                         idx_or_atomname2,
                                         idx_or_atomname3,
                                         idx_or_atomname4,
                                         doublebondstereoconstraint):

        assert isinstance(doublebondstereoconstraint,
                          DoubleBondStereoConstraint)
        if isinstance(idx_or_atomname1, str):
            try:
                idx1 = self.atom_names.index(idx_or_atomname1)
            except Exception:
                s = "Undeclared atom name: '" + idx_or_atomname1 + "'"
                raise MolQueryError(s)
        elif isinstance(idx_or_atomname1, int):
            idx1 = idx_or_atomname1
        else:
            MolQueryError('Unrecognized atom name instance')
        if isinstance(idx_or_atomname2, str):
            try:
                idx2 = self.atom_names.index(idx_or_atomname2)
            except Exception:
                s = "Undeclared atom name: '" + idx_or_atomname2 + "'"
                raise MolQueryError(s)
        elif isinstance(idx_or_atomname2, int):
            idx2 = idx_or_atomname2
        else:
            MolQueryError('Unrecognized atom name instance')
        if isinstance(idx_or_atomname3, str):
            try:
                idx3 = self.atom_names.index(idx_or_atomname3)
            except Exception:
                s = "Undeclared atom name: '" + idx_or_atomname3 + "'"
                raise MolQueryError(s)
        elif isinstance(idx_or_atomname3, int):
            idx3 = idx_or_atomname3
        else:
            MolQueryError('Unrecognized atom name instance')
        if isinstance(idx_or_atomname4, str):
            try:
                idx4 = self.atom_names.index(idx_or_atomname4)
            except Exception:
                s = "Undeclared atom name: '" + idx_or_atomname4 + "'"
                raise MolQueryError(s)
        elif isinstance(idx_or_atomname4, int):
            idx4 = idx_or_atomname4
        else:
            MolQueryError('Unrecognized atom name instance')
        # If empty, initialize

        self.double_bond_stereo_constraints.append([idx1, idx2, idx3, idx4,
                                                    doublebondstereoconstraint]
                                                   )

    def GetQueryMatches(self, mol, debug=0):
        mol = Chem.AddHs(mol)
        # test mol constraints
        for mol_constraint in self.mol_constraints:
            try:
                mol_constraint(mol)
            except Exception:
                return tuple()
        # Match mol substructure matches
        maxMatches = 10000
        rdkit_matches = mol.GetSubstructMatches(self.mol, uniquify=False,
                                                maxMatches=maxMatches)
        lenMatches = len(rdkit_matches)
        if lenMatches == maxMatches:
            print('\nMax RDKit substructure matches exceeded. All groups in',
                  'molecule may not have been determined.')
            print(Chem.MolToSmiles(Chem.RemoveHs(mol)))
        if debug:
            print('structure matches:' + str(rdkit_matches))
        # if no rdkit match
        if not rdkit_matches:
            return tuple()
        # Enhanced bond matching
        matches1 = list()
        for match_indice in rdkit_matches:
            try:
                for bond_constraint in self.bond_constraints:
                    idx1 = match_indice[bond_constraint[0]]
                    idx2 = match_indice[bond_constraint[1]]
                    bond_constraint[2](idx1, idx2, mol)
                matches1.append(match_indice)
            except Exception:
                continue
        # Enhanced atom matching
        matches2 = list()
        for match_indice in matches1:
            try:
                for i in range(0, len(match_indice)):
                    if i in self.atom_constraints:
                        # extract atom
                        # atom_index = match_indice[i]
                        atom = mol.GetAtomWithIdx(match_indice[i])
                        for atom_constraint in self.atom_constraints[i]:
                            atom_constraint(atom)
                matches2.append(match_indice)
            except Exception:
                continue
        # Enhanced bond stereo matching
        matches3 = list()
        for match_indice in matches2:
            try:
                for double_bond_stereo_constraint in\
                     self.double_bond_stereo_constraints:
                    idx1 = match_indice[double_bond_stereo_constraint[0]]
                    idx2 = match_indice[double_bond_stereo_constraint[1]]
                    idx3 = match_indice[double_bond_stereo_constraint[2]]
                    idx4 = match_indice[double_bond_stereo_constraint[3]]
                    double_bond_stereo_constraint[4](idx1, idx2,
                                                     idx3, idx4, mol)
                matches3.append(match_indice)
            except Exception:
                continue
        return tuple(matches3)
