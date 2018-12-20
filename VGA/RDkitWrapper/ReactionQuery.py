from ..Error import ReactionQueryError
from rdkit import Chem
from . MolQuery import MolQuery
from collections import defaultdict
from itertools import product as itpd
from rdkit.Chem import rdqueries
from operator import add
"""
Written by Geun Ho Gu 8/4/2016
ReactionQuery class is equipped with enhanced query comparison over
RDkit query comparison
This module uses a graph description language called R.I.N.G.
Rangarajan, S., Kaminski, T., Van Wyk, E., Bhan A., and Daoutidis, P.,
"Language-oriented rule-based reaction network generation and analysis:
    Algorithms of RING", Computers and Chemical Engineering 64 (2014)
    124, 10.1016/j.compchemeng.2014.02.007
Rangarajan, S., Bhan, A., and Daoutidis, P.,
"Language-oriented rule-based reaction network generation and analysis:
    Descrpition of RING", Computers and Chemical Engineering 45 (2012) 114,
    10.1016/j.compchemeng.2012.06.008
"""


class transformation(object):
    pass


class BondForm(transformation):
    def __init__(self, idx1, idx2, bondtype):
        self.idx1 = idx1
        self.idx2 = idx2
        self.bondtype = bondtype

    def __call__(self, comb_mol, mapped_index):
        comb_mol.AddBond(mapped_index[self.idx1],
                         mapped_index[self.idx2],
                         self.bondtype)


class BondBreak(transformation):
    def __init__(self, idx1, idx2):
        self.idx1 = idx1
        self.idx2 = idx2

    def __call__(self, comb_mol, mapped_index):
        comb_mol.RemoveBond(mapped_index[self.idx1], mapped_index[self.idx2])


class BondModify(transformation):
    def __init__(self, idx1, idx2, bondtype):
        self.idx1 = idx1
        self.idx2 = idx2
        self.bondtype = bondtype

    def __call__(self, comb_mol, mapped_index):
        comb_mol.RemoveBond(mapped_index[self.idx1], mapped_index[self.idx2])
        comb_mol.AddBond(mapped_index[self.idx1],
                         mapped_index[self.idx2], self.bondtype)


class BondIncrease(transformation):
    def __init__(self, idx1, idx2):
        self.idx1 = idx1
        self.idx2 = idx2

    def getbondtype(self, bondtype):
        if bondtype.GetBondType().__str__() == 'SINGLE':
            return Chem.BondType().DOUBLE
        elif bondtype.GetBondType().__str__() == 'DOUBLE':
            return Chem.BondType().TRIPLE
        elif bondtype.GetBondType().__str__() == 'TRIPLE':
            return Chem.BondType().QUADRUPLE
        elif bondtype.GetBondType().__str__() == 'QUADRUPLE':
            return Chem.BondType().QUINTUPLE
        elif bondtype.GetBondType().__str__() == 'QUINTUPLE':
            raise ReactionQueryError('BondIncrease: Increasing bond order',
                                     'above quintuple bond is not supported')
        elif bondtype.GetBondType().__str__() == 'AROMATIC':
            raise ReactionQueryError('BondIncrease: Increasing bond order',
                                     'of aromatic bond is not supported')
        elif bondtype.GetBondType().__str__() == 'DATIVE':
            raise ReactionQueryError('BondIncrease: Increasing bond order of',
                                     'partial bond is not supported')
        else:
            raise ReactionQueryError("BondDecrease: Unsupported bond type:'"
                                     + bondtype.GetBondType().__str__() + "'")

    def __call__(self, comb_mol, mapped_index):
        bondtype = self.getbondtype(comb_mol.
                                    GetBondBetweenAtoms(mapped_index
                                                        [self.idx1],
                                                        mapped_index
                                                        [self.idx2]))
        comb_mol.RemoveBond(mapped_index[self.idx1], mapped_index[self.idx2])
        comb_mol.AddBond(mapped_index[self.idx1],
                         mapped_index[self.idx2],
                         bondtype)


class AtomTypeModify(transformation):
    def __init__(self, idx, radical, charge, valence):
        self.idx = idx
        self.radical = radical
        self.charge = charge
        self.valence = valence

    def __call__(self, comb_mol, mapped_index):
        atom = comb_mol.GetAtomWithIdx(mapped_index[self.idx])
        if self.valence != 0:
            from rdkit.Chem import GetPeriodicTable
            atomicnum = atom.GetAtomicNum()
            atom = rdqueries.AtomNumEqualsQueryAtom(atomicnum)
            valence = GetPeriodicTable().GetDefaultValence(atomicnum)
            atom.ExpandQuery(rdqueries.
                             TotalValenceEqualsQueryAtom(valence+self.valence))
            atom.ExpandQuery(rdqueries.
                             FormalChargeEqualsQueryAtom(self.charge))
            comb_mol.ReplaceAtom(self.idx, atom)

        atom.SetNumRadicalElectrons(self.radical)
        atom.SetFormalCharge(self.charge)


class RadicalIncrease(transformation):
    def __init__(self, idx):
        self.idx = idx

    def __call__(self, comb_mol, mapped_index):
        atom = comb_mol.GetAtomWithIdx(mapped_index[self.idx])
        atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons()+1)


class RadicalDecrease(transformation):
    def __init__(self, idx):
        self.idx = idx

    def __call__(self, comb_mol, mapped_index):
        atom = comb_mol.GetAtomWithIdx(mapped_index[self.idx])
        atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons()-1)


class ChargeIncrease(transformation):
    def __init__(self, idx):
        self.idx = idx

    def __call__(self, comb_mol, mapped_index):
        atom = comb_mol.GetAtomWithIdx(mapped_index[self.idx])
        atom.SetFormalCharge(atom.GetFormalCharge()+1)


class ChargeDecrease(transformation):
    def __init__(self, idx):
        self.idx = idx

    def __call__(self, comb_mol, mapped_index):
        atom = comb_mol.GetAtomWithIdx(mapped_index[self.idx])
        atom.SetFormalCharge(atom.GetFormalCharge()-1)


class BondDecrease(transformation):
    def __init__(self, idx1, idx2):
        self.idx1 = idx1
        self.idx2 = idx2

    def getbondtype(self, bondtype):
        if bondtype.GetBondType().__str__() == 'SINGLE':
            return None
        elif bondtype.GetBondType().__str__() == 'DOUBLE':
            return Chem.BondType().SINGLE
        elif bondtype.GetBondType().__str__() == 'TRIPLE':
            return Chem.BondType().DOUBLE
        elif bondtype.GetBondType().__str__() == 'QUADRUPLE':
            return Chem.BondType().TRIPLE
        elif bondtype.GetBondType().__str__() == 'QUINTUPLE':
            return Chem.BondType().QUADRUPLE
        elif bondtype.GetBondType().__str__() == 'AROMATIC':
            raise ReactionQueryError('BondDecrease: Decreasing bond order',
                                     'of aromatic bond is not supported')
        elif bondtype.GetBondType().__str__() == 'DATIVE':
            raise ReactionQueryError('BondDecrease: Decreasing bond order',
                                     'of partial bond is not supported')
        else:
            raise ReactionQueryError("BondDecrease: Unsupported bond type:'"
                                     + bondtype.GetBondType().__str__() + "'")

    def __call__(self, comb_mol, mapped_index):
        bondtype = self.getbondtype(comb_mol.
                                    GetBondBetweenAtoms(mapped_index
                                                        [self.idx1],
                                                        mapped_index
                                                        [self.idx2]))
        comb_mol.RemoveBond(mapped_index[self.idx1],
                            mapped_index[self.idx2])
        if bondtype:
            comb_mol.AddBond(mapped_index[self.idx1],
                             mapped_index[self.idx2],
                             bondtype)


class ReactionQuery(object):
    def __init__(self):
        self.mol = Chem.RWMol(Chem.Mol())
        self.transformations = list()
        self.reactantquery = dict()
        self.constraints = defaultdict(list)
        self.atom_names = list()

    def __repr__(self):
        s = "%s(" % (self.__class__.__name__)
        if hasattr(self, 'name'):
            s += "name='%s', " % (self.name)
        s += "NumberReactants='%s, " % (len(self.reactants))
        nconstraint = len(self.mol_constraints)
        for constraint in self.atom_constraints:
            nconstraint += len(self.atom_constraints[constraint])
        s += "NumberConstraints=('%s')" % (nconstraint)
        s += ")"
        return s

    def AppendReactantQuery(self, reactant):
        if not isinstance(reactant, (Chem.Mol, MolQuery)):
            raise ReactionQueryError("Unrecognized instance'"
                                     + type(reactant)+"'")
        self.reactantquery[reactant.name] = reactant
        self.atom_names += reactant.atom_names

    def GetNumReactantTemplates(self):
        return len(self.reactantquery)

    def RunReactants(self, reactants):
        if isinstance(reactants, Chem.Mol):
            reactants = [reactants]
        elif isinstance(reactants, tuple):
            reactants = list(reactants)
        # Error test
        if len(reactants) != len(self.reactantquery):
            raise ReactionQueryError('Number of reactant mis-match')
        for i in range(0, len(reactants)):
            if not isinstance(reactants[i], Chem.Mol):
                raise ReactionQueryError("Unrecognized instance'"
                                         + type(reactants[i])+"'")
            reactants[i] = Chem.AddHs(reactants[i])
        # make matrix of matching index
        self.match_indexes = list()
        i = 0
        for reactant_name in self.reactantquery:
            if isinstance(self.reactantquery[reactant_name], Chem.Mol):
                match = reactants[i].GetSubstructMatches(self.reactantquery
                                                         [reactant_name])
                if not match:
                    return tuple()
                else:
                    self.match_indexes.append(match)
            elif isinstance(self.reactantquery[reactant_name], MolQuery):
                match = self.reactantquery[reactant_name].\
                    GetQueryMatches(reactants[i])
                if not match:
                    return tuple()
                else:
                    self.match_indexes.append(match)
            i += 1
        # make a index for combined reactants mol object
        self.combined_mol_match_index = list()
        for match in itpd(*self.match_indexes):
            row = list()
            modifier = 0
            for i in range(0, len(match)):
                if i != 0:
                    modifier = reactants[i-1].GetNumAtoms()
                row += list(map(add, list(match[i]), [modifier]*len(match[i])))
            self.combined_mol_match_index.append(row)
        # Combine reactants
        self.combined_mol = reactants[0]
        for i in range(1, len(reactants)):
            self.combined_mol = Chem.CombineMols(self.combined_mol,
                                                 reactants[i])

        self.combined_mol = Chem.RWMol(self.combined_mol)
        # Transform!!
        product_list = list()
        for matches in self.combined_mol_match_index:
            products = self.combined_mol.__copy__()
            for transform in self.transformations:
                transform(products, matches)
            products = Chem.GetMolFrags(products,
                                        asMols=True,
                                        sanitizeFrags=False)
            product_list.append(products)
        return tuple(product_list)

    def ModifyMolInPlace(self, reactants):
        # This applies transformation for all the matching pattern found in mol
        if len(self.reactantquery) != 1:
            raise ReactionQueryError('ModifyMolInPlace: Only usable for',
                                     'monomolecular transformation')
        if isinstance(reactants, Chem.Mol):
            reactants = [reactants]
        elif isinstance(reactants, tuple):
            reactants = list(reactants)
        # Error test
        if len(reactants) != len(self.reactantquery):
            raise ReactionQueryError('Number of reactant mis-match')
        for i in range(0, len(reactants)):
            if not isinstance(reactants[i], Chem.Mol):
                raise ReactionQueryError("Unrecognized instance'"
                                         + type(reactants[i])+"'")
            reactants[i] = Chem.AddHs(reactants[i])
        # make matrix of matching index
        self.match_indexes = list()
        i = 0
        for reactant_name in self.reactantquery:
            if isinstance(self.reactantquery[reactant_name], Chem.Mol):
                match = reactants[i].GetSubstructMatches(self.reactantquery
                                                         [reactant_name])
                if not match:
                    return tuple()
                else:
                    self.match_indexes.append(match)
            elif isinstance(self.reactantquery[reactant_name], MolQuery):
                match = self.reactantquery[reactant_name].\
                    GetQueryMatches(reactants[i])
                if not match:
                    return tuple()
                else:
                    self.match_indexes.append(match)
            i += 1
        # make a index for combined reactants mol object
        self.combined_mol_match_index = list()
        for match in itpd(*self.match_indexes):
            row = list()
            modifier = 0
            for i in range(0, len(match)):
                if i != 0:
                    modifier = reactants[i-1].GetNumAtoms()
                row += list(map(add,
                                list(match[i]),
                                [modifier]*len(match[i])))
            self.combined_mol_match_index.append(row)
        # Combine reactants
        self.combined_mol = reactants[0]
        for i in range(1, len(reactants)):
            self.combined_mol = Chem.CombineMols(self.combined_mol,
                                                 reactants[i])

        self.combined_mol = Chem.RWMol(self.combined_mol)
        # Transform!!        product_list = list()
        for matches in self.combined_mol_match_index:
            for transform in self.transformations:
                transform(self.combined_mol, matches)
        return self.combined_mol
