from .. Error import RINGReaderError
from .. RDkitWrapper.ReactionQuery import ReactionQuery, BondForm, \
    BondIncrease, BondDecrease, BondModify, AtomTypeModify, BondBreak,\
    RadicalIncrease, RadicalDecrease, ChargeIncrease, ChargeDecrease
from .MolQueryRead import MolQueryReader
from rdkit import Chem


class ReactionQueryReader(object):
    """
    creates enhanced reaction query that is based on RING language

    from pgradd.RINGParser.Reader import Read
    from rdkit import Chem
    s = "\
    rule increaseBO{     \
    reactant r1{     \
    C labeled c1     \
    H labeled h1 single bond to c1    \
    }     \
    increase number of radical (c1)     \
    increase number of radical (h1)     \
    break bond(c1,h1) }"
    rxnquery = Read(s)
    mol = Chem.MolFromSmiles('CC')
    products = rxnquery.RunReactants(mol)
    for product in products:
    for p in product:
        print Chem.MolToSmiles(p)
    """

    def __init__(self, tree, RINGgroups=None):
        # ast = Abstract Syntax Tree
        self.tree = tree
        self.RINGgroups = RINGgroups
        self.atom_names = list()
        self.atom_belonging_mol = list()
        self.electronbalance = list()

    def ReadBondType(self, tree):
        if tree[0] == 'single':
            bondtype = Chem.BondType.SINGLE
            balance = 1
        elif tree[0] == 'double':
            bondtype = Chem.BondType.DOUBLE
            balance = 2
        elif tree[0] == 'triple':
            bondtype = Chem.BondType.TRIPLE
            balance = 3
        elif tree[0] == 'quadruple':
            bondtype = Chem.BondType.QUADRUPLE
            balance = 4
        elif tree[0] == 'aromatic':
            bondtype = Chem.BondType.AROMATIC
            balance = 1.5
        elif tree[0] == 'partial':
            bondtype = Chem.BondType.DATIVE
            balance = 0
        else:
            raise RINGReaderError("BondType: Unsupported bond type: '"
                                  + tree[0] + "'")
        return bondtype, balance

    def ReadBondForm(self, tree, reactionquery):
        i = 0
        if tree[i][0] == 'BondType':
            bondtype, balance = self.ReadBondType(tree[i][1:])
            i += 1
        else:
            bondtype = Chem.BondType.SINGLE
            balance = 1
        try:
            idx1 = self.atom_names.index(tree[i][1])
            i += 1
        except Exception:
            raise RINGReaderError("BondForm: Unrecognized atom label: '"
                                  + tree[i][1] + "'")
        try:
            idx2 = self.atom_names.index(tree[i][1])
        except Exception:
            raise RINGReaderError("BondForm: Unrecognized atom label: '"
                                  + tree[i][1] + "'")
        self.electronbalance[idx1] += -balance
        self.electronbalance[idx2] += -balance
        reactionquery.transformations.append(BondForm(idx1, idx2, bondtype))

    def ReadBondBreak(self, tree, reactionquery):
        i = 0
        if tree[i][0] == 'BondType':
            bondtype, balance = self.ReadBondType(tree[i][1:])
            i += 1
        else:
            bondtype = Chem.BondType.SINGLE
            balance = +1

        try:
            atom_name1 = tree[i][1]
            idx1 = self.atom_names.index(atom_name1)
            reactant_name1 = self.atom_belonging_mol[idx1]
            idx1_in_query = reactionquery.reactantquery[reactant_name1].\
                atom_names.index(atom_name1)
            i += 1
        except Exception:
            raise RINGReaderError("BondBreak: Unrecognized atom label: '"
                                  + tree[i][1] + "'")
        try:
            atom_name2 = tree[i][1]
            idx2 = self.atom_names.index(atom_name2)
            reactant_name2 = self.atom_belonging_mol[idx2]
            idx2_in_query = reactionquery.reactantquery[reactant_name2].\
                atom_names.index(atom_name2)
        except Exception:
            raise RINGReaderError("BondBreak: Unrecognized atom label: '"
                                  + tree[i][1] + "'")
        if reactant_name1 != reactant_name2:
            raise RINGReaderError("BondBreak: Two supplied atom labels aren't",
                                  "in the same molecule")
        bond = reactionquery.reactantquery[reactant_name1].mol.\
            GetBondBetweenAtoms(idx1_in_query, idx2_in_query)
        if not bond:
            raise RINGReaderError("BondBreak: No bond between two atom ("
                                  + atom_name1 + "," + atom_name2 + ")")
        elif bond.GetBondType().__str__() == 'UNSPECIFIED':
            raise RINGReaderError("BondBreak: Break of unspecified bond",
                                  "is not supported ")
        elif bond.GetBondType() != bondtype:
            raise RINGReaderError("BondBreak: Bond mismatch between two",
                                  "atom (" + atom_name1 + "," + atom_name2
                                  + "). (Bond in mol:"
                                  + bond.GetBondType().__str__() +
                                  " ,bond broken:" + bondtype.__str__())

        self.electronbalance[idx1] += balance
        self.electronbalance[idx2] += balance
        reactionquery.transformations.append(BondBreak(idx1, idx2))

    def ReadBondModify(self, tree, reactionquery):
        try:
            atom_name1 = tree[0][1]
            idx1 = self.atom_names.index(atom_name1)
            reactant_name1 = self.atom_belonging_mol[idx1]
            idx1_in_query = reactionquery.reactantquery[reactant_name1].\
                atom_names.index(atom_name1)
        except Exception:
            raise RINGReaderError("BondModify: Unrecognized atom label: '"
                                  + tree[0][1] + "'")
        try:
            atom_name2 = tree[1][1]
            idx2 = self.atom_names.index(atom_name2)
            reactant_name2 = self.atom_belonging_mol[idx2]
            idx2_in_query = reactionquery.reactantquery[reactant_name2].\
                atom_names.index(atom_name2)
        except Exception:
            raise RINGReaderError("BondModify: Unrecognized atom label: '"
                                  + tree[1][1] + "'")
        if reactant_name1 != reactant_name2:
            raise RINGReaderError("BondModify: Two supplied atom labels",
                                  "aren't in the same molecule")
        bond = reactionquery.reactantquery[reactant_name1].\
            mol.GetBondBetweenAtoms(idx1_in_query, idx2_in_query)
        if not bond:
            raise RINGReaderError("BondModify: No bond between two atom ("
                                  + atom_name1 + "," + atom_name2 + ")")
        elif bond.GetBondType().__str__() == 'UNSPECIFIED':
            raise RINGReaderError("BondModify: BondModify: modification of",
                                  "unspecified bond is not supported ")
        existingbond = bond.GetBondType().__str__()
        if existingbond == 'SINGLE':
            balance1 = 1
        elif existingbond == 'DOUBLE':
            balance1 = 2
        elif existingbond == 'TRIPLE':
            balance1 = 3
        elif existingbond == 'QUADRUPLE':
            balance1 = 4
        elif existingbond == 'AROMATIC':
            balance1 = 1.5
        elif existingbond == 'DATIVE':
            balance1 = 0
        else:
            raise RINGReaderError("BondModify: Unsupported bond type for",
                                  "modifiation: '"+existingbond+"'")
        bondtype, balance = self.ReadBondType(tree[2][1:])
        self.electronbalance[idx1] -= balance - balance1
        self.electronbalance[idx2] -= balance - balance1
        reactionquery.transformations.append(BondModify(idx1, idx2, bondtype))

    def ReadBondIncrease(self, tree, reactionquery):
        try:
            idx1 = self.atom_names.index(tree[0][1])
        except Exception:
            raise RINGReaderError("BondIncrease: Unrecognized atom label: '"
                                  + tree[0][1] + "'")
        try:
            idx2 = self.atom_names.index(tree[1][1])
        except Exception:
            raise RINGReaderError("BondIncrease: Unrecognized atom label: '"
                                  + tree[1][1] + "'")

        self.electronbalance[idx1] += -1
        self.electronbalance[idx2] += -1
        reactionquery.transformations.append(BondIncrease(idx1, idx2))

    def ReadBondDecrease(self, tree, reactionquery):
        try:
            idx1 = self.atom_names.index(tree[0][1])
        except Exception:
            raise RINGReaderError("BondDecrease: Unrecognized atom label: '"
                                  + tree[0][1] + "'")
        try:
            idx2 = self.atom_names.index(tree[1][1])
        except Exception:
            raise RINGReaderError("BondDecrease: Unrecognized atom label: '"
                                  + tree[1][1] + "'")

        self.electronbalance[idx1] += 1
        self.electronbalance[idx2] += 1
        reactionquery.transformations.append(BondDecrease(idx1, idx2))

    def ReadAtomSuffix(self, tree):
        # '+','-','.',':','+.','-.','*'
        if tree[0] == '+.':
            radical = 1
            charge = 1
            valence = 0
        elif tree[0] == '-.':
            radical = 1
            charge = -1
            valence = 0
        elif tree[0] == '+':
            radical = 0
            charge = 1
            valence = 0
        elif tree[0] == '-':
            radical = 0
            charge = -1
            valence = 0
        elif tree[0] == '.':
            radical = 1
            charge = 0
            valence = 0
        elif tree[0] == ':':
            radical = 2
            charge = 0
            valence = 0
        elif tree[0] == '*':
            radical = 0
            valence = 1
            charge = 1
        else:
            s = "AtomSuffic: Unsupported atom suffic: '" + tree[0] + "'"
            raise NotImplementedError(s)
        return radical, charge, valence

    def ReadAtomType(self, tree):
        assert tree[0][0] == 'Symbols'
        symbol = tree[0][1][0]

        if len(tree) > 1:
            assert tree[1][0] == 'AtomSuffix'
            radical, charge, valence = self.ReadAtomSuffix(tree[1][1:])
        return symbol, radical, charge, valence

    def ReadAtomLabel(self, tree, reactionquery):
        try:
            atom_name = tree[0]
            idx = self.atom_names.index(atom_name)
            reactant_name = self.atom_belonging_mol[idx]
            idx_in_query = reactionquery.reactantquery[reactant_name].\
                atom_names.index(atom_name)
            atom = reactionquery.reactantquery[reactant_name].\
                mol.GetAtomWithIdx(idx_in_query)
        except Exception:
            raise RINGReaderError("AtomLabel: Unrecognized atom label: '"
                                  + tree[0] + "'")
        return atom_name, idx, reactant_name, idx_in_query, atom

    def ReadAtomTypeModify(self, tree, reactionquery):
        assert tree[0][0] == 'AtomLabel'
        _, idx, _, _, atom = self.ReadAtomLabel(tree[0][1:], reactionquery)
        assert tree[1][0] == 'AtomType'
        symbol, radical, charge, valence = self.ReadAtomType(tree[1][1:])
        if atom.GetSymbol() != symbol:
            raise NotImplementedError("AtomTypeModify: Atom Label change",
                                      "not supported")

        self.electronbalance[idx] -= radical - atom.GetNumRadicalElectrons()
        self.electronbalance[idx] += charge - atom.GetFormalCharge()
        reactionquery.transformations.append(AtomTypeModify(idx,
                                                            radical,
                                                            charge,
                                                            valence))

    def ReadRadicalModify(self, tree, reactionquery):
        assert tree[0][0] == 'AtomLabel'
        _, idx, _, _, atom = self.ReadAtomLabel(tree[0][1:], reactionquery)
        radical = tree[1]
        if radical < 0:
            raise RINGReaderError("RadicalModify: Number of radical",
                                  "electrons cannot be below 0")
        self.electronbalance[idx] -= radical - atom.GetNumRadicalElectrons()
        reactionquery.transformations.append(AtomTypeModify(idx,
                                                            radical,
                                                            0, 0))

    def ReadRadicalIncrease(self, tree, reactionquery):
        assert tree[0][0] == 'AtomLabel'
        _, idx, _, _, atom = self.ReadAtomLabel(tree[0][1:], reactionquery)

        self.electronbalance[idx] -= 1
        reactionquery.transformations.append(RadicalIncrease(idx))

    def ReadRadicalDecrease(self, tree, reactionquery):
        assert tree[0][0] == 'AtomLabel'
        atom_name, idx, _, _, atom = self.ReadAtomLabel(tree[0][1:],
                                                        reactionquery)
        self.electronbalance[idx] += 1
        reactionquery.transformations.append(RadicalDecrease(idx))

    def ReadChargeIncrease(self, tree, reactionquery):
        assert tree[0][0] == 'AtomLabel'
        _, idx, _, _, atom = self.ReadAtomLabel(tree[0][1:], reactionquery)

        self.electronbalance[idx] -= 1
        reactionquery.transformations.append(ChargeIncrease(idx))

    def ReadChargeDecrease(self, tree, reactionquery):
        assert tree[0][0] == 'AtomLabel'
        atom_name, idx, _, _, atom = self.ReadAtomLabel(tree[0][1:],
                                                        reactionquery)
        self.electronbalance[idx] += 1
        reactionquery.transformations.append(ChargeDecrease(idx))

    def ReadConnectivityChange(self, tree, reactionquery):
        if tree[0][0] == 'BondForm':
            self.ReadBondForm(tree[0][1:], reactionquery)
        elif tree[0][0] == 'BondBreak':
            self.ReadBondBreak(tree[0][1:], reactionquery)
        elif tree[0][0] == 'BondModify':
            self.ReadBondModify(tree[0][1:], reactionquery)
        elif tree[0][0] == 'BondIncrease':
            self.ReadBondIncrease(tree[0][1:], reactionquery)
        elif tree[0][0] == 'BondDecrease':
            self.ReadBondDecrease(tree[0][1:], reactionquery)
        elif tree[0][0] == 'AtomTypeModify':
            self.ReadAtomTypeModify(tree[0][1:], reactionquery)
        elif tree[0][0] == 'RadicalModify':
            self.ReadRadicalModify(tree[0][1:], reactionquery)
        elif tree[0][0] == 'RadicalIncrease':
            self.ReadRadicalIncrease(tree[0][1:], reactionquery)
        elif tree[0][0] == 'RadicalDecrease':
            self.ReadRadicalDecrease(tree[0][1:], reactionquery)
        elif tree[0][0] == 'ChargeIncrease':
            self.ReadChargeIncrease(tree[0][1:], reactionquery)
        elif tree[0][0] == 'ChargeDecrease':
            self.ReadChargeDecrease(tree[0][1:], reactionquery)

    def ReadTransformationChain(self, tree, reactionquery):
        assert tree[0][0] == 'ConnectivityChange'
        self.ReadConnectivityChange(tree[0][1:], reactionquery)
        if 1 < len(tree):
            assert tree[1][0] == 'TransformationChain'
            self.ReadTransformationChain(tree[1][1:], reactionquery)

    def LabelMapping(self, tree, labelmapping=None):
        if not labelmapping:
            labelmapping = dict()
        assert tree[0][0] == 'AtomLabel'
        assert tree[1][0] == 'AtomLabel'
        labelmapping[tree[0][1]] = tree[1][1]
        if 2 < len(tree):
            assert tree[2][0] == 'LabelMapping'
            labelmapping = self.LabelMapping(tree[2][1:], labelmapping)
        return labelmapping

    def ReadReactantGroup(self, tree, reactionquery):
        assert tree[0][0] == 'ReactantName'
        assert tree[1][0] == 'GroupName'
        assert tree[2][0] == 'LabelMapping'
        labelmapping = self.LabelMapping(tree[2][1:])
        if tree[1][1] not in self.RINGgroups:
            raise RINGReaderError("ReactantGroup: Unrecognized group name:'"
                                  + tree[1][1] + "'")

        if len(labelmapping) != len(self.RINGgroups[tree[1][1]]):
            raise RINGReaderError("ReactantGroup: Label mapping length",
                                  "doesn''t match length of atoms in group :'"
                                  + tree[1][1] + "'")
        reactionquery.reactantquery[tree[0][1]] = self.RINGgroups[tree[1][1]]
        reactionquery.reactantquery[tree[0][1]].name = tree[0][1]
        for i in range(0, len(reactionquery.reactantquery[tree[0][1]].
                              atom_names)):
            try:
                reactionquery.reactantquery[tree[0][1]].atom_names[i] = \
                    labelmapping[reactionquery.reactantquery[tree[0][1]].
                                 atom_names[i]]
            except KeyError:
                s = 'Unrecognized label '
                s += reactionquery.reactantquery[tree[0][1]].atom_names[i]
                s += ' in reactant ' +\
                    reactionquery.reactantquery[tree[0][1]].name
                s += ' for duplication opreation for group ' + tree[1][1]
                raise RINGReaderError(s)
        self.atom_names += reactionquery.reactantquery[tree[0][1]].atom_names
        self.electronbalance +=\
            [0]*len(reactionquery.reactantquery[tree[0][1]].atom_names)
        self.atom_belonging_mol +=\
            [0]*len(reactionquery.reactantquery[tree[0][1]].name)

    def ReadDuplicates(self, tree, reactionquery):
        assert tree[0][0] == 'ReactantName'
        assert tree[1][0] == 'ReactantName'
        assert tree[2][0] == 'LabelMapping'
        labelmapping = self.LabelMapping(tree[2][1:])
        if len(labelmapping) != len(reactionquery.reactantquery[tree[1][1]].
                                    atom_names):
            raise RINGReaderError('ReadDuplicates: Labelmapping length',
                                  "doesn''t match length of fragment",
                                  'duplicated')
        reactionquery.reactantquery[tree[0][1]] = reactionquery.\
            reactantquery[tree[1][1]]
        reactionquery.reactantquery[tree[0][1]].name = tree[0][1]
        for i in range(0, len(reactionquery.reactantquery[tree[0][1]].
                              atom_names)):
            try:
                reactionquery.reactantquery[tree[0][1]].atom_names[i] = \
                    labelmapping[reactionquery.reactantquery[tree[0][0]].
                                 atom_names[i]]
            except KeyError:
                s = 'Unrecognized label '
                s += reactionquery.reactantquery[tree[0][1]].atom_names[i]
                s += ' in reactant ' +\
                    reactionquery.reactantquery[tree[1][1]].name
                s += ' for duplication opreation for reactant ' + tree[0][1]
                raise RINGReaderError(s)
        self.atom_names += reactionquery.reactantquery[tree[0][1]].atom_names
        self.electronbalance +=\
            [0]*len(reactionquery.reactantquery[tree[0][1]].atom_names)
        self.atom_belonging_mol += reactionquery.reactantquery[tree[0][1]].\
            name*len(reactionquery.reactantquery[tree[0][1]].atom_names)

    def ReadReactants(self, tree, reactionquery):
        assert tree[0][0] in ['ReactantQuery', 'ReactantGroup', 'Duplicates']
        if tree[0][0] == 'ReactantQuery':
            molquery = MolQueryReader(tree[0][1:]).Read()
            reactionquery.AppendReactantQuery(molquery)
            self.atom_names += molquery.atom_names
            self.electronbalance += [0]*len(molquery.atom_names)
            self.atom_belonging_mol += [molquery.name]*len(molquery.atom_names)
        elif tree[0][0] == 'ReactantGroup':
            self.ReadReactantGroup(tree[1][1:], reactionquery)
        elif tree[0][0] == 'Duplicates':
            self.ReadDuplicates(tree[1][1:], reactionquery)
        if len(tree) == 2:
            self.ReadReactants(tree[1][1:], reactionquery)

    def Read(self):
        # initialize Enhanced ReactionQuery
        reactionquery = ReactionQuery()
        # fragment name
        i = 0
        assert self.tree[i][0] == 'ReactionName'
        reactionquery.name = self.tree[i][1]
        i += 1
        # start building reactionquery
        assert self.tree[i][0] == 'Reactants'
        self.ReadReactants(self.tree[i][1:], reactionquery)
        i += 1
        if self.tree[i][0] == 'Constraints':
            raise NotImplementedError("Constraints not supported yet")
            self.ReadConstraints(self.tree[i][1:], reactionquery)
            i += 1
        # graph transformation is done here.
        assert self.tree[i][0] == 'TransformationChain'
        self.ReadTransformationChain(self.tree[i][1:], reactionquery)
        # error check
        s = str()
        for i in range(0, len(self.electronbalance)):
            if self.electronbalance[i] != 0:
                s += "Electron balance mismatch for atom label: '" +\
                    self.atom_names[i] + "' by ("\
                    + str(self.electronbalance[i]) + ")\n"
        if s:
            raise RINGReaderError(s)
        # return
        return reactionquery
