from .. Error import RINGReaderError
from .. RDkitWrapper.MolQuery import AtomRing,AtomConnectivityAtom,\
    ConstraintNumber,MolCharge,MolAromatic,MolOlefinic,MolParaffinic,\
    MolCyclic,MolLinear,MolQuery, AtomRadical, BondQuery, BondConstraint, \
    AtomConnectivityGroup, AtomIsInRing, AtomIsAromatic, AtomIsAllylic
from rdkit.Chem import rdqueries
from rdkit import Chem
"""
Example:
from msr.RINGParser.Reader import Read
from rdkit import Chem 
s = " \
fragment a{C labeled C1 \
C labeled C2 ring bond to C1} "
mplquery = Read(s)
mol = Chem.MolFromSmiles('C1CCCCC1')
molquery.GetQueryMatches(mol)
"""
class MolQueryReader(object):
    def __init__(self, tree, RINGgroups=None):
        # ast = Abstract Syntax Tree
        self.tree = tree
        # dict of groups where the key is equal to group name
        self.RINGgroups=RINGgroups
        
    def ReadBondTypeAtomConstraint(self,tree):
        return BondQuery(tree[0])
        
    def ReadGroupName(self, tree):
        return tree[0]
    
    def ReadAtomConstraintConnectivity(self, tree):
        i = 0
        # Boolean
        if tree[i][0] == 'Boolean':
            s = tree[i][1]
            if s != '!': 
                raise NotImplementedError("Unsupported Boolean operator: '"+s+"'")
            i += 1
            negate = True
        else:
            negate = False
        # ConstraintNumber
        if tree[i][0] == 'ConstraintNumber':
            CN = ConstraintNumber(tree[i][1:])
            i += 1
        else:
            CN = ConstraintNumber('>=1')
        # Connected atom/group
        if tree[i][0] == 'AtomType':
            t = 0
            connected, constraints = self.ReadAtomType(tree[i][1:])
            # rdkit atom object
            i += 1
        elif tree[i][0] == 'GroupName':
            t = 1
            try:
                connected = self.RINGgroups[self.ReadGroupName(tree[i][1:])]
            except KeyError:
                raise RINGReaderError("Unrecognized group name :'"+self.ReadGroupName(tree[i][1:])+"'")
            i += 1
            # string
        # Bond type
        if len(tree) == i+1:
            assert tree[i][0] == 'BondType'
            bondquery = self.ReadBondTypeAtomConstraint(tree[i][1:])
        else:
            bondquery = BondQuery('single')
        if t == 0:
            return AtomConnectivityAtom(negate,CN,connected,bondquery,constraints)
        elif t == 1:
            return AtomConnectivityGroup(negate,CN,connected,bondquery)
            
    def ReadAtomConstraintRing(self, tree):
        i = 0
        # Boolean
        if tree[i][0] == 'Boolean':
            s = tree[i][1]
            if s != '!': 
                raise NotImplementedError("Unsupported Boolean operator: '"+s+"'")
            i += 1
            negate = True
        else:
            negate = False
        assert tree[i][0] == 'ConstraintNumber'
        CN = ConstraintNumber(tree[i][1:])
        return AtomRing(negate,CN)
        
    def ReadAtomConstraintRadical(self, tree):
        i = 0
        # Boolean
        if tree[i][0] == 'Boolean':
            s = tree[i][1]
            if s != '!': 
                raise NotImplementedError("Unsupported Boolean operator: '"+s+"'")
            i += 1
            negate = True
        else:
            negate = False
        assert tree[i][0] == 'ConstraintNumber'
        CN = ConstraintNumber(tree[i][1:])
        return AtomRadical(negate,CN)
    
    def ReadAtomConstraints(self, tree):
        assert tree[0][0] in ('AtomConstraintConnectivity',\
            'AtomConstraintRing','AtomConstraintRadical')
        if tree[0][0] == 'AtomConstraintConnectivity':
            return self.ReadAtomConstraintConnectivity(tree[0][1:])
        elif tree[0][0] == 'AtomConstraintRing':
            return self.ReadAtomConstraintRing(tree[0][1:])
        elif tree[0][0] == 'AtomConstraintRadical':
            return self.ReadAtomConstraintRadical(tree[0][1:])
    
    def ReadAtomConstraintChain(self, tree, molquery, idx):
        assert tree[0][0] == 'AtomConstraints'
        constraint = self.ReadAtomConstraints(tree[0][1:])
        molquery.AppendAtomConstraint(constraint,idx)
        
        if len(tree) > 1:
            assert tree[1][0] == 'AtomConstraintChain'
            self.ReadAtomConstraintChain(tree[1][1:], molquery, idx)
        
    def ReadAtomSuffix(self, tree, atom):
        constraint = None
        #'+','-','.',':','+.','-.','*'
        if tree[0] == '+.':
            constraint = AtomRadical(False,ConstraintNumber('=1'))
            atom.ExpandQuery(rdqueries.FormalChargeEqualsQueryAtom(1))
        elif tree[0] == '-.':
            constraint = AtomRadical(False,ConstraintNumber('=1'))
            atom.ExpandQuery(rdqueries.FormalChargeEqualsQueryAtom(-1))
        elif tree[0] == '+':
            atom.ExpandQuery(rdqueries.FormalChargeEqualsQueryAtom(1))
        elif tree[0] == '-':
            atom.ExpandQuery(rdqueries.FormalChargeEqualsQueryAtom(-1))
        elif tree[0] == '.':
            constraint = AtomRadical(False,ConstraintNumber('=1'))
        elif tree[0] == ':':
            constraint = AtomRadical(False,ConstraintNumber('=2'))
        elif tree[0] == ':.':
            constraint = AtomRadical(False,ConstraintNumber('=3'))
        elif tree[0] == '*':
            from rdkit.Chem import GetPeriodicTable
            #if type(atom).__name__ == 'QueryAtom':
            #    raise NotImplementedError('Onium $,&,X atoms not supported yet')
            atomicnum = atom.GetAtomicNum()
            atom = rdqueries.AtomNumEqualsQueryAtom(atomicnum)
            valence = GetPeriodicTable().GetDefaultValence(atomicnum)
            atom.ExpandQuery(rdqueries.TotalValenceEqualsQueryAtom(valence+1))
            atom.ExpandQuery(rdqueries.FormalChargeEqualsQueryAtom(1))
        elif tree[0] == '?':
            pass
        else:
            s = "Unsupported atom suffic: '" + tree[0] + "'"
            raise NotImplementedError(s)
        return constraint
            
    def ReadSymbols(self, tree):
        if tree[0] in ['any atom','$']:
            atom = rdqueries.AtomNumGreaterQueryAtom(0)
        elif tree[0] in ['heteroatom','&']:
            #N, O, P, S
            atom = rdqueries.AtomNumEqualsQueryAtom(7)
            atom.ExpandQuery(rdqueries.AtomNumEqualsQueryAtom(8),\
                how=Chem.rdchem.CompositeQueryType.COMPOSITE_OR)
            atom.ExpandQuery(rdqueries.AtomNumEqualsQueryAtom(15),\
                how=Chem.rdchem.CompositeQueryType.COMPOSITE_OR)
            atom.ExpandQuery(rdqueries.AtomNumEqualsQueryAtom(16),\
                how=Chem.rdchem.CompositeQueryType.COMPOSITE_OR)
        elif tree[0] in ['heavy atom','X']:
            # heavier than H
            atom = rdqueries.AtomNumGreaterQueryAtom(1)
        elif tree[0][0].islower():
            # aromatic molecule
            symbol = tree[0][0].upper()+tree[0][1:]
            try:
                atom = Chem.Atom(symbol)
                atom.SetIsAromatic(True)
            except RuntimeError:
                msg = 'Element aromatic ' +symbol+ ' not found'
                raise RINGReaderError(msg)
        elif tree[0] == 'M':
            # metal
            atom = rdqueries.AtomNumGreaterQueryAtom(19)
        else:
            try:
                atom = Chem.Atom(tree[0])
                atom = rdqueries.AtomNumEqualsQueryAtom(atom.GetAtomicNum())
            except RuntimeError:
                msg = 'Element ' +tree[0]+ ' not found'
                raise RINGReaderError(msg)
        return atom
    def ReadAtomPrefix(self, tree):
        if tree[0] == 'aromatic':
            return AtomIsAromatic(negate=False)
        elif tree[0] == 'nonaromatic':
            return AtomIsAromatic(negate=True)
        elif tree[0] == 'ringatom':
            return AtomIsInRing(negate=False)
        elif tree[0] == 'nonringatom':
            return AtomIsInRing(negate=True)
        elif tree[0] == 'allylic':
            return AtomIsAllylic(negate=False)
        
    def ReadAtomType(self, tree):
        i=0
        constraints = list()
        if tree[i][0] == 'AtomPrefix':
            constraints.append(self.ReadAtomPrefix(tree[i][1:]))
            i += 1
        
        assert tree[i][0] == 'Symbols'
        atom = self.ReadSymbols(tree[i][1:])
        i += 1        
        
        if len(tree) > i:
            constraint = self.ReadAtomSuffix(tree[i][1:], atom)
            if constraint:
                constraints.append(constraint)
        else:
            # no suffix means it's neutral. RDkit atom with 0 formal charge and
            # radical electron finds all variables of formal charge and radical
            # electrons, so they need to be treated.
            atom.ExpandQuery(rdqueries.FormalChargeEqualsQueryAtom(0))
            constraints.append(AtomRadical(False,ConstraintNumber('=0')))
        # return atom
        return atom, constraints
    
    def ReadAtom(self, tree, molquery):
        assert tree[0][0] == 'AtomType'
        atom, constraints = self.ReadAtomType(tree[0][1:])
        idx = molquery.mol.AddAtom(atom)
        if constraints:
            for constraint in constraints:
                molquery.AppendAtomConstraint(constraint,idx)
        assert tree[1][0] == 'AtomLabel'
        molquery.atom_names.append(tree[1][1])
        if len(tree) > 2:
            assert tree[2][0] == 'AtomConstraintChain'
            self.ReadAtomConstraintChain(tree[2][1:], molquery, idx)
            
    def ReadBondTypeBondedAtom(self,idx,idx_connected,bondtype,molquery):
        if bondtype == 'single':
            molquery.mol.AddBond(idx,idx_connected,Chem.BondType.SINGLE)
        elif bondtype == 'double':
            molquery.mol.AddBond(idx,idx_connected,Chem.BondType.DOUBLE)
        elif bondtype == 'triple':
            molquery.mol.AddBond(idx,idx_connected,Chem.BondType.TRIPLE)
        elif bondtype == 'quadruple':
            molquery.mol.AddBond(idx,idx_connected,Chem.BondType.QUADRUPLE)
        elif bondtype == 'ring':
            molquery.mol.AddBond(idx,idx_connected,Chem.BondType.UNSPECIFIED)
            molquery.AppendBondConstraint(idx,idx_connected,BondConstraint(BondQuery(bondtype)))
        elif bondtype == 'nonring':
            molquery.mol.AddBond(idx,idx_connected,Chem.BondType.UNSPECIFIED)
            molquery.AppendBondConstraint(idx,idx_connected,BondConstraint(BondQuery(bondtype)))
        elif bondtype == 'aromatic':
            molquery.mol.AddBond(idx,idx_connected,Chem.BondType.AROMATIC)
        elif bondtype == 'any':
            molquery.mol.AddBond(idx,idx_connected,Chem.BondType.UNSPECIFIED)
            return True
        elif bondtype == 'strong':
            molquery.mol.AddBond(idx,idx_connected,Chem.BondType.UNSPECIFIED)
            molquery.AppendBondConstraint(idx,idx_connected,BondConstraint(BondQuery(bondtype)))
        elif bondtype == 'partial':
            molquery.mol.AddBond(idx,idx_connected,Chem.BondType.UNSPECIFIED)
            molquery.AppendBondConstraint(idx,idx_connected,BondConstraint(BondQuery(bondtype)))
        else:
            raise NotImplementedError("Unsupported bond type: '"+bondtype+"'")
            
    def ReadBondedAtom(self, tree, molquery):
        assert tree[0][0] == 'AtomType'
        atom, constraints = self.ReadAtomType(tree[0][1:])
        idx = molquery.mol.AddAtom(atom)
        if constraints:
            for constraint in constraints:
                molquery.AppendAtomConstraint(constraint,idx)
            
        assert tree[1][0] == 'AtomLabel'
        if tree[1][0] in molquery.atom_names:
            raise RINGReaderError('Atom Label '+tree[1][0]+' is alreadyd declared!')
        molquery.atom_names.append(tree[1][1])
        assert tree[2][0] == 'BondType'
        bondtype = tree[2][1:][0]
        assert tree[3][0] == 'AtomLabel'
        try:
            idx_connected = molquery.atom_names.index(tree[3][1])
        except:
            msg = 'Atom Label '+tree[3][1]+' not found'
            raise RINGReaderError(msg)
        self.ReadBondTypeBondedAtom(idx,idx_connected,bondtype,molquery)
        
        if len(tree) > 4:
            assert tree[4][0] == 'AtomConstraintChain'
            self.ReadAtomConstraintChain(tree[4][1:], molquery, idx)
            
    def ReadRingBond(self, tree, molquery):
        assert tree[0][0] == 'AtomLabel'
        try:
            idx1 = molquery.atom_names.index(tree[0][1])
        except:
            msg = 'Atom Label '+tree[0][1]+' not found'
            raise RINGReaderError(msg)
        
        assert tree[1][0] == 'BondType'
        bondtype = tree[1][1:][0]
        assert tree[2][0] == 'AtomLabel'
        try:
            idx2 = molquery.atom_names.index(tree[2][1])
        except:
            msg = 'Atom Label '+tree[2][1]+' not found'
            raise RINGReaderError(msg)
        self.ReadBondTypeBondedAtom(idx1,idx2,bondtype,molquery)
            
    def ReadAtomChain(self, tree, molquery):
        assert tree[0][0] in ['BondedAtom','RingBond']
        if tree[0][0] == 'BondedAtom':
            self.ReadBondedAtom(tree[0][1:], molquery)
        elif tree[0][0] == 'RingBond':
            self.ReadRingBond(tree[0][1:], molquery)
        
        if len(tree) > 1:
            assert tree[1][0] == 'AtomChain'
            self.ReadAtomChain(tree[1][1:], molquery)
        
        
    def ReadMolQuery(self, tree, molquery):
        assert tree[0][0] == 'Atom'
        self.ReadAtom(tree[0][1:], molquery)
        
        if len(tree) > 1:
            assert tree[1][0] == 'AtomChain'
            self.ReadAtomChain(tree[1][1:], molquery)
        
        
    def ReadMolQueryPrefix(self,tree,molquery):
        i = 0
        if tree[i] == 'positive':
            molquery.AppendMolConstraint(MolCharge(ConstraintNumber('=1')))
            i+=1
        elif tree[i] == 'negative':
            molquery.AppendMolConstraint(MolCharge(ConstraintNumber('=-1')))
            i+=1
        elif tree[i] == 'neutral':
            molquery.AppendMolConstraint(MolCharge(ConstraintNumber('=0')))
            i+=1
            
        if i<len(tree):
            if tree[i] == 'aromatic':
                molquery.AppendMolConstraint(MolAromatic())
                i+=1
            elif tree[i] == 'olefinic':
                molquery.AppendMolConstraint(MolOlefinic())
                i+=1
            elif tree[i] == 'paraffinic':
                molquery.AppendMolConstraint(MolParaffinic())
                i+=1
                
        if i<len(tree):
            if tree[i] == 'cyclic':
                molquery.AppendMolConstraint(MolCyclic())
            elif tree[i] == 'linear':
                molquery.AppendMolConstraint(MolLinear())
            else:
                raise NotImplementedError('Unsupported mol prefix: %s'%tree[0])
        
    def Read(self):
        # initialize Enhanced MolQuery
        molquery = MolQuery()
        # prefix
        assert self.tree[0][0] == 'Prefix'
        if len(self.tree[0]) != 1:
            self.ReadMolQueryPrefix(self.tree[0][1:],molquery)
        # fragment name
        assert self.tree[1][0] in ['FragmentName','ReactantName','GroupName']
        molquery.name = self.tree[1][1]
        # start building molecule
        assert self.tree[2][0] == 'MolQuery'
        self.ReadMolQuery(self.tree[2][1:], molquery)
        # return
        return molquery

