from rdkit import Chem
from ..RDkitWrapper.GenRxnNet import GenerateRxnNet
from ..RINGParser.Reader import Read


def enumerate_res_struc(mol):
    """
    Generates resonance structures. Structure with minimum number of radical
    electrons is outputted

    Arguments:
    - mol:     Input. It can be smiles or mol

    Return:
    - list of resonance structures

    Example:
    from rdkit import Chem
    from msr.GroupAdd.FindResStruc import enumerate_res_struc
    A = enumerate_res_struc('[CH][C][C][O]')
    for product in A:
        NRE = 0
        for atoms in product.GetAtoms():
            NRE += atoms.GetNumRadicalElectrons()
        print(NRE, Chem.MolToSmiles(product))

    """
    # Initialize
    # Does not accept list
    if isinstance(mol, list):
        raise ValueError('Only single molecule accepted')
    # Convert to rdkit mol
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)

    # Delete all double/triple bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
            bond.SetBondType(Chem.rdchem.BondType.SINGLE)
            atom = bond.GetBeginAtom()
            atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons()+2)
            atom = bond.GetEndAtom()
            atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons()+2)
        elif bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            bond.SetBondType(Chem.rdchem.BondType.SINGLE)
            atom = bond.GetBeginAtom()
            atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons()+1)
            atom = bond.GetEndAtom()
            atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons()+1)
        elif bond.GetBondType() == Chem.rdchem.BondType.QUADRUPLE:
            bond.SetBondType(Chem.rdchem.BondType.SINGLE)
            atom = bond.GetBeginAtom()
            atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons()+3)
            atom = bond.GetEndAtom()
            atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons()+3)

    # Chem.AssignRadicals(mol) # This is broken. I don't know how it works.
    # Frustrating
    print('Resonance Structure')
    # Get all possible resonance structures
    # set up reactions
    reaction_rules = list()
    # Conjugated bond formation
    temp = """
    rule increaseBO{
        reactant r1{
        $? labeled a1 {has >0 radical electrons}
        $? labeled a2 any bond to a1 {has >0 radical electrons}
        }
        increase bond order (a1,a2)
        decrease number of radical (a1)
        decrease number of radical (a2)
    }
    """
    reaction_rules.append(Read(temp))
    temp = """
    rule increaseBO{
        reactant r1{
        $+ labeled a1
        $ labeled a2 any bond to a1 {has >0 radical electrons}
        }
        increase number of radical (a1)
        decrease number of radical (a2)
        decrease formal charge (a1)
        increase formal charge (a2)
    }
    """
    reaction_rules.append(Read(temp))

    # run
    resonance_structures = GenerateRxnNet(mol, reaction_rules)
    """
    Structure must:
    obey as much as possible the octet rule (8 valence electrons around each
    atom rather than having deficiencies or surplus) have a maximum number
    of covalent bonds carry a minimum of charged atoms. If unlike charges
    are present their separation must be least while for like charges the
    separation must be maximum.
    """
    # Octet rule
    # here this is achieved by finding the minimum number of radical electron.
    # Find the minimum number of radicals
    TotalNumRadicalElectron = list()
    for resonance_structure in resonance_structures:
        NTRE = 0  # Number of Total Radical Electorn
        for atoms in resonance_structure.GetAtoms():
            NTRE += atoms.GetNumRadicalElectrons()
        TotalNumRadicalElectron.append(NTRE)
    MinTotalNumRadicalElectron = min(TotalNumRadicalElectron)
    # Remove resonance structure with number of radical
    # electron more than the min
    for i in range(len(TotalNumRadicalElectron)-1, -1, -1):
        if TotalNumRadicalElectron[i] > MinTotalNumRadicalElectron:
            del resonance_structures[i]

    return resonance_structures
