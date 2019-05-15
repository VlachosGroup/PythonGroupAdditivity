""" This is a RDKIT wrapper code that generate reaction network"""

from rdkit import Chem
from rdkit.Chem.AllChem import ReactionFromSmarts
from rdkit.Chem.rdchem import GetPeriodicTable
from rdkit.Chem.rdchem import PeriodicTable
from itertools import product as itpd
from .. RINGParser import Read


def GenerateRxnNet(initial_reactant, reaction_rules):
    """
    Generates reaction network following the algorithm from
    (Ind. Eng. Chem. Res. 2010, 49 (21), 10459-10470)

    Arguments:
    - initial_reactant:     can be smiles or mol
    - reaction_rules:       can be smarts string or rdkit.Chem.rdChemReactions.
                            ChemicalReaction object

    Return:
    - list of reaction intermediates

    Example:
    Ethane C-H scission and C-C scission
    A = generate_rxn_net('CC',['[C:1][H:2]>>[C:1].[H:2]',
                         '[C:1][C:2]>>[C:1].[C:2]'])
    for product in A:
        print(Chem.MolToSmiles(product))

    Tip:
    For dehydrogenation, you must number carbon and hydrogen i.e.
    [C][H]>>[C].[H]             (x)
    [C:1][H:1]>>[C:1].[H:2]     (o)
    Same goes for changing bond order
    [C][C]>>[C]=[C]             (x)
    [C:1][C:2]>>[C:1]=[C:2]     (o)
    Aromatization and Kekulize has some problem with several species. So, we
    don't set aromatization for sanitize, and try kekulize.

    TODO:
    - record reactions as well (make it an option as it's expensive)

    """

    # set-up reactants
    if not isinstance(initial_reactant, list):
        initial_reactant = [initial_reactant]
    if isinstance(initial_reactant[0], str):
        for i in range(0, len(initial_reactant)):
            initial_reactant[i] = Chem.MolFromSmiles(initial_reactant[i],
                                                     sanitize=False)
            # sanitize everything except aromatization set
            _sanitize_except_aromatization(initial_reactant[i])
    # Treatment necessary for radicals
    # (https://github.com/rdkit/rdkit/issues/69)
    for i in range(0, len(initial_reactant)):
        # print Chem.MolToSmiles(initial_reactant[i])
        initial_reactant[i] = Chem.AddHs(initial_reactant[i])
        # sanitize everything except aromatization set
        _sanitize_except_aromatization(initial_reactant[i])
        # Chem.SanitizeMol(initial_reactant[i])
        # Chem.Kekulize(initial_reactant[i])
        for atoms in initial_reactant[i].GetAtoms():
            atoms.SetNoImplicit(True)
        Chem.AssignRadicals(initial_reactant[i])

    # set up reactions
    if not isinstance(reaction_rules, list):
        reaction_rules = [reaction_rules]
    if isinstance(reaction_rules[0], str):
        for i in range(0, len(reaction_rules)):
            try:
                reaction_rules[i] = Read(reaction_rules[i])
            except Exception:
                reaction_rules[i] = ReactionFromSmarts(reaction_rules[i])

    # generator main algorithm
    unprocessed = initial_reactant
    processed = []
    while unprocessed:
        # Pop a molecule and put it in a processed list
        reactant0 = unprocessed[0]
        processed.insert(0, unprocessed[0])
        del unprocessed[0]
        # go through all reactions
        for reaction_rule in reaction_rules:
            # set up reactant list.
            # Generate combinatorial product list of reactants if reaction
            # requires several reactants
            reactant_list = itpd([list(range(1, len(processed)))],
                                 repeat=reaction_rule.
                                 GetNumReactantTemplates()-1)
            # go through each set of reactants
            for reactant_indexes in reactant_list:
                # Reaction
                # Make the reactant mol tuple (Runreactants only accept tuple)
                reactants = (reactant0,)
                for reactant_index in reactant_indexes:
                    reactants += (processed[reactant_index],)
                # React
                ele_reactions = reaction_rule.RunReactants(reactants)

                # Record reactions (TODO)

                # Pre-processing products
                # Go through reactiosn and make a single list of products
                products = []
                for ele_reaction in ele_reactions:
                    for mol in ele_reaction:
                        products.append(mol)
                # Treatment necessary for radicals
                # (https://github.com/rdkit/rdkit/issues/69)
                for mol in products:
                    for atoms in mol.GetAtoms():
                        atoms.SetNoImplicit(True)
                        atoms.UpdatePropertyCache(strict=False)
                    Chem.AssignRadicals(mol)
                    # Remove molecule with atoms with over valence
                for i in range(len(products)-1, -1, -1):
                    for atoms in products[i].GetAtoms():
                        if PeriodicTable.GetDefaultValence(GetPeriodicTable(),
                                                           atoms.GetAtomicNum()
                                                           ) < \
                             atoms.GetTotalValence():
                            del products[i]
                            break
                # remove duplicates
                # TODO. This removes also species with different charges
                for i in range(len(products)-1, -1, -1):
                    for j in range(0, i):
                        if products[i].GetNumAtoms() ==\
                            products[j].GetNumAtoms() and \
                            products[i].GetNumAtoms() ==\
                                len(products[i].
                                    GetSubstructMatch(products[j])):

                            del products[i]
                            break
                # update unprocessed molecule list
                # check for duplicate and append to unprocessed_list if missing
                for mol1 in products:
                    inthelist = 0
                    for mol2 in processed:
                        # first check the nubmer of atoms and then
                        # look for substructure match
                        if mol1.GetNumAtoms() == mol2.GetNumAtoms() and \
                            mol1.GetNumAtoms() == len(mol1.GetSubstructMatch
                                                      (mol2)):
                            # if it's in processed list, break
                            inthelist = 1
                            break
                    # not in the processed list. append to unprocessed
                    if inthelist == 0:
                        unprocessed.insert(0, mol1)
    # Prettify
    for i in range(0, len(processed)):
        # print Chem.MolToSmiles(processed[i])
        processed[i] = Chem.RemoveHs(processed[i], sanitize=False)
        _sanitize_except_aromatization(processed[i])
        # print Chem.MolToSmiles(processed[i])
    return processed


def _sanitize_except_aromatization(mol):
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.
                         SANITIZE_ADJUSTHS)
        Chem.SanitizeMol(mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.
                         SANITIZE_CLEANUP)
        Chem.SanitizeMol(mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.
                         SANITIZE_CLEANUPCHIRALITY)
        Chem.SanitizeMol(mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.
                         SANITIZE_FINDRADICALS)
        try:
            Chem.SanitizeMol(mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.
                             SANITIZE_KEKULIZE)
        except Exception:
            pass
        Chem.SanitizeMol(mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.
                         SANITIZE_PROPERTIES)
        Chem.SanitizeMol(mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.
                         SANITIZE_SETCONJUGATION)
        Chem.SanitizeMol(mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.
                         SANITIZE_SETHYBRIDIZATION)
        Chem.SanitizeMol(mol, sanitizeOps=Chem.rdmolops.SanitizeFlags.
                         SANITIZE_SYMMRINGS)
