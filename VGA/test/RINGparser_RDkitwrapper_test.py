# -*- coding: utf-8 -*-
from VGA.RINGParser.Reader import Read
from rdkit import Chem
# molquery test
## Basic connectivity test
testmol = Chem.MolFromSmiles('CCC')
s = """
fragment a{
    C labeled c1
    C labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1), (1, 0), (1, 2), (2, 1))
## Double bond, triple bond test
testmol = Chem.MolFromSmiles('C=C-C#C')
s = """
fragment a{
    C labeled c1
    C labeled c2 double bond to c1
    C labeled c3 single bond to c2
    C labeled c4 triple bond to c3
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0,1,2,3),)
## aromatic bond test
testmol = Chem.MolFromSmiles('c1ccccc1')
s = """
fragment a{
    C labeled c1
    C labeled c2 aromatic bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1), (0, 5), (1, 0), (1, 2), (2, 1), (2, 3), (3, 2), (3, 4), (4, 3), (4, 5), (5, 0), (5, 4))
## ring bond test
testmol = Chem.MolFromSmiles('CCC1CCC1')
s = """
fragment a{
    C labeled c1
    C labeled c2 ring bond to c1
    C labeled c3 ring bond to c2
    C labeled c4 ring bond to c3
    ringbond c4 ring bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((2, 3, 4, 5), (2, 5, 4, 3), (3, 2, 5, 4), (3, 4, 5, 2), (4, 3, 2, 5), (4, 5, 2, 3), (5, 2, 3, 4), (5, 4, 3, 2))
## ring bond test
testmol = Chem.MolFromSmiles('CCC1CCC1')
s = """
fragment a{
    C labeled c1
    C labeled c2 ring bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((2, 3), (2, 5), (3, 2), (3, 4), (4, 3), (4, 5), (5, 2), (5, 4))
## non-ring bond test
testmol = Chem.MolFromSmiles('CCC1CCC1')
s = """
fragment a{
    C labeled c1
    C labeled c2 nonring bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1), (1, 0), (1, 2), (2, 1))
## any bond test
testmol = Chem.MolFromSmiles('CC=CC#C')
s = """
fragment a{
    C labeled c1
    C labeled c2 any bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1), (1, 0), (1, 2), (2, 1), (2, 3), (3, 2), (3, 4), (4, 3))
## strong bond test
testmol = Chem.MolFromSmiles('CC=CC#C')
s = """
fragment a{
    C labeled c1
    C labeled c2 strong bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((1, 2), (2, 1), (3, 4), (4, 3))
## fragment prefix test
testmol = Chem.MolFromSmiles('C[CH2+]')
s = """
positive fragment a{
    C labeled c1
    C labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1), (1, 0))
testmol = Chem.MolFromSmiles('C[CH2-]')
s = """
negative fragment a{
    C labeled c1
    C labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1), (1, 0))
testmol = Chem.MolFromSmiles('C[CH2-]')
s = """
positive fragment a{
    C labeled c1
    C labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ()
testmol = Chem.MolFromSmiles('C[CH2+]')
s = """
negative fragment a{
    C labeled c1
    C labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ()
testmol = Chem.MolFromSmiles('C=CC')
s = """
olefinic fragment a{
    C labeled c1
    C labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((1,2),(2,1))
testmol = Chem.MolFromSmiles('C=CC')
s = """
paraffinic fragment a{
    C labeled c1
    C labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ()
testmol = Chem.MolFromSmiles('CCC')
s = """
paraffinic fragment a{
    C labeled c1
    C labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1), (1, 0), (1, 2), (2, 1))
testmol = Chem.MolFromSmiles('CCC')
s = """
linear fragment a{
    C labeled c1
    C labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1), (1, 0), (1, 2), (2, 1))
testmol = Chem.MolFromSmiles('CCC')
s = """
cyclic fragment a{
    C labeled c1
    C labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ()
testmol = Chem.MolFromSmiles('C1CCC1C')
s = """
cyclic fragment a{
    C labeled c1
    C labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1), (0, 3), (1, 0), (1, 2), (2, 1), (2, 3), (3, 0),\
    (3, 2), (3, 4), (4, 3))
# symbols and atom suffix test
testmol = Chem.MolFromSmiles('CCC')
s = """
fragment a{
    $ labeled c1
    $ labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1), (0, 3), (0, 4), (0, 5), (1, 0), (1, 2), (1, 6),\
    (1, 7), (2, 1), (2, 8), (2, 9), (2, 10), (3, 0), (4, 0), (5, 0), (6, 1),\
    (7, 1), (8, 2), (9, 2), (10, 2))
testmol = Chem.MolFromSmiles('CCO')
s = """
fragment a{
    X labeled c1
    X labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1), (1, 0), (1, 2), (2, 1))
testmol = Chem.MolFromSmiles('CCS')
s = """
fragment a{
    X labeled c1
    & labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((1, 2),)
testmol = Chem.MolFromSmiles('[CH][CH][C+][C-][C][CH+][CH-][CH4+]')
s = """
fragment a{
    X. labeled c1
    X: labeled c2 single bond to c1
    C+. labeled c3 single bond to c2
    C-. labeled c4 single bond to c3
    C: labeled c5 single bond to c4
    C+ labeled c6 single bond to c5
    C- labeled c7 single bond to c6
    C* labeled c8 single bond to c7
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1, 2, 3, 4, 5, 6, 7),)
# atom constraint test
testmol = Chem.MolFromSmiles('CCC')
s = """
fragment a{
    C labeled c1
    C labeled c2 single bond to c1 {connected to =2 C}
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
#print match_index
assert match_index == ((0, 1), (2, 1))
testmol = Chem.MolFromSmiles('CCC')
s = """
fragment a{
    C labeled c1
    C labeled c2 single bond to c1 {connected to =1 C}
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((1, 0), (1, 2))
testmol = Chem.MolFromSmiles('CC=C')
s = """
fragment a{
    C labeled c1
    C labeled c2 single bond to c1 {connected to =1 C with double bond}
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1),)
testmol = Chem.MolFromSmiles('CC=C')
s = """
fragment a{
    C labeled c1 {connected to >1 C with any bond}
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((1,),)
testmol = Chem.MolFromSmiles('CC=C')
s = """
fragment a{
    C labeled c1 {!connected to >1 C with any bond}
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0,),(2,))
testmol = Chem.MolFromSmiles('CC1CCC1')
s = """
fragment a{
    C labeled c1 {!in ring of size >0}
    C labeled c2 single bond to c1 {in ring of size >0}
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1),)
# atom prefix
testmol = Chem.MolFromSmiles('CC1CCC1')
s = """
fragment a{
    nonringatom C labeled c1
    ringatom C labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1),)
testmol = Chem.MolFromSmiles('Cc1ccccc1')
s = """
fragment a{
    nonaromatic C labeled c1
    aromatic C labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1),)
testmol = Chem.MolFromSmiles('CC=C')
s = """
fragment a{
    C labeled c1
    allylic C labeled c2 single bond to c1
}
"""
molquery = Read(s)
match_index = molquery.GetQueryMatches(testmol)
assert match_index == ((0, 1),)
# Reaction
## Break Bond, Form bond, Decrease bond order
reactants = list()
reactants.append(Chem.MolFromSmiles('C#C'))
reactants.append(Chem.MolFromSmiles('[H][H]'))
s = """
rule test{
    reactant a{
        C labeled c1
        C labeled c2 strong bond to c1}
    reactant b{
        H labeled h1
        H labeled h2 single bond to h1}
    break single bond (h1,h2)
    form single bond (c1,h1)
    form single bond (c2,h2)
    decrease bond order (c1,c2)
}
"""
reaction_query = Read(s)
pl = reaction_query.RunReactants(reactants)
for products in pl:
    assert Chem.MolToSmiles(products[0]) == '[H]C([H])=C([H])[H]'
## Break Bond, Form bond, Modify Bond
reactants = list()
reactants.append(Chem.MolFromSmiles('C#C'))
reactants.append(Chem.MolFromSmiles('[H][H]'))
s = """
rule test{
    reactant a{
        C labeled c1
        C labeled c2 triple bond to c1}
    reactant b{
        H labeled h1
        H labeled h2 single bond to h1}
    break single bond (h1,h2)
    form single bond (c1,h1)
    form single bond (c2,h2)
    modify bond (c1,c2,double)
}
"""
reaction_query = Read(s)
pl = reaction_query.RunReactants(reactants)
for products in pl:
    assert Chem.MolToSmiles(products[0]) == '[H]C([H])=C([H])[H]'
## Break Bond, Form bond, Increase bond order
reactants = list()
reactants.append(Chem.MolFromSmiles('CC'))
s = """
rule test{
    reactant a{
        C labeled c1
        C labeled c2 single bond to c1
        H labeled h1 single bond to c1
        H labeled h2 single bond to c2}
    break single bond (h1,c1)
    break single bond (h2,c2)
    form single bond (h1,h2)
    increase bond order (c1,c2)
}
"""
reaction_query = Read(s)
pl = reaction_query.RunReactants(reactants)
for products in pl:
    assert Chem.MolToSmiles(products[0]) == '[H]C([H])=C([H])[H]'
    assert Chem.MolToSmiles(products[1]) == '[H][H]'
## Break Bond, Form bond, Atom Modify
reactants = list()
reactants.append(Chem.MolFromSmiles('C#C'))
s = """
rule test{
    reactant a{
        C labeled c1
        C labeled c2 triple bond to c1}
    decrease bond order (c1,c2)
    modify atomtype(c1,C.)
    modify atomtype(c2,C.)
}
"""
reaction_query = Read(s)
pl = reaction_query.RunReactants(reactants)
for products in pl:
    assert Chem.MolToSmiles(products[0]) == '[H][C]=[C][H]'
## Break Bond, Form bond, Increase Radical electrons
reactants = list()
reactants.append(Chem.MolFromSmiles('C#C'))
s = """
rule test{
    reactant a{
        C labeled c1
        C labeled c2 triple bond to c1}
    decrease bond order (c1,c2)
    increase number of radical(c1)
    increase number of radical(c2)
}
"""
reaction_query = Read(s)
pl = reaction_query.RunReactants(reactants)
for products in pl:
    assert Chem.MolToSmiles(products[0]) == '[H][C]=[C][H]'
## Break Bond, Form bond, Increase Radical electrons
reactants = list()
reactants.append(Chem.MolFromSmiles('[CH][CH]'))
s = """
rule test{
    reactant a{
        C: labeled c1
        C: labeled c2 single bond to c1}
    increase bond order (c1,c2)
    decrease number of radical(c1)
    decrease number of radical(c2)
}
"""
reaction_query = Read(s)
pl = reaction_query.RunReactants(reactants)
for products in pl:
    assert Chem.MolToSmiles(products[0]) == '[H][C]=[C][H]'
## Break Bond, Form bond, modify Radical electrons
reactants = list()
reactants.append(Chem.MolFromSmiles('[CH][CH]'))
s = """
rule test{
    reactant a{
        C: labeled c1
        C: labeled c2 single bond to c1}
    increase bond order (c1,c2)
    modify number of radical(c1,1)
    modify number of radical(c2,1)
}
"""
reaction_query = Read(s)
pl = reaction_query.RunReactants(reactants)
for products in pl:
    assert Chem.MolToSmiles(products[0]) == '[H][C]=[C][H]'