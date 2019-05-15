# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 11:22:49 2018

@author: gerha
"""
from pGrAdd.RINGParser.Reader import Read
from rdkit import Chem


## fragment prefix test
testmol = Chem.MolFromSmiles('C[CH2+]')
s = """
positive fragment a{
    C labeled c1
    C labeled c2 single bond to c1
}
"""
molquery = Read(s)


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
    
# Break Bond, Form bond, Increase Radical electrons
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

# Break Bond, Form bond, Increase Radical electrons
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

# Break Bond, Form bond, modify Radical electrons
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