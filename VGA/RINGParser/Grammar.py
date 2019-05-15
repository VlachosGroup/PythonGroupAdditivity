#  coding: utf-8
from .Parser import All                 # This parse class will go through them
#                                         until it runs out.
from .Parser import Either              # This parse class will require at
#                                         least one attribute inside to occur.
#                                         (other wise error)
from .Parser import Optional            # attributes of this class does
#                                         is optional
from .Parser import Literals            # literal string comaprison (several
#                                         string) with sort. So that Cl is
#                                         checked for C.
from .Parser import Digit               # 1 digit
from .Parser import Number              # any length number
from .Parser import String              # any length string
from .Parser import Filler              # Fillers are for checking if correct
#                                         syntanx are given. Are not recorded.
"""
This file contains Abstract Syntax Tree of the RING input. Parser grabs the
tree info from here, and interpret string.
strict_grammar follows the original RING paper abstract tree

Reference:

Rangarajan, S., Kaminski, T., Van Wyk, E., Bhan A., and Daoutidis, P.,
"Language-oriented rule-based reaction network generation and analysis:
    Algorithms of RING", Computers and Chemical Engineering 64 (2014) 124,
    10.1016/j.compchemeng.2014.02.007
Rangarajan, S., Bhan, A., and Daoutidis, P.,
"Language-oriented rule-based reaction network generation and analysis:
    Descrpition of RING", Computers and Chemical Engineering 45 (2012) 114,
    10.1016/j.compchemeng.2012.06.008
"""
strict_grammar = ('RINGInput', {
    'RINGInput': Either('Fragment', 'ReactionRule'),
    'Fragment': All('Prefix', Filler('fragment'), 'FragmentName',
                    Filler('{'), 'MolQuery', Filler('}')),
    'FragmentName': String(),
    'MolQuery': All('Atom', Optional('AtomChain')),
    'Prefix': All(Optional(Literals(['positive', 'negative', 'neutral'])),
                  Optional(Literals(['aromatic', 'olefinic', 'paraffinic'])),
                  Optional(Literals(['cyclic', 'linear']))),
    'Atom': All('AtomType', Filler('labeled'), 'AtomLabel',
                Optional(All(Filler('{'), 'AtomConstraintChain',
                             Filler('}')))),
    'AtomType': All(Optional('AtomPrefix'), 'Symbols', Optional('AtomSuffix')),
    'AtomPrefix': Literals(['aromatic', 'nonaromatic', 'ringatom',
                            'nonringatom', 'allylic']),
    'Symbols': Either(Literals(['any atom', '$', 'heteroatom',
                                '&', 'heavy atom', 'X']), String()),
    'AtomSuffix': Literals(['+', '-', '.', ':', '+.', '-.', '*']),
    'AtomLabel': String(),
    'AtomConstraintChain': All('AtomConstraints',
                               Optional(All(Filler(','),
                                            'AtomConstraintChain'))),
    'AtomConstraints': Either('AtomConstraintConnectivity',
                              'AtomConstraintRing'),
    'AtomConstraintConnectivity': All(Optional('Boolean'),
                                      Filler('connected to'),
                                      Optional('ConstraintNumber'),
                                      Either('GroupName', 'AtomType'),
                                      Optional(All(Filler('with'),
                                                   'BondType',
                                                   Filler('bond')))),
    'AtomConstraintRing': All(Optional('Boolean'),
                              Filler('in ring of size'),
                              'ConstraintNumber'),
    'Boolean': Literals(['!', '||', '&&', '+', '-']),
    'ConstraintNumber': All(Optional(Literals(['>', '=', '<', '>=', '<='])),
                            Digit()),
    'GroupName': All(Filler('group'), String()),
    'BondType': Literals(['single', 'double', 'triple', 'ring', 'nonring',
                          'aromatic', 'any', 'strong', 'partial']),
    'AtomChain': All(Either('BondedAtom', 'RingBond'),
                     Optional('AtomChain')),
    'BondedAtom': All('AtomType', Filler('labeled'), 'AtomLabel', 'BondType',
                      Filler('bond to'), 'AtomLabel',
                      Optional(All(Filler('{'), 'AtomConstraintChain',
                                   Filler('}')))),
    'RingBond': All(Filler('ringbond'), 'AtomLabel', 'BondType',
                    Filler('bond to'), 'AtomLabel'),
    'ReactionRule': All(Filler('rule'), 'ReactionName',
                        Filler('{'), 'Reactants',
                        Optional('Constraints'),
                        'TransformationChain', Filler('}')),
    'ReactionName': String(),
    'Reactants': All(Either('ReactantQuery', 'ReactantGroup', 'Duplicates'),
                     Optional('Reactants')),
    'ReactantQuery': All('Prefix', Filler('reactant'), 'ReactantName',
                         Filler('{'), 'MolQuery', Filler('}')),
    'ReactantName': String(),
    'ReactantGroup': All(Filler('reactant'), 'ReactantName', 'GroupName',
                         Filler('('), 'LabelMapping', Filler(')')),
    'Duplicates': All(Filler('reactant'), 'ReactantName',
                      Filler('duplicates'), 'ReactantName',
                      Filler('('), 'LabelMapping', Filler(')')),
    'LabelMapping': All('AtomLabel', Filler('=>'), 'AtomLabel',
                        Optional(All(Filler(','), 'LabelMapping'))),
    'Constraints': All(Filler('constraints{'), Optional('FragmentChain'),
                       'ConstraintChain', Filler('}')),
    'FragmentChain': All('Fragment', Optional('FragmentChain')),
    'ConstraintChain': All(Optional('Boolean'),
                           Either('BranchConstraint', 'Constraint'),
                           Optional(All('Boolean', 'ConstraintChain'))),
    'BranchConstraint': All(Filler('('),
                            Either('ConstraintChain', 'Constraint'),
                            Filler(')')),
    'Constraint': Either('C_Size', 'C_Charge', 'C_Cyclic', 'C_Characteristic',
                         'C_Fragment', 'C_Group'),
    'C_Size': All('SizeChain', 'ConstraintNumber'),
    'SizeChain': All('Size', Optional(All('Boolean', 'SizeChain'))),
    'Size': All('ReactantName', Filler('.size')),
    'C_Charge': All('ChargeChain', 'ConstraintNumber'),
    'ChargeChain': All('Charge', Optional(All('Boolean', 'ChargeChain'))),
    'Charge': All('ReactantName', Filler('.charge')),
    'C_Cylic': All('ReactantName', Filler('is cyclic')),
    'C_Characteristic': Either('C_Aromatic', 'C_Oxygenate', 'C_Heteroaromatic',
                               'C_Bridged', 'C_DeclaredCharacteristic',
                               'C_Smiles', 'C_Formula'),
    'C_Aromatic': All('ReactantName', Filler('is aromatic')),
    'C_Oxygenate': All('ReactantName', Filler('is oxygenate')),
    'C_Heteroaromatic': All('ReactantName', Filler('is heteroaromatic')),
    'C_Bridged': All('ReactantName', Filler('is bridged')),
    'C_Declaredcharacteristic': All('ReactantName', Filler('is'),
                                    'CharacteristicName'),
    'CharacteristicName': String(),
    'C_Smiles': All('ReactantName', Filler('is'), 'Smiles'),
    'Smiles': String(),
    'C_Formula': All('ReactantName', Filler('.formula is'),
                     'MolecularFormulaChain'),
    'MolecularFormulaChain': All('ElementSymbol', Optional(Number()),
                                 Optional('MolecularFormulaChain')),
    'C_Fragment': All('ReactantName', Filler('contains'),
                      Optional(All('ConstraintNumber', Filler('of'))),
                      'FragmentName'),
    'C_Group': All('ReactantName', Filler('contains'),
                   Optional(All('ConstraintNumber', Filler('of'))),
                   Optional(Filler('group')), 'GroupName'),
    'TransformationChain': All('ConnectivityChange',
                               Optional('TransformationChain')),
    'ConnectivityChange': Either('BondForm', 'BondBreak', 'BondModify',
                                 'BondDecrease', 'BondIncrease',
                                 'AtomTypeModify'),
    'BondForm': All(Filler('form'), Optional('BondType'), Filler('bond'),
                    Filler('('), 'AtomLabel', Filler(','), 'AtomLabel',
                    Filler(')')),
    'BondBreak': All(Filler('break'), Optional('BondType'), Filler('bond'),
                     Filler('('), 'AtomLabel', Filler(','), 'AtomLabel',
                     Filler(')')),
    'BondModify': All(Filler('modify bond'), Filler('('), 'AtomLabel',
                      Filler(','), 'AtomLabel', Filler(','),
                      'BondType', Filler(')')),
    'BondIncrease': All(Filler('increase bond order'), Filler('('),
                        'AtomLabel', Filler(','), 'AtomLabel',
                        Filler(')')),
    'BondDecrease': All(Filler('decrease bond order'), Filler('('),
                        'AtomLabel', Filler(','), 'AtomLabel',
                        Filler(')')),
    'AtomTypeModify': All(Filler('modify atomtype'), Filler('('),
                          'AtomLabel', Filler(','),
                          'AtomType', Filler(')'))
})

# Enhanced for radical electron related functions
enhanced_grammar = ('RINGInput', strict_grammar[1].copy())
enhanced_grammar[1]['BondType'] = Literals(['single', 'double', 'triple',
                                            'quadruple', 'ring', 'nonring',
                                            'aromatic', 'any', 'strong',
                                            'partial'])
enhanced_grammar[1]['AtomConstraints'] =\
    Either('AtomConstraintConnectivity', 'AtomConstraintRing',
           'AtomConstraintRadical', 'AtomConstraintNRing')
enhanced_grammar[1]['AtomConstraintRadical'] =\
    All(Optional('Boolean'), Filler('has'), 'ConstraintNumber',
        Filler('radical electrons'))
enhanced_grammar[1]['ConnectivityChange'] =\
    Either('BondForm', 'BondBreak', 'BondModify', 'BondDecrease',
           'BondIncrease', 'AtomTypeModify', 'RadicalModify',
           'RadicalIncrease', 'RadicalDecrease', 'ChargeIncrease',
           'ChargeDecrease')
enhanced_grammar[1]['RadicalModify'] =\
    All(Filler('modify number of radical'), Filler('('), 'AtomLabel',
        Filler(','), Number(), Filler(')'))
enhanced_grammar[1]['RadicalIncrease'] =\
    All(Filler('increase number of radical'), Filler('('),
        'AtomLabel', Filler(')'))
enhanced_grammar[1]['RadicalDecrease'] =\
    All(Filler('decrease number of radical'), Filler('('),
        'AtomLabel', Filler(')'))
enhanced_grammar[1]['ChargeIncrease'] =\
    All(Filler('increase formal charge'), Filler('('),
        'AtomLabel', Filler(')'))
enhanced_grammar[1]['ChargeDecrease'] =\
    All(Filler('decrease formal charge'),
        Filler('('), 'AtomLabel', Filler(')'))
# For ?, any type of charge, and radicals are accepted.
enhanced_grammar[1]['AtomSuffix'] = Literals(['+', '-', '.', ':', '+.',
                                              '-.', '*', '?', ':.'])
# Stereochemsitry
enhanced_grammar[1]['AtomChain'] =\
    All(Either('BondedAtom', 'RingBond', 'StereoDoubleBond'),
        Optional('AtomChain'))
enhanced_grammar[1]['StereoDoubleBond'] =\
    All(Filler('stereo double bond'), 'AtomLabel',
        Optional('Boolean'), 'DoubleBondStereoType',
        Filler('to'), 'AtomLabel', Filler('for double bond between'),
        'AtomLabel', Filler('and'), 'AtomLabel')
enhanced_grammar[1]['DoubleBondStereoType'] =\
    Literals(['cis', 'trans', 'notspecified'])
# Number of ring constraint
enhanced_grammar[1]['AtomConstraintNRing'] =\
    All(Optional('Boolean'), Filler('in'), 'ConstraintNumber', Filler('ring'))


def update_names(rules):
    for name in rules:
        rules[name].set_name(name)


update_names(strict_grammar[1])
update_names(enhanced_grammar[1])
