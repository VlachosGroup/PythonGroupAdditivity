# -*- coding: utf-8 -*-
"""
pmutt.test_constants
Tests for SMILES_to_thermo_examples file
Created on Tues Oct 23 2018
"""
import unittest
from pgradd.GroupAdd.Library import GroupLibrary
import pgradd.ThermoChem


class TestExamples(unittest.TestCase):

    def test_BensonGA(self):
        lib = GroupLibrary.Load('BensonGA')
        descriptors = lib.GetDescriptors('C1CO1')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'C(C)(H)2(O)': 2, 'O(C)2': 1, 'Oxirane': 1}
        HoRT = thermochem.get_HoRT(298.15)
        self.assertAlmostEqual(HoRT, -21.09467743150278)
        self.assertEqual(descriptors, GroupDict)

    def test_Salciccioli(self):
        lib = GroupLibrary.Load('SalciccioliGA2012')
        descriptors = lib.GetDescriptors('C([Pt])C[Pt]')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'C(C)(H)2(Pt)': 2, 'surface-ring strain': 0.217}
        H = thermochem.get_H(298.15, units='kcal/mol')
        self.assertAlmostEqual(H, -11.307743997749277)
        self.assertEqual(descriptors, GroupDict)

    def test_GU_Aq(self):
        lib = GroupLibrary.Load('GuSolventGA2017Aq')
        descriptors = lib.GetDescriptors('C(=O)([Pt])O')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'CO(O)(Pt)+O(CO)(H)': 1.0}
        HoRT = thermochem.get_HoRT(500)
        self.assertAlmostEqual(HoRT, -109.86212002776878)
        self.assertEqual(descriptors, GroupDict)

    def test_Wittreich_Surface_ND(self):
        lib = GroupLibrary.Load('GRWSurface2018')
        descriptors = lib.GetDescriptors('[Pt]C([Pt])C([Pt])([Pt])C=O')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'C(C)(H)(Pt)2': 1, 'C(C)(CO)(Pt)2': 1, 'CO(C)(H)': 1,
                     'CPt2CPt2': 1, 'CCPt2': 1, 'surface-ring strain': 0.392}
        HoRT = thermochem.get_HoRT(750)
        self.assertAlmostEqual(HoRT, -13.423119203382337)
        self.assertEqual(descriptors, GroupDict)

    def test_Wittreich_Surface_Dim(self):
        lib = GroupLibrary.Load('GRWSurface2018')
        descriptors = lib.GetDescriptors('[Pt]C([Pt])C([Pt])([Pt])C=O')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'C(C)(H)(Pt)2': 1, 'C(C)(CO)(Pt)2': 1, 'CO(C)(H)': 1,
                     'CPt2CPt2': 1, 'CCPt2': 1, 'surface-ring strain': 0.392}
        H = thermochem.get_H(750, 'kcal/mol')
        self.assertAlmostEqual(H, -20.005853103142883)
        self.assertEqual(descriptors, GroupDict)

    def test_Wittreich_Aqueous_ND(self):
        lib = GroupLibrary.Load('GRWAqueous2018')
        descriptors = lib.GetDescriptors('C(=O)([Pt])O')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'CO(O)(Pt)+O(CO)(H)': 1.0}
        HoRT = thermochem.get_HoRT(500)
        self.assertAlmostEqual(HoRT, -107.57909464133714)
        self.assertEqual(descriptors, GroupDict)

    def test_Wittreich_Aqueous_Dim(self):
        lib = GroupLibrary.Load('GRWAqueous2018')
        descriptors = lib.GetDescriptors('C(=O)([Pt])O')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'CO(O)(Pt)+O(CO)(H)': 1.0}
        H = thermochem.get_H(500, 'kJ/mol')
        self.assertAlmostEqual(H, -447.23102885789655)
        self.assertEqual(descriptors, GroupDict)

    def test_Xie_Ru0001_Example1_ND(self):
        lib = GroupLibrary.Load('XieGA2022')
        descriptors = lib.GetDescriptors('[Ru]C([Ru])C([Ru])([Ru])C')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'C(C)(H)(Ru)2': 1, 'C(C)2(Ru)2': 1, 'C(C)(H)3': 1,
                     'CRu2CRu2': 1}
        HoRT = thermochem.get_HoRT(500)
        self.assertAlmostEqual(HoRT, -35.040312149773726)
        self.assertEqual(descriptors, GroupDict)

    def test_Xie_Ru0001_Example1_Dim(self):
        lib = GroupLibrary.Load('XieGA2022')
        descriptors = lib.GetDescriptors('[Ru]C([Ru])C([Ru])([Ru])C')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'C(C)(H)(Ru)2': 1, 'C(C)2(Ru)2': 1, 'C(C)(H)3': 1,
                     'CRu2CRu2': 1}
        H = thermochem.get_H(500, 'kJ/mol')
        self.assertAlmostEqual(H, -145.6706333743726)
        self.assertEqual(descriptors, GroupDict)

    def test_Xie_Ru0001_Example2_ND(self):
        lib = GroupLibrary.Load('XieGA2022')
        descriptors = lib.GetDescriptors('CCC')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'C(C)(H)3': 2, 'C(C)2(H)2': 1}
        HoRT = thermochem.get_HoRT(500)
        self.assertAlmostEqual(HoRT, -41.49969417868688)
        self.assertEqual(descriptors, GroupDict)

    def test_Xie_Ru0001_Example2_Dim(self):
        lib = GroupLibrary.Load('XieGA2022')
        descriptors = lib.GetDescriptors('CCC')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'C(C)(H)3': 2, 'C(C)2(H)2': 1}
        H = thermochem.get_H(500, 'kJ/mol')
        self.assertAlmostEqual(H, -172.52376948049303)
        self.assertEqual(descriptors, GroupDict)

    def test_Selements_Example_Dim(self):
        lib = GroupLibrary.Load('BensonGA')
        descriptors = lib.GetDescriptors('CCCCCC')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'C(C)(H)3': 2, 'C(C)2(H)2': 4}
        G = thermochem.get_G(T=298.15, units='kJ/mol', S_elements=True)
        self.assertAlmostEqual(G, -7.732446702038452)
        self.assertEqual(descriptors, GroupDict)


if __name__ == '__main__':
    unittest.main()
