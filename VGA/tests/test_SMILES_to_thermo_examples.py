# -*- coding: utf-8 -*-
"""
pMuTT.test_constants
Tests for SMILES_to_thermo_examples file
Created on Tues Oct 23 2018
"""
import unittest
from VGA.GroupAdd.Library import GroupLibrary
import VGA.ThermoChem


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
        HoRT = thermochem.get_HoRT(298.15)
        self.assertAlmostEqual(HoRT, 37.62494617247582)
        self.assertEqual(descriptors, GroupDict)

    def test_GU_Aq(self):
        lib = GroupLibrary.Load('GuSolventGA2017Aq')
        descriptors = lib.GetDescriptors('C(=O)([Pt])O')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'CO(O)(Pt)+O(CO)(H)': 1.0}
        HoRT = thermochem.get_HoRT(500)
        self.assertAlmostEqual(HoRT, -109.86212002776878)
        self.assertEqual(descriptors, GroupDict)

    def test_Wittreich_Surface(self):
        lib = GroupLibrary.Load('GRWSurface2018')
        descriptors = lib.GetDescriptors('[Pt]C([Pt])C([Pt])([Pt])C=O')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'C(C)(H)(Pt)2': 1, 'C(C)(CO)(Pt)2': 1, 'CO(C)(H)': 1,
                     'CPt2CPt2': 1, 'CCPt2': 1, 'surface-ring strain': 0.392}
        HoRT = thermochem.get_HoRT(750)
        self.assertAlmostEqual(HoRT, -13.423119203382337)
        self.assertEqual(descriptors, GroupDict)

    def test_Wittreich_Aqueous(self):
        lib = GroupLibrary.Load('GRWAqueous2018')
        descriptors = lib.GetDescriptors('C(=O)([Pt])O')
        thermochem = lib.Estimate(descriptors, 'thermochem')
        GroupDict = {'CO(O)(Pt)+O(CO)(H)': 1.0}
        HoRT = thermochem.get_HoRT(500)
        self.assertAlmostEqual(HoRT, -107.57909464133714)
        self.assertEqual(descriptors, GroupDict)


if __name__ == '__main__':
    unittest.main()
