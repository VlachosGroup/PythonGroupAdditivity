# -*- coding: utf-8 -*-
from pgradd.RINGParser.Reader import Read
from rdkit import Chem
import unittest


class TestRINGParser(unittest.TestCase):

    def test_molquery(self):
        testmol = Chem.MolFromSmiles('CCC')
        s = """
        fragment a{
            C labeled c1
            C labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((1, 0), (0, 1), (1, 2), (2, 1))))

    def test_double_triple_bond(self):
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
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1, 2, 3),)))

    def test_aromatic_bond(self):
        testmol = Chem.MolFromSmiles('c1ccccc1')
        s = """
        fragment a{
            C labeled c1
            C labeled c2 aromatic bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1), (0, 5), (1, 0), (1, 2), (2, 1),
                                     (2, 3), (3, 2), (3, 4), (4, 3), (4, 5),
                                     (5, 4), (5, 0))))

    def test_ring_bond1(self):
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
        self.assertListEqual(sorted(match_index),
		                     sorted(((2, 3, 4, 5), (2, 5, 4, 3),
                                     (3, 2, 5, 4), (3, 4, 5, 2),
                                     (4, 3, 2, 5), (4, 5, 2, 3),
                                     (5, 2, 3, 4), (5, 4, 3, 2))))

    def test_ring_bond2(self):
        testmol = Chem.MolFromSmiles('CCC1CCC1')
        s = """
        fragment a{
            C labeled c1
            C labeled c2 ring bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((2, 3), (2, 5), (3, 2), (3, 4), (4, 3),
                                     (4, 5), (5, 4), (5, 2))))

    def test_non_ring_bond(self):
        testmol = Chem.MolFromSmiles('CCC1CCC1')
        s = """
        fragment a{
            C labeled c1
            C labeled c2 nonring bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1), (1, 0), (1, 2), (2, 1))))

    def test_any_bond1(self):
        testmol = Chem.MolFromSmiles('CC=CC#C')
        s = """
        fragment a{
            C labeled c1
            C labeled c2 any bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1), (1, 0), (1, 2), (2, 1),
                                     (2, 3), (3, 2), (3, 4), (4, 3))))

    def test_any_bond2(self):
        testmol = Chem.MolFromSmiles('CC=CC#C')
        s = """
        fragment a{
            C labeled c1
            C labeled c2 any bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1), (1, 0), (1, 2), (2, 1),
                                     (2, 3), (3, 2), (3, 4), (4, 3))))

    def test_strong_bond(self):
        testmol = Chem.MolFromSmiles('CC=CC#C')
        s = """
        fragment a{
            C labeled c1
            C labeled c2 strong bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((1, 2), (2, 1), (3, 4), (4, 3))))

    def test_other_bond1(self):
        testmol = Chem.MolFromSmiles('C[CH2-]')
        s = """
        positive fragment a{
            C labeled c1
            C labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertEqual(match_index, ())

    def test_other_bond2(self):
        testmol = Chem.MolFromSmiles('C[CH2+]')
        s = """
        negative fragment a{
            C labeled c1
            C labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertEqual(match_index, ())

    def test_other_bond3(self):
        testmol = Chem.MolFromSmiles('C=CC')
        s = """
        olefinic fragment a{
            C labeled c1
            C labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((1, 2), (2, 1))))

    def test_other_bond4(self):
        testmol = Chem.MolFromSmiles('C=CC')
        s = """
        paraffinic fragment a{
            C labeled c1
            C labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertEqual(match_index, ())

    def test_other_bond5(self):
        testmol = Chem.MolFromSmiles('CCC')
        s = """
        paraffinic fragment a{
            C labeled c1
            C labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1), (1, 0), (1, 2), (2, 1))))

    def test_other_bond6(self):
        testmol = Chem.MolFromSmiles('CCC')
        s = """
        linear fragment a{
            C labeled c1
            C labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1), (1, 0), (1, 2), (2, 1))))

    def test_other_bond7(self):
        testmol = Chem.MolFromSmiles('CCC')
        s = """
        cyclic fragment a{
            C labeled c1
            C labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertEqual(match_index, ())

    def test_other_bond8(self):
        testmol = Chem.MolFromSmiles('C1CCC1C')
        s = """
        cyclic fragment a{
            C labeled c1
            C labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1), (0, 3), (1, 0), (1, 2),
                                     (2, 1), (2, 3), (3, 2),
                                     (3, 4), (3, 0), (4, 3))))

    def test_symbol_atomsuffix(self):
        testmol = Chem.MolFromSmiles('CCC')
        s = """
        fragment a{
            $ labeled c1
            $ labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1), (0, 3), (0, 4), (0, 5), (1, 0),
                                     (1, 2), (1, 6), (1, 7), (2, 1), (2, 8),
                                     (2, 9), (2, 10), (3, 0), (4, 0),
                                     (5, 0), (6, 1), (7, 1), (8, 2), (9, 2),
                                     (10, 2))))

    def test_other_bond9(self):
        testmol = Chem.MolFromSmiles('CCO')
        s = """
        fragment a{
            X labeled c1
            X labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1), (1, 0), (1, 2), (2, 1))))

    def test_other_bond10(self):
        testmol = Chem.MolFromSmiles('CCS')
        s = """
        fragment a{
            X labeled c1
            & labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((1, 2),)))

    def test_atom_constraint1(self):
        testmol = Chem.MolFromSmiles('CCC')
        s = """
        fragment a{
            C labeled c1
            C labeled c2 single bond to c1 {connected to =2 C}
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1), (2, 1))))

    def test_atom_constraint2(self):
        testmol = Chem.MolFromSmiles('CCC')
        s = """
        fragment a{
            C labeled c1
            C labeled c2 single bond to c1 {connected to =1 C}
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((1, 0), (1, 2))))

    def test_atom_constraint3(self):
        testmol = Chem.MolFromSmiles('CC=C')
        s = """
        fragment a{
            C labeled c1
            C labeled c2 single bond to c1 {connected to =1 C with double bond}
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1),)))

    def test_atom_constraint4(self):
        testmol = Chem.MolFromSmiles('CC=C')
        s = """
        fragment a{
            C labeled c1 {connected to >1 C with any bond}
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((1,),)))

    def test_atom_constraint5(self):
        testmol = Chem.MolFromSmiles('CC=C')
        s = """
        fragment a{
            C labeled c1 {!connected to >1 C with any bond}
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, ), (2, ))))

    def test_atom_constraint6(self):
        testmol = Chem.MolFromSmiles('CC1CCC1')
        s = """
        fragment a{
            C labeled c1 {!in ring of size >0}
            C labeled c2 single bond to c1 {in ring of size >0}
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1),)))

    def test_atom_prefix1(self):
        testmol = Chem.MolFromSmiles('CC1CCC1')
        s = """
        fragment a{
            nonringatom C labeled c1
            ringatom C labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1),)))

    def test_atom_prefix2(self):
        testmol = Chem.MolFromSmiles('Cc1ccccc1')
        s = """
        fragment a{
            nonaromatic C labeled c1
            aromatic C labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1),)))

    def test_atom_prefix3(self):
        testmol = Chem.MolFromSmiles('CC=C')
        s = """
        fragment a{
            C labeled c1
            allylic C labeled c2 single bond to c1
        }
        """
        molquery = Read(s)
        match_index = molquery.GetQueryMatches(testmol)
        self.assertListEqual(sorted(match_index),
                             sorted(((0, 1),)))


if __name__ == '__main__':
    unittest.main()
