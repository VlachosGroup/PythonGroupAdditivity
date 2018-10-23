# -*- coding: utf-8 -*-
from VGA.RINGParser.Reader import Read
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
        self.assertEqual(match_index, ((0, 1), (1, 0), (1, 2), (2, 1)))

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
        self.assertEqual(match_index, ((0, 1, 2, 3),))

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
        self.assertEqual(match_index, ((0, 1), (0, 5), (1, 0), (1, 2), (2, 1),
                                       (2, 3), (3, 2), (3, 4), (4, 3), (4, 5),
                                       (5, 0), (5, 4)))

    def test_ring_bond(self):
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
        self.assertEqual(match_index, ((2, 3, 4, 5), (2, 5, 4, 3),
                                       (3, 2, 5, 4), (3, 4, 5, 2),
                                       (4, 3, 2, 5), (4, 5, 2, 3),
                                       (5, 2, 3, 4), (5, 4, 3, 2)))


if __name__ == '__main__':
    unittest.main()
