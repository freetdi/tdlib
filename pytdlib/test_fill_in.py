import tdlib
import unittest

from graphs import *

class TestTdLib(unittest.TestCase):
    def test_fillIn_decomp_peter(self):
        V, E, w = tdlib.fillIn_decomp(V_Petersen, E_Petersen)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E), True)
        self.assertEqual(w, 4)

    def test_fillIn_decomp_papp(self):
        V, E, w = tdlib.fillIn_decomp(V_Pappus, E_Pappus)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V, E), True)
        self.assertEqual(w, 6)

    def test_fillIn_ordering(self):
        O = tdlib.fillIn_ordering(V_Petersen, E_Petersen)
        self.assertEqual(O, [0, 2, 6, 3, 5, 1, 4, 7, 8, 9])
        O = tdlib.fillIn_ordering(V_Pappus, E_Pappus)
        self.assertEqual(O, [0, 2, 4, 7, 9, 11, 13, 15, 17, 1, \
                             10, 3, 5, 6, 8, 12, 14, 16])

if __name__ == '__main__':
    unittest.main()
