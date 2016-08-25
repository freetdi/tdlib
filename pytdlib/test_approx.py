import tdlib
import unittest

from graphs import *

class TestTdLib(unittest.TestCase):
    def test_seperator_algorithm(self):
        V, E, w = tdlib.seperator_algorithm(V_P6, E_P6)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E) == 0, True)
        self.assertEqual(w, 3)

        V, E, w = tdlib.seperator_algorithm(V_K5, E_K5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V, E) == 0, True)
        self.assertEqual(w, 4)

        V, E, lb = tdlib.seperator_algorithm(V_Pappus, E_Pappus)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V, E) == 0, True)
        self.assertEqual(lb, 8)

    def test_minDegree_decomp(self):
        V, E, w = tdlib.minDegree_decomp(V_Petersen, E_Petersen)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E) == 0, True)
        self.assertEqual(w >= 4, True)
        self.assertEqual(w <= 5, True)
        V, E, w = tdlib.minDegree_decomp(V_Pappus, E_Pappus)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V, E) == 0, True)
        self.assertEqual(w, 6)


if __name__ == '__main__':
    unittest.main()
