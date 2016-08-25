import tdlib
import unittest

from graphs import *

class TestTdLib_post(unittest.TestCase):

    def test_MSVS(self):
        V, E = tdlib.trivial_decomposition(V_P6, E_P6)
        V, E, w = tdlib.MSVS(V_P6, E_P6, V, E)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E) == 0, True)
        self.assertEqual(w, 1)

        V, E = tdlib.trivial_decomposition(V_Petersen, E_Petersen)
        V, E, w = tdlib.MSVS(V_Petersen, E_Petersen, V, E)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E) == 0, True)
        self.assertEqual(w, 4)

        V, E = tdlib.trivial_decomposition(V_Pappus, E_Pappus)
        V, E, w = tdlib.MSVS(V_Pappus, E_Pappus, V, E)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V, E) == 0, True)
        self.assertEqual(w, 6)

    def test_minimalChordal(self):
        V1, E1, lb1 = tdlib.ordering_to_treedec(V_Pappus, E_Pappus, [i for i in range(0, len(V_Pappus))])
        self.assertEqual(lb1, 8)

        O2 = tdlib.minimalChordal_ordering(V_Pappus, E_Pappus, [0, 1, 2, 3, 4, 5, 6, 7, 8, \
                                                  9, 10, 11, 12, 13, 14, 15, 16, 17])
        V_T2, E_T2, w2 = tdlib.ordering_to_treedec(V_Pappus, E_Pappus, O2)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V_T2, E_T2) == 0, True)
        self.assertEqual(w2, 6)


if __name__ == '__main__':
    unittest.main()
