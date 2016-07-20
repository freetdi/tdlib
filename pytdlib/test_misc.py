import tdlib
import unittest

from graphs import *

class TestTdLib_misc(unittest.TestCase):

    def test_treedec_to_ordering(self):
        V, E, lb = tdlib.seperator_algorithm(V_P6, E_P6)
        O = tdlib.treedec_to_ordering(V, E)
        self.assertEqual(O, [0, 2, 1, 3, 4, 5])

        V = [["a", "d", "c"], ["c", "b", "e"]]
        E = [0, 1]
        O = tdlib.treedec_to_ordering(V, E)
        self.assertEqual(O, ["a", "d", "c", "b", "e"])

    def test_ordering_to_treedec(self):
        V, E, lb = tdlib.ordering_to_treedec(V_P6, E_P6, [1,3,0,2,4,5])
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E) == 0, True)

    def test_trivial_decomposition(self):
        V, E = tdlib.trivial_decomposition(["a", "b", "c"], [])
        self.assertEqual(V, [["a", "b", "c"]])

        V, E = tdlib.trivial_decomposition(V_P6, E_P6)
        self.assertEqual(V, [[0,1,2,3,4,5]])

        V, E = tdlib.trivial_decomposition(V_K5, E_K5)
        self.assertEqual(V, [[0,1,2,3,4]])

    def test_is_valid_treedecomposition(self):
        V, E = tdlib.trivial_decomposition(V_P6, E_P6)

        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E) == 0, True)

        status = tdlib.is_valid_treedecomposition(V_K5, E_K5, V, E)
        self.assertEqual(status, -6)

        V, E, tw = tdlib.exact_decomposition_cutset(V_Petersen, E_Petersen)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E) == 0, True)

        V, E, tw = tdlib.exact_decomposition_cutset(V_Wagner, E_Wagner)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V, E) == 0, True)


if __name__ == '__main__':
    unittest.main()
