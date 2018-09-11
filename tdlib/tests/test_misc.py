import base
import sys
import tdlib
import unittest

from graphs import *

#don't confuse python unittest
sys.argv=sys.argv[:1]

class TestTdLib_misc(unittest.TestCase):
    def test_conversion_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            T, w = tdlib.PP_MD(G)
            O = tdlib.treedec_to_ordering(T)
            T, w = tdlib.ordering_to_treedec(G, O)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_conversion_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.PP_MD(G)
        O = tdlib.treedec_to_ordering(T)
        T, w = tdlib.ordering_to_treedec(G, O)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_conversion_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.PP_MD(G)
        O = tdlib.treedec_to_ordering(T)
        T, w = tdlib.ordering_to_treedec(G, O)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_conversion_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.PP_MD(G)
        O = tdlib.treedec_to_ordering(T)
        T, w = tdlib.ordering_to_treedec(G, O)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_conversion_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.PP_MD(G)
        O = tdlib.treedec_to_ordering(T)
        T, w = tdlib.ordering_to_treedec(G, O)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_conversion_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.PP_MD(G)
        O = tdlib.treedec_to_ordering(T)
        T, w = tdlib.ordering_to_treedec(G, O)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_conversion_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.PP_MD(G)
        O = tdlib.treedec_to_ordering(T)
        T, w = tdlib.ordering_to_treedec(G, O)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_conversion_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            G = Graph(V, E)
            T, w = tdlib.PP_MD(G)
            self.assertTrue(tdlib.is_valid_treedecomposition(G, T))
            O = tdlib.treedec_to_ordering(T)
            T, w = tdlib.ordering_to_treedec(G, O)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)


    def test_trivial_decomposition_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            T, w = tdlib.trivial_decomposition(G)
            self.assertEqual(T.vertices(), [list(V)])

    def test_trivial_decomposition_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.trivial_decomposition(G)
        self.assertEqual(T.vertices(), [list(V_P6)])

    def test_trivial_decomposition_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.trivial_decomposition(G)
        self.assertEqual(T.vertices(), [list(V_K5)])

    def test_trivial_decomposition_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.trivial_decomposition(G)
        self.assertEqual(T.vertices(), [list(V_Petersen)])

    def test_trivial_decomposition_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.trivial_decomposition(G)
        self.assertEqual(T.vertices(), [list(V_Petersen_double)])

    def test_trivial_decomposition_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.trivial_decomposition(G)
        self.assertEqual(T.vertices(), [list(V_Wagner)])

    def test_trivial_decomposition_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.trivial_decomposition(G)
        self.assertEqual(T.vertices(), [list(V_Pappus)])

    def test_trivial_decomposition_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.trivial_decomposition(G)
        self.assertEqual(T.vertices(), [list(V_Grid_5_5)])

    def test_trivial_decomposition_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            G = Graph(V, E)
            T, w = tdlib.trivial_decomposition(G)
            self.assertEqual(T.vertices(), [V])

    """ cython error
    def test_trivial_decomposition_9(self):
        G = Graph(["a", "b", "c"], [])
        T, w = tdlib.trivial_decomposition(G)
        self.assertEqual(T.vertices(), [G.vertices()])
    """

    def test_is_valid_treedecomposition_0(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.trivial_decomposition(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_is_valid_treedecomposition_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.trivial_decomposition(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(Graph(V_K5, E_K5), T, False), False)

    def test_is_valid_treedecomposition_2(self):
        G = Graph(V_Petersen, E_Petersen)
        T, tw = tdlib.exact_decomposition_cutset(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_is_valid_treedecomposition_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, tw = tdlib.exact_decomposition_cutset(G)
        del T.vertices()[-1]
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T, False), False)

    def test_is_valid_treedecomposition_4(self):
        G = Graph(V_Wagner, E_Wagner)
        T, tw = tdlib.exact_decomposition_cutset(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_is_valid_treedecomposition_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, tw = tdlib.exact_decomposition_cutset(G)
        del T.vertices()[-1]
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T, False), False)

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
