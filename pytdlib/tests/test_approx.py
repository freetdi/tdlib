import test_base
import tdlib
import unittest
import util

from graphs import *
class TestTdLib(unittest.TestCase):

    def test_minDegree_decomp_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            T, w = tdlib.minDegree_decomp(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_minDegree_decomp_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 1)

    def test_minDegree_decomp_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_minDegree_decomp_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 5) # could be 4?

    def test_minDegree_decomp_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 5) # could be 4?

    def test_minDegree_decomp_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_minDegree_decomp_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 6)

    def test_minDegree_decomp_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 5)

    def test_minDegree_decomp_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            G = Graph(V, E)
            T, w = tdlib.minDegree_decomp(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_minDegree_ordering_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            O = tdlib.minDegree_ordering(G)
            self.assertEqual(len(O), len(G.vertices()))

    def test_minDegree_ordering_1(self):
        G = Graph(V_P6, E_P6)
        O = tdlib.minDegree_ordering(G)
        self.assertEqual(len(O), len(G.vertices()))

    def test_minDegree_ordering_2(self):
        G = Graph(V_K5, E_K5)
        O = tdlib.minDegree_ordering(G)
        self.assertEqual(len(O), len(G.vertices()))

    def test_minDegree_ordering_3(self):
        G = Graph(V_Petersen, E_Petersen)
        O = tdlib.minDegree_ordering(G)
        self.assertEqual(util.is_permutation(O), True)

    def test_minDegree_ordering_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        O = tdlib.minDegree_ordering(G)
        self.assertEqual(len(O), len(G.vertices()))

    def test_minDegree_ordering_5(self):
        G = Graph(V_Wagner, E_Wagner)
        O = tdlib.minDegree_ordering(G)
        self.assertEqual(len(O), len(G.vertices()))

    def test_minDegree_ordering_6(self):
        G = Graph(V_Pappus, E_Pappus)
        O = tdlib.minDegree_ordering(G)
        self.assertEqual(util.is_permutation(O), True)

    def test_minDegree_ordering_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        O = tdlib.minDegree_ordering(G)
        self.assertEqual(len(O), len(G.vertices()))

    def test_minDegree_ordering_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            G = Graph(V, E)
            O = tdlib.minDegree_ordering(G)
            self.assertEqual(len(O), len(G.vertices()))

if __name__ == '__main__':
    unittest.main()
