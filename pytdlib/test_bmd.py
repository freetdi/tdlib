import tdlib
import unittest

from graphs import *
class TestTdLib(unittest.TestCase):

    def test_boost_minDegree_decomp(self):
        G = Graph(V_RandomGNM_250_1000, E_RandomGNM_250_1000)
        T, w = tdlib.boost_minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 110)

    def test_boost_minDegree_decomp_decomp_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            T, w = tdlib.boost_minDegree_decomp(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_boost_minDegree_decomp_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.boost_minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 1)

    def test_boost_minDegree_decomp_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.boost_minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_boost_minDegree_decomp_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.boost_minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 5)

    def test_boost_minDegree_decomp_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.boost_minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 5)

    def test_boost_minDegree_decomp_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.boost_minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_boost_minDegree_decomp_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.boost_minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 6)

    def test_boost_minDegree_decomp_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.boost_minDegree_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 5)

    def test_boost_minDegree_decomp_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            G = Graph(V, E)
            T, w = tdlib.boost_minDegree_decomp(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
