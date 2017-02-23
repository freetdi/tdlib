import base
import tdlib
import unittest

from graphs import *

class TestTdLib_post(unittest.TestCase):
    def test_MSVS_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            T, w = tdlib.trivial_decomposition(G)
            T, w = tdlib.MSVS(G, T)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_MSVS_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.trivial_decomposition(G)
        T, w = tdlib.MSVS(G, T)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 1)

    def test_MSVS_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.trivial_decomposition(G)
        T, w = tdlib.MSVS(G, T)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_MSVS_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.trivial_decomposition(G)
        T, w = tdlib.MSVS(G, T)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_MSVS_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.trivial_decomposition(G)
        T, w = tdlib.MSVS(G, T)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_MSVS_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.trivial_decomposition(G)
        T, w = tdlib.MSVS(G, T)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_MSVS_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.trivial_decomposition(G)
        T, w = tdlib.MSVS(G, T)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 6)

    def test_MSVS_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.trivial_decomposition(G)
        T, w = tdlib.MSVS(G, T)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 5)

    def test_MSVS_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            G = Graph(V, E)
            T, w = tdlib.trivial_decomposition(G)
            T, w = tdlib.MSVS(G, T)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_minimalChordal_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            T, w = tdlib.trivial_decomposition(G)
            T, w = tdlib.minimalChordal_decomp(G, T)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_minimalChordal_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.trivial_decomposition(G)
        T, w = tdlib.minimalChordal_decomp(G, T)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 1)

    def test_minimalChordal_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.trivial_decomposition(G)
        T, w = tdlib.minimalChordal_decomp(G, T)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_minimalChordal_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.trivial_decomposition(G)
        T, w = tdlib.minimalChordal_decomp(G, T)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_minimalChordal_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.trivial_decomposition(G)
        T, w = tdlib.minimalChordal_decomp(G, T)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_minimalChordal_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.trivial_decomposition(G)
        T, w = tdlib.minimalChordal_decomp(G, T)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_minimalChordal_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.trivial_decomposition(G)
        T, w = tdlib.minimalChordal_decomp(G, T)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 6)

    def test_minimalChordal_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.trivial_decomposition(G)
        T, w = tdlib.minimalChordal_decomp(G, T)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 5)

    def test_minimalChordal_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            G = Graph(V, E)
            T, w = tdlib.trivial_decomposition(G)
            T, w = tdlib.minimalChordal_decomp(G, T)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
