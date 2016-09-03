import tdlib
import unittest

from graphs import *

class TestTdLib(unittest.TestCase):
    def test_fillIn_decomp_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            T, w = tdlib.fillIn_decomp(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_fillIn_decomp_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.fillIn_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 1)

    def test_fillIn_decomp_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.fillIn_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_fillIn_decomp_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.fillIn_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_fillIn_decomp_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.fillIn_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_fillIn_decomp_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.fillIn_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_fillIn_decomp_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.fillIn_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 6)

    def test_fillIn_decomp_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.fillIn_decomp(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 5)

    def test_fillIn_decomp_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            G = Graph(V, E)
            T, w = tdlib.fillIn_decomp(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_fillIn_ordering_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            O = tdlib.fillIn_ordering(G)
            self.assertEqual(len(O), len(G.vertices()))

    def test_fillIn_ordering_1(self):
        G = Graph(V_P6, E_P6)
        O = tdlib.fillIn_ordering(G)
        self.assertEqual(len(O), len(G.vertices()))

    def test_fillIn_ordering_2(self):
        G = Graph(V_K5, E_K5)
        O = tdlib.fillIn_ordering(G)
        self.assertEqual(len(O), len(G.vertices()))

    def test_fillIn_ordering_3(self):
        G = Graph(V_Petersen, E_Petersen)
        O = tdlib.fillIn_ordering(G)
        self.assertEqual(O, [0, 2, 6, 3, 5, 1, 4, 7, 8, 9])

    def test_fillIn_ordering_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        O = tdlib.fillIn_ordering(G)
        self.assertEqual(len(O), len(G.vertices()))

    def test_fillIn_ordering_5(self):
        G = Graph(V_Wagner, E_Wagner)
        O = tdlib.fillIn_ordering(G)
        self.assertEqual(len(O), len(G.vertices()))

    def test_fillIn_ordering_6(self):
        G = Graph(V_Pappus, E_Pappus)
        O = tdlib.fillIn_ordering(G)
        self.assertEqual(O, [0, 2, 4, 7, 9, 11, 13, 15, 17, 1, \
                             10, 3, 5, 6, 8, 12, 14, 16])

    def test_fillIn_ordering_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        O = tdlib.fillIn_ordering(G)
        self.assertEqual(len(O), len(G.vertices()))

    def test_fillIn_ordering_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            G = Graph(V, E)
            O = tdlib.fillIn_ordering(G)
            self.assertEqual(len(O), len(G.vertices()))



if __name__ == '__main__':
    unittest.main()
