import base
import sys
import tdlib
import unittest

from graphs import *

#don't confuse python unittest
sys.argv=sys.argv[:1]

class TestTdLib(unittest.TestCase):
    def test_seperator_algorithm_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            T, w = tdlib.seperator_algorithm(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_seperator_algorithm_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.seperator_algorithm(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w >= 1 and w <= 5, True)

    def test_seperator_algorithm_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.seperator_algorithm(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w >= 4 and w <= 17, True)

    def test_seperator_algorithm_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.seperator_algorithm(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w >= 4 and w <= 17, True)

    def test_seperator_algorithm_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.seperator_algorithm(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w >= 4 and w <= 17, True)

    def test_seperator_algorithm_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.seperator_algorithm(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w >= 4 and w <= 17, True)

    def test_seperator_algorithm_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.seperator_algorithm(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w >= 6 and w <= 25, True)

    def test_seperator_algorithm_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.seperator_algorithm(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w >= 5 and w <= 21, True)

    def test_seperator_algorithm_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            G = Graph(V, E)
            T, w = tdlib.seperator_algorithm(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
