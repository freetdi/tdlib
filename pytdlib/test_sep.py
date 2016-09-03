import tdlib
import unittest

from graphs import *

class TestTdLib(unittest.TestCase):
    def test_seperator_algorithm_0(self):
        for V, E in cornercases:
            V_T, E_T, w = tdlib.seperator_algorithm(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    def test_seperator_algorithm_1(self):
        V, E, w = tdlib.seperator_algorithm(V_P6, E_P6)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E), True)
        self.assertEqual(w >= 1 and w <= 5, True)

    def test_seperator_algorithm_2(self):
        V, E, w = tdlib.seperator_algorithm(V_K5, E_K5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V, E), True)
        self.assertEqual(w >= 4 and w <= 17, True)

    def test_seperator_algorithm_3(self):
        V, E, w = tdlib.seperator_algorithm(V_Petersen, E_Petersen)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E), True)
        self.assertEqual(w >= 4 and w <= 17, True)

    def test_seperator_algorithm_4(self):
        V, E, w = tdlib.seperator_algorithm(V_Petersen_double, E_Petersen_double)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen_double, E_Petersen_double, V, E), True)
        self.assertEqual(w >= 4 and w <= 17, True)

    def test_seperator_algorithm_5(self):
        V, E, w = tdlib.seperator_algorithm(V_Wagner, E_Wagner)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V, E), True)
        self.assertEqual(w >= 4 and w <= 17, True)

    def test_seperator_algorithm_6(self):
        V, E, w = tdlib.seperator_algorithm(V_Pappus, E_Pappus)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V, E), True)
        self.assertEqual(w >= 6 and w <= 25, True)

    def test_seperator_algorithm_7(self):
        V, E, w = tdlib.seperator_algorithm(V_Grid_5_5, E_Grid_5_5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Grid_5_5, E_Grid_5_5, V, E), True)
        self.assertEqual(w >= 5 and w <= 21, True)

    def test_seperator_algorithm_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            V_T, E_T, w = tdlib.seperator_algorithm(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
