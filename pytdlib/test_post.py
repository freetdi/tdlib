import tdlib
import unittest

from graphs import *

class TestTdLib_post(unittest.TestCase):
    def test_MSVS_0(self):
        for V, E in cornercases:
            V_T, E_T = tdlib.trivial_decomposition(V, E)
            V_T, E_T, w = tdlib.MSVS(V, E, V_T, E_T)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    def test_MSVS_1(self):
        V_T, E_T = tdlib.trivial_decomposition(V_P6, E_P6)
        V_T, E_T, w = tdlib.MSVS(V_P6, E_P6, V_T, E_T)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V_T, E_T), True)
        self.assertEqual(w, 1)

    def test_MSVS_2(self):
        V_T, E_T = tdlib.trivial_decomposition(V_K5, E_K5)
        V_T, E_T, w = tdlib.MSVS(V_K5, E_K5, V_T, E_T)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V_T, E_T), True)
        self.assertEqual(w, 4)

    def test_MSVS_3(self):
        V_T, E_T = tdlib.trivial_decomposition(V_Petersen, E_Petersen)
        V_T, E_T, w = tdlib.MSVS(V_Petersen, E_Petersen, V_T, E_T)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V_T, E_T), True)
        self.assertEqual(w, 4)

    def test_MSVS_4(self):
        V, E = tdlib.trivial_decomposition(V_Petersen_double, E_Petersen_double)
        V, E, w = tdlib.MSVS(V_Petersen_double, E_Petersen_double, V, E)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen_double, E_Petersen_double, V, E), True)
        self.assertEqual(w, 4)

    def test_MSVS_5(self):
        V, E = tdlib.trivial_decomposition(V_Wagner, E_Wagner)
        V, E, w = tdlib.MSVS(V_Wagner, E_Wagner, V, E)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V, E), True)
        self.assertEqual(w, 4)

    def test_MSVS_6(self):
        V_T, E_T = tdlib.trivial_decomposition(V_Pappus, E_Pappus)
        V_T, E_T, w = tdlib.MSVS(V_Pappus, E_Pappus, V_T, E_T)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V_T, E_T), True)
        self.assertEqual(w, 6)

    def test_MSVS_7(self):
        V_T, E_T = tdlib.trivial_decomposition(V_Grid_5_5, E_Grid_5_5)
        V_T, E_T, w = tdlib.MSVS(V_Grid_5_5, E_Grid_5_5, V_T, E_T)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Grid_5_5, E_Grid_5_5, V_T, E_T), True)
        self.assertEqual(w, 5)

    def test_MSVS_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            V_T, E_T = tdlib.trivial_decomposition(V, E)
            V_T, E_T, w = tdlib.MSVS(V, E, V_T, E_T)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    def test_minimalChordal_0(self):
        for V, E in cornercases:
            V_T, E_T = tdlib.trivial_decomposition(V, E)
            V_T, E_T, w = tdlib.minimalChordal_decomp(V, E, V_T, E_T)
            #O = tdlib.minimalChordal_decomp(V, E, range(0, len(V)))
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    def test_minimalChordal_1(self):
        V_T, E_T = tdlib.trivial_decomposition(V_P6, E_P6)
        V_T, E_T, w = tdlib.minimalChordal_decomp(V_P6, E_P6, V_T, E_T)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V_T, E_T), True)
        self.assertEqual(w, 1)

    def test_minimalChordal_2(self):
        V_T, E_T = tdlib.trivial_decomposition(V_K5, E_K5)
        V_T, E_T, w = tdlib.minimalChordal_decomp(V_K5, E_K5, V_T, E_T)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V_T, E_T), True)
        self.assertEqual(w, 4)

    def test_minimalChordal_3(self):
        V_T, E_T = tdlib.trivial_decomposition(V_Petersen, E_Petersen)
        V_T, E_T, w = tdlib.minimalChordal_decomp(V_Petersen, E_Petersen, V_T, E_T)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V_T, E_T), True)
        self.assertEqual(w, 4)

    def test_minimalChordal_4(self):
        V, E = tdlib.trivial_decomposition(V_Petersen_double, E_Petersen_double)
        V, E, w = tdlib.minimalChordal_decomp(V_Petersen_double, E_Petersen_double, V, E)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen_double, E_Petersen_double, V, E), True)
        self.assertEqual(w, 4)

    def test_minimalChordal_5(self):
        V, E = tdlib.trivial_decomposition(V_Wagner, E_Wagner)
        V, E, w = tdlib.minimalChordal_decomp(V_Wagner, E_Wagner, V, E)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V, E), True)
        self.assertEqual(w, 4)

    def test_minimalChordal_6(self):
        V_T, E_T = tdlib.trivial_decomposition(V_Pappus, E_Pappus)
        V_T, E_T, w = tdlib.minimalChordal_decomp(V_Pappus, E_Pappus, V_T, E_T)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V_T, E_T), True)
        self.assertEqual(w, 6)

    def test_minimalChordal_7(self):
        V_T, E_T = tdlib.trivial_decomposition(V_Grid_5_5, E_Grid_5_5)
        V_T, E_T, w = tdlib.minimalChordal_decomp(V_Grid_5_5, E_Grid_5_5, V_T, E_T)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Grid_5_5, E_Grid_5_5, V_T, E_T), True)
        self.assertEqual(w, 5)

    def test_minimalChordal_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            V_T, E_T = tdlib.trivial_decomposition(V, E)
            V_T, E_T, w = tdlib.minimalChordal_decomp(V, E, V_T, E_T)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)


if __name__ == '__main__':
    unittest.main()
