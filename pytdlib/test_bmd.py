import tdlib
import unittest

from graphs import *
class TestTdLib(unittest.TestCase):

    def test_boost_minDegree_decomp(self):
        V, E, w = tdlib.boost_minDegree_decomp(V_Petersen, E_Petersen)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E), True)
        self.assertEqual(w >= 4, True)
        self.assertEqual(w <= 5, True)
        V, E, w = tdlib.boost_minDegree_decomp(V_Pappus, E_Pappus)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V, E), True)
        self.assertEqual(w, 6)
        V, E, w = tdlib.boost_minDegree_decomp(V_RandomGNM_250_1000, E_RandomGNM_250_1000)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_RandomGNM_250_1000, E_RandomGNM_250_1000, V, E), True)
        self.assertEqual(w, 110)

    def test_boost_minDegree_decomp_decomp_0(self):
        for V, E in cornercases:
            V_T, E_T, w = tdlib.boost_minDegree_decomp(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    def test_boost_minDegree_decomp_1(self):
        V, E, w = tdlib.boost_minDegree_decomp(V_P6, E_P6)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E), True)
        self.assertEqual(w, 1)

    def test_boost_minDegree_decomp_2(self):
        V, E, w = tdlib.boost_minDegree_decomp(V_K5, E_K5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V, E), True)
        self.assertEqual(w, 4)

    def test_boost_minDegree_decomp_3(self):
        V, E, w = tdlib.boost_minDegree_decomp(V_Petersen, E_Petersen)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E), True)
        self.assertEqual(w, 5)

    def test_boost_minDegree_decomp_4(self):
        V, E, w = tdlib.boost_minDegree_decomp(V_Petersen_double, E_Petersen_double)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen_double, E_Petersen_double, V, E), True)
        self.assertEqual(w, 5)

    def test_boost_minDegree_decomp_5(self):
        V, E, w = tdlib.boost_minDegree_decomp(V_Wagner, E_Wagner)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V, E), True)
        self.assertEqual(w, 4)

    def test_boost_minDegree_decomp_6(self):
        V, E, w = tdlib.boost_minDegree_decomp(V_Pappus, E_Pappus)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V, E), True)
        self.assertEqual(w, 6)

    def test_boost_minDegree_decomp_7(self):
        V, E, w = tdlib.boost_minDegree_decomp(V_Grid_5_5, E_Grid_5_5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Grid_5_5, E_Grid_5_5, V, E), True)
        self.assertEqual(w, 5)

    def test_boost_minDegree_decomp_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            V_T, E_T, w = tdlib.boost_minDegree_decomp(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
