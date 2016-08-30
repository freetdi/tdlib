import tdlib
import unittest

from graphs import *

class TestTdLib(unittest.TestCase):
    """ error
    def test_seperator_algorithm_0(self):
        for V, E in cornercases:
            V_T, E_T, w = tdlib.seperator_algorithm(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

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

    def test_minDegree_decomp_0(self):
        for V, E in cornercases:
            V_T, E_T, w = tdlib.minDegree_decomp(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    def test_minDegree_decomp_1(self):
        V, E, w = tdlib.minDegree_decomp(V_P6, E_P6)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E), True)
        self.assertEqual(w, 1)

    def test_minDegree_decomp_2(self):
        V, E, w = tdlib.minDegree_decomp(V_K5, E_K5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V, E), True)
        self.assertEqual(w, 4)

    def test_minDegree_decomp_3(self):
        V, E, w = tdlib.minDegree_decomp(V_Petersen, E_Petersen)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E), True)
        self.assertEqual(w, 4)

    def test_minDegree_decomp_4(self):
        V, E, w = tdlib.minDegree_decomp(V_Petersen_double, E_Petersen_double)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen_double, E_Petersen_double, V, E), True)
        self.assertEqual(w, 4)

    def test_minDegree_decomp_5(self):
        V, E, w = tdlib.minDegree_decomp(V_Wagner, E_Wagner)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V, E), True)
        self.assertEqual(w, 4)

    def test_minDegree_decomp_6(self):
        V, E, w = tdlib.minDegree_decomp(V_Pappus, E_Pappus)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V, E), True)
        self.assertEqual(w, 6)

    def test_minDegree_decomp_7(self):
        V, E, w = tdlib.minDegree_decomp(V_Grid_5_5, E_Grid_5_5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Grid_5_5, E_Grid_5_5, V, E), True)
        self.assertEqual(w, 5)

    def test_minDegree_decomp_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            V_T, E_T, w = tdlib.minDegree_decomp(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    """ segfault
    def test_boost_minDegree_decomp_decomp_0(self):
        for V, E in cornercases:
            V_T, E_T, w = tdlib.boost_minDegree_decomp(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

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

    def test_minDegree_ordering_0(self):
        for V, E in cornercases:
            O = tdlib.minDegree_ordering(V, E)
            self.assertEqual(len(O), len(V))

    def test_minDegree_ordering_1(self):
        O = tdlib.minDegree_ordering(V_P6, E_P6)
        self.assertEqual(len(O), len(V_P6))

    def test_minDegree_ordering_2(self):
        O = tdlib.minDegree_ordering(V_K5, E_K5)
        self.assertEqual(len(O), len(V_K5))

    def test_minDegree_ordering_3(self):
        O = tdlib.minDegree_ordering(V_Petersen, E_Petersen)
        self.assertEqual(O, [0, 2, 6, 3, 5, 1, 4, 7, 8, 9])

    def test_minDegree_ordering_4(self):
        O = tdlib.minDegree_ordering(V_Petersen_double, E_Petersen_double)
        self.assertEqual(len(O), len(V_Petersen_double))

    def test_minDegree_ordering_5(self):
        O = tdlib.minDegree_ordering(V_Wagner, E_Wagner)
        self.assertEqual(len(O), len(V_Wagner))

    def test_minDegree_ordering_6(self):
        O = tdlib.minDegree_ordering(V_Pappus, E_Pappus)
        self.assertEqual(O, [0, 2, 4, 7, 9, 11, 13, 15, 17, 1, 10, 3, \
                             5, 6, 8, 12, 14, 16])

    def test_minDegree_ordering_7(self):
        O = tdlib.minDegree_ordering(V_Grid_5_5, E_Grid_5_5)
        self.assertEqual(len(O), len(V_Grid_5_5))

    def test_minDegree_ordering_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            O = tdlib.minDegree_ordering(V, E)
            self.assertEqual(len(O), len(V))


if __name__ == '__main__':
    unittest.main()
