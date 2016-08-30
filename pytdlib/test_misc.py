import tdlib
import unittest

from graphs import *

class TestTdLib_misc(unittest.TestCase):

    def test_conversion_0(self):
        for V, E in cornercases:
            V_T, E_T, w = tdlib.PP_MD(V, E)
            O = tdlib.treedec_to_ordering(V_T, E_T)
            V_T, E_T, w = tdlib.ordering_to_treedec(V, E, O)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    def test_conversion_1(self):
        V_T, E_T, w = tdlib.PP_MD(V_P6, E_P6)
        O = tdlib.treedec_to_ordering(V_T, E_T)
        V_T, E_T, w = tdlib.ordering_to_treedec(V_P6, E_P6, O)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V_T, E_T), True)

    def test_conversion_2(self):
        V_T, E_T, w = tdlib.PP_MD(V_K5, E_K5)
        O = tdlib.treedec_to_ordering(V_T, E_T)
        V_T, E_T, w = tdlib.ordering_to_treedec(V_K5, E_K5, O)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V_T, E_T), True)

    def test_conversion_3(self):
        V_T, E_T, w = tdlib.PP_MD(V_Petersen, E_Petersen)
        O = tdlib.treedec_to_ordering(V_T, E_T)
        V_T, E_T, w = tdlib.ordering_to_treedec(V_Petersen, E_Petersen, O)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V_T, E_T), True)

    def test_conversion_4(self):
        V_T, E_T, w = tdlib.PP_MD(V_Petersen_double, E_Petersen_double)
        O = tdlib.treedec_to_ordering(V_T, E_T)
        V_T, E_T, w = tdlib.ordering_to_treedec(V_Petersen_double, E_Petersen_double, O)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen_double, E_Petersen_double, V_T, E_T), True)

    def test_conversion_5(self):
        V_T, E_T, w = tdlib.PP_MD(V_Wagner, E_Wagner)
        O = tdlib.treedec_to_ordering(V_T, E_T)
        V_T, E_T, w = tdlib.ordering_to_treedec(V_Wagner, E_Wagner, O)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V_T, E_T), True)

    def test_conversion_6(self):
        V_T, E_T, w = tdlib.PP_MD(V_Pappus, E_Pappus)
        O = tdlib.treedec_to_ordering(V_T, E_T)
        V_T, E_T, w = tdlib.ordering_to_treedec(V_Pappus, E_Pappus, O)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V_T, E_T), True)

    def test_conversion_7(self):
        V_T, E_T, w = tdlib.PP_MD(V_Grid_5_5, E_Grid_5_5)
        O = tdlib.treedec_to_ordering(V_T, E_T)
        V_T, E_T, w = tdlib.ordering_to_treedec(V_Grid_5_5, E_Grid_5_5, O)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Grid_5_5, E_Grid_5_5, V_T, E_T), True)

    def test_conversion_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            V_T, E_T, w = tdlib.PP_MD(V, E)
            O = tdlib.treedec_to_ordering(V_T, E_T)
            V_T, E_T, w = tdlib.ordering_to_treedec(V, E, O)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)


    def test_trivial_decomposition_0(self):
        for V, E in cornercases:
            V_T, E_T = tdlib.trivial_decomposition(V, E)
            self.assertEqual(V_T, [V])

    def test_trivial_decomposition_1(self):
        V, E = tdlib.trivial_decomposition(V_P6, E_P6)
        self.assertEqual(V, [V_P6])

    def test_trivial_decomposition_2(self):
        V, E = tdlib.trivial_decomposition(V_K5, E_K5)
        self.assertEqual(V, [V_K5])

    def test_trivial_decomposition_3(self):
        V, E = tdlib.trivial_decomposition(V_Petersen, E_Petersen)
        self.assertEqual(V, [V_Petersen])

    def test_trivial_decomposition_4(self):
        V, E = tdlib.trivial_decomposition(V_Petersen_double, E_Petersen_double)
        self.assertEqual(V, [V_Petersen_double])

    def test_trivial_decomposition_5(self):
        V, E = tdlib.trivial_decomposition(V_Wagner, E_Wagner)
        self.assertEqual(V, [V_Wagner])

    def test_trivial_decomposition_6(self):
        V, E = tdlib.trivial_decomposition(V_Pappus, E_Pappus)
        self.assertEqual(V, [V_Pappus])

    def test_trivial_decomposition_7(self):
        V, E = tdlib.trivial_decomposition(V_Grid_5_5, E_Grid_5_5)
        self.assertEqual(V, [V_Grid_5_5])

    def test_trivial_decomposition_8(self):
        for i in range(0, 10):
            V, E = randomGNP(20, 0.2)
            V_T, E_T = tdlib.trivial_decomposition(V, E)
            self.assertEqual(V_T, [V])

    def test_trivial_decomposition_9(self):
        V, E = tdlib.trivial_decomposition(["a", "b", "c"], [])
        self.assertEqual(V, [["a", "b", "c"]])


    def test_is_valid_treedecomposition_0(self):
        V, E = tdlib.trivial_decomposition(V_P6, E_P6)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E), True)

    def test_is_valid_treedecomposition_1(self):
        V, E = tdlib.trivial_decomposition(V_P6, E_P6)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V, E, False), False)

    def test_is_valid_treedecomposition_2(self):
        V, E, tw = tdlib.exact_decomposition_cutset(V_Petersen, E_Petersen)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E), True)

    def test_is_valid_treedecomposition_3(self):
        V, E, tw = tdlib.exact_decomposition_cutset(V_Petersen, E_Petersen)
        del V[-1]
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E, False), False)

    def test_is_valid_treedecomposition_4(self):
        V, E, tw = tdlib.exact_decomposition_cutset(V_Wagner, E_Wagner)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V, E), True)

    def test_is_valid_treedecomposition_5(self):
        V, E, tw = tdlib.exact_decomposition_cutset(V_Wagner, E_Wagner)
        del V[0]
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V, E, False), False)


if __name__ == '__main__':
    unittest.main()
