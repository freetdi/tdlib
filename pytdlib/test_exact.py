import tdlib
import unittest

from graphs import *

class TestTdLib_exact(unittest.TestCase):
    def test_exact_decomposition_cutset_0(self):
        for V, E in cornercases:
            V_T, E_T, w = tdlib.exact_decomposition_cutset(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    def test_exact_decomposition_cutset_1(self):
        V, E, w = tdlib.exact_decomposition_cutset(V_P6, E_P6)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E), True)
        self.assertEqual(w, 1)

    def test_exact_decomposition_cutset_2(self):
        V, E, w = tdlib.exact_decomposition_cutset(V_K5, E_K5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V, E), True)
        self.assertEqual(w, 4)

    def test_exact_decomposition_cutset_3(self):
        V, E, w = tdlib.exact_decomposition_cutset(V_Petersen, E_Petersen)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E), True)
        self.assertEqual(w, 4)

    def test_exact_decomposition_cutset_4(self):
        V, E, w = tdlib.exact_decomposition_cutset(V_Petersen_double, E_Petersen_double)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen_double, E_Petersen_double, V, E), True)
        self.assertEqual(w, 4)

    def test_exact_decomposition_cutset_5(self):
        V, E, w = tdlib.exact_decomposition_cutset(V_Wagner, E_Wagner)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V, E), True)
        self.assertEqual(w, 4)

    """ takes too much time
    def test_exact_decomposition_cutset_6(self):
        V, E, w = tdlib.exact_decomposition_cutset(V_Pappus, E_Pappus)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V, E), True)
        self.assertEqual(w, 6)
    """

    def test_exact_decomposition_cutset_7(self):
        V, E, w = tdlib.exact_decomposition_cutset(V_Grid_5_5, E_Grid_5_5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Grid_5_5, E_Grid_5_5, V, E), True)
        self.assertEqual(w, 5)

    def test_exact_decomposition_cutset_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                V_T, E_T, w = tdlib.exact_decomposition_cutset(V, E)
                self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    def test_exact_decomposition_cutset_9(self):
        status = True
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.3)
                Q, R, w2 = tdlib.PP_FI_TM(V, E)
                isleq = tdlib.exact_decomposition_cutset_decision(V, E, w2)
                if(not isleq):
                    print("error [validate width], graph: " + str(V) + ", " + str(E))
                    print("width_PP_FI_TM: " + str(w2))
                    status = False

        self.assertEqual(status, True)

    """ segfault
    def test_exact_decomposition_cutset_0(self):
        for V, E in cornercases:
            V_T, E_T, w = tdlib.exact_decomposition_dynamic(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    def test_exact_decomposition_dynamic_1(self):
        V, E, tw = tdlib.exact_decomposition_dynamic(V_P6, E_P6)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E), True)
        self.assertEqual(tw, 1)

    def test_exact_decomposition_dynamic_2(self):
        V, E, tw = tdlib.exact_decomposition_dynamic(V_K5, E_K5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V, E), True)
        self.assertEqual(tw, 4)

if __name__ == '__main__':
    unittest.main()

# vim:ts=4:sw=4:et
