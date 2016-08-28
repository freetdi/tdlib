import tdlib
import unittest

from graphs import *

class TestTdLib_exact(unittest.TestCase):
    def test_exact_decomposition_cutset(self):
        N, M, tw = tdlib.exact_decomposition_cutset(V_Grid_5_5, E_Grid_5_5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Grid_5_5, E_Grid_5_5, N, M), True)
        self.assertEqual(tw, 5)

        V, E, tw = tdlib.exact_decomposition_cutset(V_P6, E_P6)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E), True)
        self.assertEqual(tw, 1)

        V, E, tw = tdlib.exact_decomposition_cutset(V_K5, E_K5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V, E), True)
        self.assertEqual(tw, 4)

        V, E, tw = tdlib.exact_decomposition_cutset(V_Petersen, E_Petersen)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E), True)
        self.assertEqual(tw, 4)

        V, E, tw = tdlib.exact_decomposition_cutset(V_Wagner, E_Wagner)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V, E), True)
        self.assertEqual(tw, 4)

    def test_exact_decomposition_dynamic(self):
        V, E, tw = tdlib.exact_decomposition_dynamic(V_P6, E_P6)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E), True)
        self.assertEqual(tw, 1)

        V, E, tw = tdlib.exact_decomposition_dynamic(V_K5, E_K5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V, E), True)
        self.assertEqual(tw, 4)

    def test_random_valid_treedecomposition(self):
        status_list = list()
        correct_status = list()
        for n in range(1, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.3)
                N, M, tw = tdlib.exact_decomposition_cutset(V, E)
                status = tdlib.is_valid_treedecomposition(V, E, N, M, message=True)
                status_list.append(status)
                correct_status.append(0)
                if(status < 0):
                    print("error in [validate decompositions], graph: " + str(V) + ", " + str(E))

    def test_random_validate_width(self):
        status = True
        for n in range(0, 13):
            for i in range(0, 100):
                V, E = randomGNP(n, 0.3)
                Q, R, w2 = tdlib.PP_FI_TM(V, E)
                isleq = tdlib.exact_decomposition_cutset_decision(V, E, w2)
                if(not isleq):
                    print("error [validate width], graph: " + str(V) + ", " + str(E))
                    print("width_PP_FI_TM: " + str(w2))
                    status = False

        self.assertEqual(status, True)

if __name__ == '__main__':
    unittest.main()

# vim:ts=4:sw=4:et
