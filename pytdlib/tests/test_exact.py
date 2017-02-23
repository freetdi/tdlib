import base
import tdlib
import unittest

from graphs import *

class TestTdLib_exact(unittest.TestCase):
    def test_exact_decomposition_cutset_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            T, w = tdlib.exact_decomposition_cutset(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_exact_decomposition_cutset_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.exact_decomposition_cutset(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 1)

    def test_exact_decomposition_cutset_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.exact_decomposition_cutset(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_exact_decomposition_cutset_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.exact_decomposition_cutset(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_exact_decomposition_cutset_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.exact_decomposition_cutset(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_exact_decomposition_cutset_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.exact_decomposition_cutset(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    """ takes too long
    def test_exact_decomposition_cutset_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.exact_decomposition_cutset(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 6)
    """

    def test_exact_decomposition_cutset_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.exact_decomposition_cutset(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 5)

    def test_exact_decomposition_cutset_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                G = Graph(V, E)
                T, w = tdlib.exact_decomposition_cutset(G)
                self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_exact_decomposition_cutset_9(self):
        status = True
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.3)
                G = Graph(V, E)
                T, w2 = tdlib.PP_FI_TM(G)
                isleq = tdlib.exact_decomposition_cutset_decision(G, w2)
                if(not isleq):
                    print("error [validate width], graph: " + str(G.vertices()) + ", " + str(G.edges()))
                    print("width_PP_FI_TM: " + str(w2))
                    N, M, width = tdlib.exact_decomposition_cutset(G)
                    hrgl = tdlib.check_treedec(V, E, N, M, message=True)
                    print("proper width_PP_FI_TM: " + str(width) + " " + str(hrgl))
                    self.assertEqual(hrgl, 0)
                    status = False

        self.assertEqual(status, True)

    def test_exact_decomposition_dynamic_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            T, w = tdlib.exact_decomposition_dynamic(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_exact_decomposition_dynamic_1(self):
        G = Graph(V_P6, E_P6)
        T, tw = tdlib.exact_decomposition_dynamic(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(tw, 1)

    def test_exact_decomposition_dynamic_2(self):
        G = Graph(V_K5, E_K5)
        T, tw = tdlib.exact_decomposition_dynamic(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(tw, 4)

if __name__ == '__main__':
    unittest.main()

# vim:ts=4:sw=4:et
