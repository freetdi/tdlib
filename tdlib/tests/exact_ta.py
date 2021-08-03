import base
import sys
import tdlib
import unittest

from graphs import *

#don't confuse python unittest
sys.argv=sys.argv[:1]

class TestTdLib_exact(unittest.TestCase):
    def test_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            T, w = tdlib.exact_decomposition_ex17(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_1(self):
        Gs = [ Graph(V_P6, E_P6), Graph(V_Petersen, E_Petersen),
            Graph(V_Petersen_double, E_Petersen_double),
            Graph(V_Wagner, E_Wagner), Graph(V_Pappus, E_Pappus),
            Graph(V_Grid_5_5, E_Grid_5_5) ]

        for G in Gs:
            T, w = tdlib.exact_decomposition_ex17(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
            print("TA", w)

    def test_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                G = Graph(V, E)
                T, w = tdlib.exact_decomposition_ex17(G)
                self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    """ not yet
    def test_9(self):
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
    """

if __name__ == '__main__':
    unittest.main()

# vim:ts=4:sw=4:et
