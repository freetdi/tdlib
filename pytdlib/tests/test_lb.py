import base
import tdlib
import unittest

from graphs import *

class TestTdLib(unittest.TestCase):


##############################################################
############ LOWER BOUNDS ####################################

    def test_lower_bounds_0a(self):
        V, E = cornercases[0]
        G = Graph(V, E)
        lb = tdlib.lower_bound(G, "deltaC_min_d")
        self.assertEqual(lb, -1)
        lb = tdlib.lower_bound(G, "deltaC_max_d")
        self.assertEqual(lb, -1)
        lb = tdlib.lower_bound(G, "deltaC_least_c")
        self.assertEqual(lb, -1)
        lb = tdlib.lower_bound(G, "LBN_deltaC")
        self.assertEqual(lb, -1)
        lb = tdlib.lower_bound(G, "LBNC_deltaC")
        self.assertEqual(lb, -1)
        lb = tdlib.lower_bound(G, "LBP_deltaC")
        self.assertEqual(lb, -1)
        lb = tdlib.lower_bound(G, "LBPC_deltaC")
        self.assertEqual(lb, -1)

    def test_lower_bounds_0b(self):
        V, E = cornercases[1]
        G = Graph(V, E)
        lb = tdlib.lower_bound(G, "deltaC_min_d")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(G, "deltaC_max_d")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(G, "deltaC_least_c")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(G, "LBN_deltaC")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(G, "LBNC_deltaC")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(G, "LBP_deltaC")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(G, "LBPC_deltaC")
        self.assertEqual(lb, 0)

    def test_lower_bounds_0c(self):
        V, E = cornercases[2]
        G = Graph(V, E)
        lb = tdlib.lower_bound(G, "deltaC_min_d")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(G, "deltaC_max_d")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(G, "deltaC_least_c")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(G, "LBN_deltaC")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(G, "LBNC_deltaC")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(G, "LBP_deltaC")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(G, "LBPC_deltaC")
        self.assertEqual(lb, 0)

    def test_lower_bounds_0d(self):
        V, E = cornercases[3]
        G = Graph(V, E)
        lb = tdlib.lower_bound(G, "deltaC_min_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "deltaC_max_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "deltaC_least_c")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBN_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBNC_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBP_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBPC_deltaC")
        self.assertEqual(lb, 4)

    def test_lower_bounds_1(self):
        G = Graph(V_P6, E_P6)
        lb = tdlib.lower_bound(G, "deltaC_min_d")
        self.assertEqual(lb, 1)
        lb = tdlib.lower_bound(G, "deltaC_max_d")
        self.assertEqual(lb, 1)
        lb = tdlib.lower_bound(G, "deltaC_least_c")
        self.assertEqual(lb, 1)
        lb = tdlib.lower_bound(G, "LBN_deltaC")
        self.assertEqual(lb, 1)
        lb = tdlib.lower_bound(G, "LBNC_deltaC")
        self.assertEqual(lb, 1)
        lb = tdlib.lower_bound(G, "LBP_deltaC")
        self.assertEqual(lb, 1)
        lb = tdlib.lower_bound(G, "LBPC_deltaC")
        self.assertEqual(lb, 1)

    def test_lower_bounds_2(self):
        G = Graph(V_K5, E_K5)
        lb = tdlib.lower_bound(G, "deltaC_min_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "deltaC_max_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "deltaC_least_c")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBN_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBNC_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBP_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBPC_deltaC")
        self.assertEqual(lb, 4)

    def test_lower_bounds_3(self):
        G = Graph(V_Petersen, E_Petersen)
        lb = tdlib.lower_bound(G, "deltaC_min_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "deltaC_max_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(G, "deltaC_least_c")
        self.assertEqual(lb, 4) # ??
        lb = tdlib.lower_bound(G, "LBN_deltaC")
        self.assertEqual(lb, 4) # ??
        lb = tdlib.lower_bound(G, "LBNC_deltaC")
        self.assertEqual(lb, 4) # ??
        lb = tdlib.lower_bound(G, "LBP_deltaC")
        self.assertEqual(lb, 4) # ??
        lb = tdlib.lower_bound(G, "LBPC_deltaC")
        self.assertEqual(lb, 4) # ??

    def test_lower_bounds_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        lb = tdlib.lower_bound(G, "deltaC_min_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "deltaC_max_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(G, "deltaC_least_c")
        self.assertEqual(lb, 4) # ??
        lb = tdlib.lower_bound(G, "LBN_deltaC")
        self.assertEqual(lb, 4) # ??
        lb = tdlib.lower_bound(G, "LBNC_deltaC")
        self.assertEqual(lb, 4) # ??
        lb = tdlib.lower_bound(G, "LBP_deltaC")
        self.assertEqual(lb, 4) # ??
        lb = tdlib.lower_bound(G, "LBPC_deltaC")
        self.assertEqual(lb, 4) # ??

    def test_lower_bounds_5(self):
        G = Graph(V_Wagner, E_Wagner)
        lb = tdlib.lower_bound(G, "deltaC_min_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(G, "deltaC_max_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(G, "deltaC_least_c")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(G, "LBN_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(G, "LBNC_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(G, "LBP_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(G, "LBPC_deltaC")
        self.assertEqual(lb, 3)

    def test_lower_bounds_6(self):
        G = Graph(V_Pappus, E_Pappus)
        lb = tdlib.lower_bound(G, "deltaC_min_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "deltaC_max_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(G, "deltaC_least_c")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBN_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBNC_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBP_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBPC_deltaC")
        self.assertEqual(lb, 4)

    def test_lower_bounds_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        lb = tdlib.lower_bound(G, "deltaC_min_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "deltaC_max_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(G, "deltaC_least_c")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBN_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBNC_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBP_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(G, "LBPC_deltaC")
        self.assertEqual(lb, 4)

    def test_lower_bounds_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                G = Graph(V, E)
                lb = tdlib.lower_bound(G, "deltaC_min_d")
                lb = tdlib.lower_bound(G, "deltaC_max_d")
                lb = tdlib.lower_bound(G, "deltaC_least_c")
                lb = tdlib.lower_bound(G, "LBN_deltaC")
                lb = tdlib.lower_bound(G, "LBNC_deltaC")
                lb = tdlib.lower_bound(G, "LBP_deltaC")
                lb = tdlib.lower_bound(G, "LBPC_deltaC")

if __name__ == '__main__':
    unittest.main()

