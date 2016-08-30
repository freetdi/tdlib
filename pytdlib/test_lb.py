import tdlib
import unittest

from graphs import *

class TestTdLib(unittest.TestCase):


##############################################################
############ LOWER BOUNDS ####################################

    def test_lower_bounds_0a(self):
        V, E = cornercases[0]
        lb = tdlib.lower_bound(V, E, "deltaC_min_d")
        self.assertEqual(lb, -1)
        lb = tdlib.lower_bound(V, E, "deltaC_max_d")
        self.assertEqual(lb, -1)
        lb = tdlib.lower_bound(V, E, "deltaC_least_c")
        self.assertEqual(lb, -1)
        lb = tdlib.lower_bound(V, E, "LBN_deltaC")
        self.assertEqual(lb, -1)
        lb = tdlib.lower_bound(V, E, "LBNC_deltaC")
        self.assertEqual(lb, -1)
        lb = tdlib.lower_bound(V, E, "LBP_deltaC")
        self.assertEqual(lb, -1)
        lb = tdlib.lower_bound(V, E, "LBPC_deltaC")
        self.assertEqual(lb, -1)

    def test_lower_bounds_0b(self):
        V, E = cornercases[1]
        lb = tdlib.lower_bound(V, E, "deltaC_min_d")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(V, E, "deltaC_max_d")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(V, E, "deltaC_least_c")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(V, E, "LBN_deltaC")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(V, E, "LBNC_deltaC")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(V, E, "LBP_deltaC")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(V, E, "LBPC_deltaC")
        self.assertEqual(lb, 0)

    def test_lower_bounds_0c(self):
        V, E = cornercases[2]
        lb = tdlib.lower_bound(V, E, "deltaC_min_d")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(V, E, "deltaC_max_d")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(V, E, "deltaC_least_c")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(V, E, "LBN_deltaC")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(V, E, "LBNC_deltaC")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(V, E, "LBP_deltaC")
        self.assertEqual(lb, 0)
        lb = tdlib.lower_bound(V, E, "LBPC_deltaC")
        self.assertEqual(lb, 0)

    def test_lower_bounds_0d(self):
        V, E = cornercases[3]
        lb = tdlib.lower_bound(V, E, "deltaC_min_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V, E, "deltaC_max_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V, E, "deltaC_least_c")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V, E, "LBN_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V, E, "LBNC_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V, E, "LBP_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V, E, "LBPC_deltaC")
        self.assertEqual(lb, 4)

    def test_lower_bounds_1(self):
        lb = tdlib.lower_bound(V_P6, E_P6, "deltaC_min_d")
        self.assertEqual(lb, 1)
        lb = tdlib.lower_bound(V_P6, E_P6, "deltaC_max_d")
        self.assertEqual(lb, 1)
        lb = tdlib.lower_bound(V_P6, E_P6, "deltaC_least_c")
        self.assertEqual(lb, 1)
        lb = tdlib.lower_bound(V_P6, E_P6, "LBN_deltaC")
        self.assertEqual(lb, 1)
        lb = tdlib.lower_bound(V_P6, E_P6, "LBNC_deltaC")
        self.assertEqual(lb, 1)
        lb = tdlib.lower_bound(V_P6, E_P6, "LBP_deltaC")
        self.assertEqual(lb, 1)
        lb = tdlib.lower_bound(V_P6, E_P6, "LBPC_deltaC")
        self.assertEqual(lb, 1)

    def test_lower_bounds_2(self):
        lb = tdlib.lower_bound(V_K5, E_K5, "deltaC_min_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_K5, E_K5, "deltaC_max_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_K5, E_K5, "deltaC_least_c")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_K5, E_K5, "LBN_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_K5, E_K5, "LBNC_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_K5, E_K5, "LBP_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_K5, E_K5, "LBPC_deltaC")
        self.assertEqual(lb, 4)

    def test_lower_bounds_3(self):
        lb = tdlib.lower_bound(V_Petersen, E_Petersen, "deltaC_min_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Petersen, E_Petersen, "deltaC_max_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Petersen, E_Petersen, "deltaC_least_c")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Petersen, E_Petersen, "LBN_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Petersen, E_Petersen, "LBNC_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Petersen, E_Petersen, "LBP_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Petersen, E_Petersen, "LBPC_deltaC")
        self.assertEqual(lb, 3)

    def test_lower_bounds_4(self):
        lb = tdlib.lower_bound(V_Petersen_double, E_Petersen_double, "deltaC_min_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Petersen_double, E_Petersen_double, "deltaC_max_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Petersen_double, E_Petersen_double, "deltaC_least_c")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Petersen_double, E_Petersen_double, "LBN_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Petersen_double, E_Petersen_double, "LBNC_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Petersen_double, E_Petersen_double, "LBP_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Petersen_double, E_Petersen_double, "LBPC_deltaC")
        self.assertEqual(lb, 3)

    def test_lower_bounds_5(self):
        lb = tdlib.lower_bound(V_Wagner, E_Wagner, "deltaC_min_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Wagner, E_Wagner, "deltaC_max_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Wagner, E_Wagner, "deltaC_least_c")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Wagner, E_Wagner, "LBN_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Wagner, E_Wagner, "LBNC_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Wagner, E_Wagner, "LBP_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Wagner, E_Wagner, "LBPC_deltaC")
        self.assertEqual(lb, 3)

    def test_lower_bounds_6(self):
        lb = tdlib.lower_bound(V_Pappus, E_Pappus, "deltaC_min_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Pappus, E_Pappus, "deltaC_max_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Pappus, E_Pappus, "deltaC_least_c")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Pappus, E_Pappus, "LBN_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Pappus, E_Pappus, "LBNC_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Pappus, E_Pappus, "LBP_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Pappus, E_Pappus, "LBPC_deltaC")
        self.assertEqual(lb, 4)

    def test_lower_bounds_7(self):
        lb = tdlib.lower_bound(V_Grid_5_5, E_Grid_5_5, "deltaC_min_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Grid_5_5, E_Grid_5_5, "deltaC_max_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Grid_5_5, E_Grid_5_5, "deltaC_least_c")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Grid_5_5, E_Grid_5_5, "LBN_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Grid_5_5, E_Grid_5_5, "LBNC_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Grid_5_5, E_Grid_5_5, "LBP_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Grid_5_5, E_Grid_5_5, "LBPC_deltaC")
        self.assertEqual(lb, 4)

    def test_lower_bounds_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                lb = tdlib.lower_bound(V, E, "deltaC_min_d")
                lb = tdlib.lower_bound(V, E, "deltaC_max_d")
                lb = tdlib.lower_bound(V, E, "deltaC_least_c")
                lb = tdlib.lower_bound(V, E, "LBN_deltaC")
                lb = tdlib.lower_bound(V, E, "LBNC_deltaC")
                lb = tdlib.lower_bound(V, E, "LBP_deltaC")
                lb = tdlib.lower_bound(V, E, "LBPC_deltaC")

if __name__ == '__main__':
    unittest.main()

