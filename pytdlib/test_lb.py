import tdlib
import unittest

from graphs import *

class TestTdLib(unittest.TestCase):


##############################################################
############ LOWER BOUNDS ####################################

    def test_lower_bounds(self):
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

if __name__ == '__main__':
    unittest.main()

