import base
import sys
import tdlib
import unittest

from graphs import *

#don't confuse python unittest
sys.argv=sys.argv[:1]

class TestTdLib(unittest.TestCase):

    def test_lower_bounds_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        lb = tdlib.lower_bound(G, "deltaC_least_c")
        self.assertEqual(lb, 4)
        print("pass")

if __name__ == '__main__':
    unittest.main()

