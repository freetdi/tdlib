import base
import sys
import tdlib
import unittest
import util

import Dimacs
from graphs import *

PREFIX = "Dimacs"
COUNT = 81

#don't confuse python unittest
sys.argv=sys.argv[:1]

class TestTdLib(unittest.TestCase):
    def test_0(self):
        print("---cornercases---")
        for V, E in cornercases:
            G = Graph(V, E)
            tdlib.generic_elimination_search1(G)

    def test_1(self):
        print("---P6---")
        G = Graph(V_P6, E_P6)
        tdlib.generic_elimination_search1(G)

    def test_3(self):
        print("---Petersen---")
        G = Graph(V_Petersen, E_Petersen)
        tdlib.generic_elimination_search1(G)

    def test_4(self):
        print("---Petersen_double---")
        G = Graph(V_Petersen_double, E_Petersen_double)
        tdlib.generic_elimination_search1(G)

    def test_5(self):
        print("---Wagner---")
        G = Graph(V_Wagner, E_Wagner)
        tdlib.generic_elimination_search1(G)

    """ takes too long on CFG1
    def test_6(self):
        G = Graph(V_Pappus, E_Pappus)
        tdlib.generic_elimination_search(G)
    """

    """ takes too long on CFG_t1
    def test_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        tdlib.generic_elimination_search(G)
    """

if __name__ == '__main__':
    unittest.main()
