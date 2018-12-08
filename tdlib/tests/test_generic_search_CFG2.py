import base
import sys
import tdlib
import unittest
import util

import Dimacs
from graphs import *

PREFIX = "Dimacs"
COUNT = 81

MAX_NODES = 2500
MAX_ORDERINGS = 1000

#don't confuse python unittest
sys.argv=sys.argv[:1]

class TestTdLib(unittest.TestCase):
    def test_0(self):
        print("---cornercases---")
        for V, E in cornercases:
            G = Graph(V, E)
            tdlib.generic_elimination_search2(G, MAX_NODES, MAX_ORDERINGS)

    def test_1(self):
        print("---P6---")
        G = Graph(V_P6, E_P6)
        tdlib.generic_elimination_search2(G, MAX_NODES, MAX_ORDERINGS)

    def test_3(self):
        print("---Petersen---")
        G = Graph(V_Petersen, E_Petersen)
        tdlib.generic_elimination_search2(G, MAX_NODES, MAX_ORDERINGS)

    def test_4(self):
        print("---Petersen_double---")
        G = Graph(V_Petersen_double, E_Petersen_double)
        tdlib.generic_elimination_search2(G, MAX_NODES, MAX_ORDERINGS)

    def test_5(self):
        print("---Wagner---")
        G = Graph(V_Wagner, E_Wagner)
        tdlib.generic_elimination_search2(G, MAX_NODES, MAX_ORDERINGS)

    def test_6(self):
        print("---Pappus---")
        G = Graph(V_Pappus, E_Pappus)
        tdlib.generic_elimination_search2(G, MAX_NODES, MAX_ORDERINGS)

    def test_7(self):
        print("---Grid-5-5---")
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        tdlib.generic_elimination_search2(G, MAX_NODES, MAX_ORDERINGS)

    def test_8(self):
        print("---Dimacs7---")
        G = Graph(Dimacs.V_7, Dimacs.E_7)

        T, w = tdlib.minDegree_decomp(G)
        print("MD_width: " + str(w))

        tdlib.generic_elimination_search2(G, MAX_NODES, MAX_ORDERINGS)

    def test_9(self):
        print("---Dimacs58---")
        G = Graph(Dimacs.V_58, Dimacs.E_58)

        T, w = tdlib.minDegree_decomp(G)
        print("MD_width: " + str(w))

        tdlib.generic_elimination_search2(G, MAX_NODES, MAX_ORDERINGS)

if __name__ == '__main__':
    unittest.main()
