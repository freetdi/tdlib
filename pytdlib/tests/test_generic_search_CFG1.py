import base
import sys
import tdlib
import unittest
import util

from graphs import *

MAX_NODES1 = 1000000
MAX_NODES2 = 10000
MAX_ORDERINGS = 1000

#don't confuse python unittest
sys.argv=sys.argv[:1]

class TestTdLib(unittest.TestCase):
    def test_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            tdlib.generic_elimination_search1(G, MAX_NODES1, MAX_ORDERINGS)

    def test_1(self):
        G = Graph(V_P6, E_P6)
        tdlib.generic_elimination_search1(G, MAX_NODES1, MAX_ORDERINGS)

    def test_3(self):
        G = Graph(V_Petersen, E_Petersen)
        tdlib.generic_elimination_search1(G, MAX_NODES1, MAX_ORDERINGS)

    def test_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        tdlib.generic_elimination_search1(G, MAX_NODES1, MAX_ORDERINGS)

    def test_5(self):
        G = Graph(V_Wagner, E_Wagner)
        tdlib.generic_elimination_search1(G, MAX_NODES1, MAX_ORDERINGS)

    def test_6(self):
        G = Graph(V_Pappus, E_Pappus)
        tdlib.generic_elimination_search1(G, MAX_NODES2, MAX_ORDERINGS)

    def test_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        tdlib.generic_elimination_search1(G, MAX_NODES2, MAX_ORDERINGS)

if __name__ == '__main__':
    unittest.main()
