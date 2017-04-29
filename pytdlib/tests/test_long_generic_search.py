import base
import sys
import tdlib
import unittest
import util

if(len(sys.argv)<2 or sys.argv[1]!="long"):
    sys.exit(77)

import Dimacs
from graphs import *

PREFIX_DIMACS = "Dimacs"
COUNT_DIMACS = 81

PREFIX_MAXSAT = "Maxsat_small"
COUNT_MAXSAT = 120

MAX_NODES = 5000
MAX_ORDERINGS = 100

#don't confuse python unittest
sys.argv=sys.argv[:1]

class TestTdLib(unittest.TestCase):
    #"""
    def test_Dimacs37(self):
        G = Graph(Dimacs.V_37, Dimacs.E_37)

        T, w = tdlib.minDegree_decomp(G)
        print("MD_width: " + str(w+1))

        T, w = tdlib.fillIn_decomp(G)
        print("FI_width: " + str(w+1))

        tdlib.generic_elimination_search4(G, MAX_NODES, MAX_ORDERINGS)
    #"""

    """
    def test_long1(self):
        for i in range(0, COUNT_DIMACS+1):
            if base.skip(PREFIX_DIMACS, i, lambda x,y: x > 100 or y > 2000):
                continue

            base.print_graph_name(PREFIX_DIMACS, i)
            
            G = Graph(eval(PREFIX_DIMACS+".V_"+str(i)), eval(PREFIX_DIMACS+".E_"+str(i)))

            tdlib.generic_elimination_search4(G, MAX_NODES, MAX_ORDERINGS)
            print("")
    """

if __name__ == '__main__':
    unittest.main()
