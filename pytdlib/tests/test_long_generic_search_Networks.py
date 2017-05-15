import base
import sys
import tdlib
import unittest
import util

if(len(sys.argv)<2 or sys.argv[1]!="long"):
    sys.exit(77)

import Networks
from graphs import *

PREFIX_NETWORKS = "Networks"
COUNT_NETWORKS = 21

MAX_NODES = 1000000
MAX_ORDERINGS = 100

#don't confuse python unittest
sys.argv=sys.argv[:1]

class TestTdLib(unittest.TestCase):
    def test_long1(self):
        print("FILL config")
        for i in range(0, COUNT_NETWORKS+1):
            #if base.skip(PREFIX_NETWORKS, i, lambda x,y: x > 100 or y > 2000):
            #    continue
#            if i != 3:
#                continue

            base.print_graph_name(PREFIX_NETWORKS, i)

            G = Graph(eval(PREFIX_NETWORKS+".V_"+str(i)), eval(PREFIX_NETWORKS+".E_"+str(i)))

            tdlib.generic_elimination_search_p17_jumper(G, MAX_NODES, MAX_ORDERINGS)
            print("")

if __name__ == '__main__':
    unittest.main()
