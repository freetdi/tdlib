import base
import tdlib
import unittest
import sys

if(len(sys.argv)<2 or sys.argv[1]!="long"):
    sys.exit(77)

from graphs import *
import Dimacs

#don't confuse python unittest
sys.argv=sys.argv[:1]

PREFIX = "Dimacs"
COUNT = 81

class TestTdLib_packages(unittest.TestCase):
    def test_FI(self):
        print("---FI---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.fillIn_decomp(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
            break
    """
    def test_PP_FI(self):
        print("---PP+FI---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.PP_FI(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_PP_FI_TM(self):
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.PP_FI_TM(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
    """

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
