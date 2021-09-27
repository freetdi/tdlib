import base
import tdlib
import unittest
import sys

COUNT = 82

if(len(sys.argv)==2 and sys.argv[1]=="short"):
    COUNT = 6
elif(len(sys.argv)<2 or sys.argv[1]!="long"):
    sys.exit(77)

from graphs import *
import Dimacs

#don't confuse python unittest
sys.argv=sys.argv[:1]

PREFIX = "Dimacs"

class TestTdLib_packages(unittest.TestCase):
    def test_FI(self):
        print("---FI---")
        for i in range(0, COUNT):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.fillIn_decomp(G)
            print(i, w)
            assert(tdlib.is_valid_treedecomposition(G, T))

    def test_PP_FI(self):
        print("---PP+FI---")
        for i in range(0, COUNT):
            print(i, end=" ")
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.PP_FI(G)
            print(i, w)
            assert(tdlib.is_valid_treedecomposition(G, T))

   # too long
    # def test_PP_FI_TM(self):
    #     print("---PP+FI+TM---")
    #     for i in range(0, COUNT):
    #         base.print_graph_name(PREFIX, i)
    #         G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
    #         T, w = tdlib.PP_FI_TM(G)
    #         print(i, w)
    #         assert(tdlib.is_valid_treedecomposition(G, T))

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
