import base
import sys
import tdlib
import unittest

if(len(sys.argv)<2 or sys.argv[1]!="long"):
    sys.exit(77)

from graphs import *
import CFGs

#don't confuse python unittest
sys.argv=sys.argv[:1]

PREFIX = "CFGs"
COUNT = 1816

class TestTdLib_packages(unittest.TestCase):
    #TODO: validation (is_clique, is_IS, is_VC,..) in tdlib?
    def test_max_clique(self):
        print("---maxClique--")
        for i in range(0, COUNT+1):
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            base.print_graph_name(PREFIX, i)
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.max_clique_with_treedecomposition(G, T)

    def test_max_independent_set(self):
        print("---maxIndependentSet---")
        for i in range(0, COUNT+1):
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            base.print_graph_name(PREFIX, i)
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.max_independent_set_with_treedecomposition(G, T)

    def test_min_vertex_cover(self):
        print("---minVertexCover--")
        for i in range(0, COUNT+1):
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            base.print_graph_name(PREFIX, i)
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.min_vertex_cover_with_treedecomposition(G, T)

    def test_min_dominating_set(self):
        print("---minDominatingSet--")
        for i in range(0, COUNT+1):
            if i == 999: #huge graph
                continue;

            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            base.print_graph_name(PREFIX, i)
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.min_dominating_set_with_treedecomposition(G, T)

    def test_min_coloring(self):
        print("---minColoring--")
        for i in range(0, COUNT+1):
            if i == 999: #huge graph
                continue

            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            base.print_graph_name(PREFIX, i)
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.min_coloring_with_treedecomposition(G, T)

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
