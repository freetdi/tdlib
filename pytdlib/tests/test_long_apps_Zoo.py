import base
import sys
import tdlib
import unittest

if(len(sys.argv)<2 or sys.argv[1]!="long"):
    sys.exit(77)

from graphs import *
import Zoo

#don't confuse python unittest
sys.argv=sys.argv[:1]

PREFIX = "Zoo"
COUNT = 149

class TestTdLib_packages(unittest.TestCase):
    #TODO: validation (is_clique, is_IS, is_VC,..) in tdlib?
    """
    def test_max_clique(self):
        print("---max_clique---")
        for i in range(0, COUNT+1):
            if i == 91:
                 continue

            if base.skip(PREFIX, i, lambda x,y: y > 2000):
                continue

            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))

            T, w = tdlib.minDegree_decomp(G)

            print("... width " + str(w))

            if w > 32:
                continue

            S = tdlib.max_clique_with_treedecomposition(G, T)

    def test_max_independent_set(self):
        print("---independent_set---")
        for i in range(0, COUNT+1):
            if i == 91:
                 continue

            if base.skip(PREFIX, i, lambda x,y: y > 2000):
                continue

            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            if w > 17:
                print("...tw > 17, skipping...")
                continue

            S = tdlib.max_independent_set_with_treedecomposition(G, T)
    """
    def test_min_vertex_cover(self):
        print("---vertex_cover---")
        for i in range(0, COUNT+1):
            if i == 91:
                 continue

            if base.skip(PREFIX, i, lambda x,y: y > 10000):
                continue

            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            if w > 25:
                print("...tw > 25, skipping...")
                continue

            S = tdlib.min_vertex_cover_with_treedecomposition(G, T)
            print(str(S))
    """
    def test_min_dominating_set(self):
        print("---dominating_set---")
        for i in range(0, COUNT+1):
            if i == 91:
                 continue

            if base.skip(PREFIX, i, lambda x,y: y > 2000):
                continue

            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            if w > 8:
                print("...tw > 8, skipping...")
                continue

            S = tdlib.min_dominating_set_with_treedecomposition(G, T)

    def test_min_coloring(self):
        print("---coloring---")
        for i in range(0, COUNT+1):
            if i == 97 or i == 91:
                 continue

            if base.skip(PREFIX, i, lambda x,y: y > 2000):
                continue

            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            if w > 8:
                print("...tw > 8, skipping...")
                continue

            S = tdlib.min_coloring_with_treedecomposition(G, T)
    """


if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
