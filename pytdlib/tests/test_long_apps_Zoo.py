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

def dump_td_as_dot(V_T, E_T, outname):
    fout = open(outname, 'w')
    fout.write("digraph G{\n")
    for i in range(0, len(V_T)):
        fout.write(str(i) + "[label=\"")
        for v in V_T[i]:
            fout.write(str(v) + " ")
        fout.write("\"];\n")
    for i in range(0, len(E_T)-1, 2):
        fout.write(str(E_T[i]) + " -> " + str(E_T[i+1]) + ";\n")
    fout.write("}\n")
    fout.close()

class TestTdLib_packages(unittest.TestCase):
    #TODO: validation (is_clique, is_IS, is_VC,..) in tdlib?
    def test_max_clique(self):
        print("---max_clique---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.max_clique_with_treedecomposition(G, T)

    def test_max_independent_set(self):
        print("---independent_set---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.max_independent_set_with_treedecomposition(G, T)

    def test_min_vertex_cover(self):
        print("---vertex_cover---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.min_vertex_cover_with_treedecomposition(G, T)

    def test_min_dominating_set(self):
        print("---dominating_set---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.min_dominating_set_with_treedecomposition(G, T)

    def test_min_coloring(self):
        print("---coloring---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.min_coloring_with_treedecomposition(G, T)

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
