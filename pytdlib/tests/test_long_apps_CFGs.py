import base
import tdlib
import unittest
import sys
from graphs import *
import CFGs

PREFIX = "CFGs"
COUNT = 1816

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
        for i in range(0, COUNT+1):
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.max_clique_with_treedecomposition(G, T)

    def test_max_independent_set(self):
        for i in range(0, COUNT+1):
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.max_independent_set_with_treedecomposition(G, T)

    def test_min_vertex_cover(self):
        for i in range(0, COUNT+1):
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.min_vertex_cover_with_treedecomposition(G, T)

    def test_min_dominating_set(self):
        for i in range(0, COUNT+1):
            if i == 999: #huge graph
                continue;

            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.min_dominating_set_with_treedecomposition(G, T)

    def test_min_coloring(self):
        for i in range(0, COUNT+1):
            if i == 999: #huge graph
                continue

            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            S = tdlib.min_coloring_with_treedecomposition(G, T)

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
