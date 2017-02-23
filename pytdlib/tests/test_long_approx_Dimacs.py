import base
import tdlib
import unittest
import sys
from graphs import *
import Dimacs

PREFIX = "Dimacs"
COUNT = 81

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
    def test_MD(self):
        #print("---MD---")
        for i in range(0, COUNT+1):
            #name = eval(PREFIX+".name_"+str(i))
            #print("     " + name)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_boost_MD(self):
        #print("---boost::MD---")
        for i in range(0, COUNT+1):
            #name = eval(PREFIX+".name_"+str(i))
            #print("     " + name)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.boost_minDegree_decomp(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_PP(self):
       #print("---PP---")
       for i in range(0, COUNT+1):
            #name = eval(PREFIX+".name_"+str(i))
            #print("     " + name)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            G_, B, lb = tdlib.preprocessing(G)
            if G.vertices() is []:
                self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_PP_MD(self):
        #print("---PP+MD---")
        for i in range(0, COUNT+1):
            #name = eval(PREFIX+".name_"+str(i))
            #print("     " + name)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.PP_MD(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_minimalChordal_trivial(self):
        #print("---minimalChordal (trivial)---")
        for i in range(0, COUNT+1):
            #name = eval(PREFIX+".name_"+str(i))
            #print("     " + name)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T1, w1 = tdlib.trivial_decomposition(G)
            T2, w2 = tdlib.minimalChordal_decomp(G, T1)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T2), True)

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
