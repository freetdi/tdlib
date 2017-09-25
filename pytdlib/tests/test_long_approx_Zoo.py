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
    def test_MD(self):
        print("---MD---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.minDegree_decomp(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_boost_MD(self):
	print("---boost_MD---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.boost_minDegree_decomp(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_FI(self):
	print("---FI---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.fillIn_decomp(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_PP(self):
	print("---PP---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            G_, B, lb = tdlib.preprocessing(G)
            if G.vertices() is []:
                self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_PP_MD(self):
	print("---PP_MD---")
        for i in range(0, COUNT+1):
	    base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.PP_MD(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_PP_FI(self):
	print("---PP_FI---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.PP_FI(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_PP_FI_TM(self):
	print("---PP_FI_TM---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.PP_FI_TM(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_MSVS_trivial(self):
	print("---MSVS_trivial---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T1, w1 = tdlib.trivial_decomposition(G)
            T2, w2 = tdlib.MSVS(G, T1)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T2), True)

    """
    def test_seperator_algorithm(self):
	print("---seperator_algorithm---")
        for i in range(0, COUNT+1):
            if base.skip(PREFIX, i, lambda x,y: x > 100 or y > 300):
                continue
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T, w = tdlib.seperator_algorithm(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
    """

    def test_conversion(self):
	print("---conversion---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            T1, w1 = tdlib.minDegree_decomp(G)
            O = tdlib.treedec_to_ordering(T1)
            T2, w2 = tdlib.ordering_to_treedec(G, O)
            self.assertEqual(w1, w2)

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
