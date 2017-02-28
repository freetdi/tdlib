import base
import tdlib
import unittest
import sys

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
    #indirect test
    def test_LB1(self):
        print("---deltaC_min_d---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            tdlib.lower_bound(G, "deltaC_min_d")

    #indirect test
    def test_LB2(self):
        print("---deltaC_max_d---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            tdlib.lower_bound(G, "deltaC_max_d")

    #indirect test
    def test_LB3(self):
        print("---deltaC_least_c---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            tdlib.lower_bound(G, "deltaC_least_c")

    #indirect test
    def test_LB4(self):
        print("---LBN_deltaC---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            tdlib.lower_bound(G, "LBN_deltaC")

    #indirect test
    def test_LB5(self):
        print("---LBNC_deltaC---")
        for i in range(0, COUNT+1):
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            tdlib.lower_bound(G, "LBNC_deltaC")

    #indirect test
    def test_LB6(self):
        print("---LBP_deltaC---")
        for i in range(0, COUNT+1):
            if base.skip(PREFIX, i, lambda x,y: x > 1000 or y > 3000):
                continue
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            tdlib.lower_bound(G, "LBP_deltaC")

    #indirect test
    def test_LB7(self):
        print("---LBPC_deltaC---")
        for i in range(0, COUNT+1):
            if base.skip(PREFIX, i, lambda x,y: x > 500 or y > 1000):
                continue
            base.print_graph_name(PREFIX, i)
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            tdlib.lower_bound(G, "LBPC_deltaC")

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
