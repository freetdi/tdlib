import tdlib
import unittest

from graphs import *
import CFGs
import dimacs

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
    #indirect test
    def test_CFGs_LB1(self):
        for i in range(0, COUNT+1):
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            tdlib.lower_bound(G, "deltaC_min_d")

    #indirect test
    def test_CFGs_LB2(self):
        for i in range(0, COUNT+1):
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            tdlib.lower_bound(G, "deltaC_max_d")

    #indirect test
    def test_CFGs_LB3(self):
        for i in range(0, COUNT+1):
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            tdlib.lower_bound(G, "deltaC_least_c")

    #indirect test
    def test_CFGs_LB4(self):
        for i in range(0, COUNT+1):
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            tdlib.lower_bound(G, "LBN_deltaC")

    #indirect test
    def test_CFGs_LB5(self):
        for i in range(0, COUNT+1):
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            tdlib.lower_bound(G, "LBNC_deltaC")

    #indirect test
    def test_CFGs_LB6(self):
        for i in range(0, COUNT+1):
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            tdlib.lower_bound(G, "LBP_deltaC")

    #indirect test
    def test_CFGs_LB7(self):
        for i in range(0, COUNT+1):
            G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)))
            tdlib.lower_bound(G, "LBPC_deltaC")

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
