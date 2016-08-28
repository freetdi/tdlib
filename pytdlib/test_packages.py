import tdlib
import unittest

from graphs import *
import CFGs
import dimacs

CFGS_count = 1816
DIMACS_count = 81
VERBOSE = True


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
    """ works...
    def test_CFGs_MD(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.minDegree_decomp(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ works...
    def test_CFGs_boost_MD(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.boost_minDegree_decomp(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ works...
    def test_CFGs_FI(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.fillIn_decomp(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ works...
    def test_CFGs_PP(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, B, lb = tdlib.preprocessing(V, E)
            if V_T is []:
                self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ works....
    def test_CFGs_PP_MD(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.PP_MD(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ works....
    def test_CFGs_PP_FI(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.PP_FI(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ works....
    def test_CFGs_PP_FI_TM(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.PP_FI_TM(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ works...
    def test_CFGs_MSVS_trivial(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T = tdlib.trivial_decomposition(V, E)
            V_T, E_T, w = tdlib.MSVS(V, E, V_T, E_T)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ works...
    def test_CFGs_minimalChordal_trivial(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T = tdlib.trivial_decomposition(V, E)
            V_T, E_T, w = tdlib.minimalChordal_decomp(V, E, V_T, E_T)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ works...
    def test_CFGs_seperator_algorithm(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.seperator_algorithm(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ works...
    def test_CFGs_conversion(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w1 = tdlib.minDegree_decomp(V, E)
            O = tdlib.treedec_to_ordering(V_T, E_T)
            V_T, E_T, w2 = tdlib.ordering_to_treedec(V, E, O)
            self.assertEqual(w1, w2)
    """

    """ LBs """

    """ works...
    #indirect test
    def test_CFGs_LB1(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            tdlib.lower_bound(V, E, "deltaC_min_d")

    #indirect test
    def test_CFGs_LB2(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            tdlib.lower_bound(V, E, "deltaC_max_d")

    #indirect test
    def test_CFGs_LB3(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            tdlib.lower_bound(V, E, "deltaC_least_c")

    #indirect test
    def test_CFGs_LB3(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            tdlib.lower_bound(V, E, "LBN_deltaC")

    #indirect test
    def test_CFGs_LB3(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            tdlib.lower_bound(V, E, "LBNC_deltaC")

    #indirect test
    def test_CFGs_LB3(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            tdlib.lower_bound(V, E, "LBP_deltaC")

    #indirect test
    def test_CFGs_LB3(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            tdlib.lower_bound(V, E, "LBPC_deltaC")
    """

    """ takes long time...
    def test_CFGs_exact_cutset(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.exact_decomposition_cutset(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """


    """ APPs """

    #validation (is_clique, is_IS, is_VC,..) in tdlib?

    """ works...
    def test_CFGs_max_clique(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.minDegree_decomp(V, E)
            S = tdlib.max_clique_with_treedecomposition(V, E, V_T, E_T)
    """

    """ works...
    def test_CFGs_max_independent_set(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.minDegree_decomp(V, E)
            S = tdlib.max_independent_set_with_treedecomposition(V, E, V_T, E_T)
    """

    """ works...
    def test_CFGs_min_vertex_cover(self):
        for i in range(0, CFGS_count+1):
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.minDegree_decomp(V, E)
            S = tdlib.min_vertex_cover_with_treedecomposition(V, E, V_T, E_T)
    """

    """ works...
    def test_CFGs_min_dominating_set(self):
        for i in range(0, CFGS_count+1):
            if i == 999: #huge graph
                continue;
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.minDegree_decomp(V, E)
            S = tdlib.min_dominating_set_with_treedecomposition(V, E, V_T, E_T)
    """

    """ works...
    def test_CFGs_min_coloring(self):
        for i in range(0, CFGS_count+1):
            if i == 999: #huge graph
                continue
            V = eval("CFGs.V_" + str(i))
            E = eval("CFGs.E_" + str(i))
            name = eval("CFGs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.minDegree_decomp(V, E)
            S = tdlib.min_coloring_with_treedecomposition(V, E, V_T, E_T)
    """

    """ works...
    def test_DIMACS_MD(self):
        for i in range(0, DIMACS_count+1):
            V = eval("dimacs.V_" + str(i))
            E = eval("dimacs.E_" + str(i))
            name = eval("dimacs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.minDegree_decomp(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ works...
    def test_DIMACS_boost_MD(self):
        for i in range(0, DIMACS_count+1):
            V = eval("dimacs.V_" + str(i))
            E = eval("dimacs.E_" + str(i))
            name = eval("dimacs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.boost_minDegree_decomp(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ ?
    def test_DIMACS_FI(self):
        for i in range(0, DIMACS_count+1):
            V = eval("dimacs.V_" + str(i))
            E = eval("dimacs.E_" + str(i))
            name = eval("dimacs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.fillIn_decomp(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ ?
    def test_DIMACS_PP_MD(self):
        for i in range(0, DIMACS_count+1):
            V = eval("dimacs.V_" + str(i))
            E = eval("dimacs.E_" + str(i))
            name = eval("dimacs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.PP_MD(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ ?
    def test_DIMACS_PP_FI(self):
        for i in range(0, DIMACS_count+1):
            V = eval("dimacs.V_" + str(i))
            E = eval("dimacs.E_" + str(i))
            name = eval("dimacs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.PP_FI(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

    """ ?
    def test_DIMACS_PP_FI_TM(self):
        for i in range(0, DIMACS_count+1):
            V = eval("dimacs.V_" + str(i))
            E = eval("dimacs.E_" + str(i))
            name = eval("dimacs.name_" + str(i))

            if VERBOSE:
                print(str(i) + ": " + name)
                print("n: " + str(len(V)))
                print("e: " + str(len(E)))

            V_T, E_T, w = tdlib.PP_FI_TM(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)
    """

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
