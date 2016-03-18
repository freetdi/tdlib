import tdlib
import unittest

#some graphs
V_P6 = [0,1,2,3,4,5]
E_P6 = [(0,1),(1,2),(2,3),(3,4),(4,5)]

V_K5 = [0,1,2,3,4]
E_K5 = [(0,1),(0,2),(0,3),(0,4),(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]

V_Petersen = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
E_Petersen = [(0, 1), (0, 4), (0, 5), (1, 2), (1, 6), (2, 3), (2, 7), \
              (3, 4), (3, 8), (4, 9), (5, 7), (5, 8), (6, 8), (6, 9), (7, 9)]

V_Wagner = [0, 1, 2, 3, 4, 5, 6, 7]
E_Wagner = [(0, 1),(0, 4), (0, 7), (1, 2), (1, 5), (2, 3), (2, 6), (3, 4), \
            (3, 7), (4, 5), (5, 6), (6, 7)]

V_Pappus = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
E_Pappus = [(0, 1), (0, 5), (0, 6), (1, 2), (1, 7), (2, 3), (2, 8), (3, 4), \
            (3, 9), (4, 5), (4, 10), (5, 11), (6, 13), (6, 17), (7, 12), \
            (7, 14), (8, 13), (8, 15), (9, 14), (9, 16), (10, 15), (10, 17), \
            (11, 12), (11, 16), (12, 15), (13, 16), (14, 17)]


#random graphs
def randomGNP(n, p):
    import random
    V = []
    E = []
    for i in range(0, n):
        V.append(i)
        for j in range(i+1,n):
            if random.randint(0, 100) <= p*100:
                E.append((i, j))
    return V, E


class TestTdLib(unittest.TestCase):
    def test_preprocessing(self):
        V_, E_, B, lb = tdlib.preprocessing(V_P6, E_P6)
        self.assertEqual(V_, [])
        self.assertEqual(E_, [])
        self.assertEqual(B, [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5]])
        self.assertEqual(lb, 1)

        V_, E_, B, lb = tdlib.preprocessing(V_K5, E_K5)
        self.assertEqual(V_, [])
        self.assertEqual(E_, [])
        self.assertEqual(B, [[0, 1, 2, 3, 4], [1, 2, 3, 4], [2, 3, 4], [3, 4]])
        self.assertEqual(lb, 4)

        V_, E_, B, lb = tdlib.preprocessing(V_Petersen, E_Petersen)
        self.assertEqual(V_, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        self.assertEqual(E_, [0,1,0,4,0,5,1,2,1,6,2,3,2,7,3, \
                              4,3,8,4,9,5,7,5,8,6,8,6,9,7,9])
        self.assertEqual(B, [])
        self.assertEqual(lb, 4)

    def test_PP_MD(self):
        V, E, lb = tdlib.PP_MD(V_P6, E_P6)
        self.assertEqual(V, [[4, 5], [3, 4], [2, 3], [1, 2], [0, 1]])
        self.assertEqual(E, [0, 1, 1, 2, 2, 3, 3, 4])
        self.assertEqual(lb, 1)

        V, E, lb = tdlib.PP_MD(V_K5, E_K5)
        self.assertEqual(V, [[0, 1, 2, 3, 4]])
        self.assertEqual(E, [])
        self.assertEqual(lb, 4)

        V, E, lb = tdlib.PP_MD(V_Petersen, E_Petersen)
        self.assertEqual(V, [[1, 4, 7, 8, 9], [1, 4, 5, 7, 8], [1, 3, 4, 7, 8], \
                             [1, 6, 8, 9], [1, 2, 3, 7], [0, 1, 4, 5]])
        self.assertEqual(E, [0, 1, 0, 2, 0, 3, 2, 4, 1, 5])
        self.assertEqual(lb, 4)

        V, E, lb = tdlib.PP_MD(V_Pappus, E_Pappus)
        self.assertEqual(V, [[3, 5, 6, 8, 12, 14, 16], [3, 5, 6, 8, 10, 12, 14], \
                             [1, 3, 5, 6, 8, 12, 14], [6, 10, 14, 17], \
                             [8, 10, 12, 15], [6, 8, 13, 16], [5, 11, 12, 16], \
                             [3, 9, 14, 16], [1, 7, 12, 14], [3, 4, 5, 10], \
                             [1, 2, 3, 8], [0, 1, 5, 6]])
        self.assertEqual(E, [0, 1, 0, 2, 1, 3, 1, 4, 0, 7, 2, 8, 1, 9, 2, \
                             10, 2, 11, 0, 6, 0, 5])
        self.assertEqual(lb, 6)

    def test_PP_FI_TM(self):
        V, E, lb = tdlib.PP_FI_TM(V_P6, E_P6)
        self.assertEqual(V, [[4, 5], [3, 4], [2, 3], [1, 2], [0, 1]])
        self.assertEqual(E, [0, 1, 1, 2, 2, 3, 3, 4])
        self.assertEqual(lb, 1)

        V, E, lb = tdlib.PP_FI_TM(V_K5, E_K5)
        self.assertEqual(V, [[0, 1, 2, 3, 4]])
        self.assertEqual(E, [])
        self.assertEqual(lb, 4)

        V, E, lb = tdlib.PP_FI_TM(V_Petersen, E_Petersen)
        self.assertEqual(V, [[0, 1, 4, 5], [1, 4, 5, 7, 8], [1, 3, 4, 7, 8], \
                             [1, 4, 7, 8, 9], [1, 2, 3, 7], [1, 6, 8, 9]])
        self.assertEqual(E, [1, 2, 1, 3, 2, 4, 3, 5, 1, 0])
        self.assertEqual(lb, 4)

        V, E, lb = tdlib.PP_FI_TM(V_Pappus, E_Pappus)
        self.assertEqual(V, [[0, 1, 5, 6], [1, 3, 5, 6, 8, 12, 14], \
                             [3, 5, 6, 8, 10, 12, 14], [3, 5, 6, 8, 12, 14, 16], \
                             [1, 2, 3, 8], [3, 4, 5, 10], [1, 7, 12, 14], \
                             [3, 9, 14, 16], [5, 11, 12, 16], [6, 8, 13, 16], \
                             [8, 10, 12, 15], [6, 10, 14, 17]])
        self.assertEqual(E, [1, 2, 1, 3, 2, 5, 1, 6, 3, 7, 3, 8, 3, 9, 2, \
                             10, 2, 11, 1, 4, 1, 0])
        self.assertEqual(lb, 6)

    def test_lower_bounds(self):
        lb = tdlib.lower_bound(V_Wagner, E_Wagner, "deltaC_min_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Wagner, E_Wagner, "deltaC_max_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Wagner, E_Wagner, "deltaC_least_c")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Wagner, E_Wagner, "LBN_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Wagner, E_Wagner, "LBNC_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Wagner, E_Wagner, "LBP_deltaC")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Wagner, E_Wagner, "LBPC_deltaC")
        self.assertEqual(lb, 3)

        lb = tdlib.lower_bound(V_Pappus, E_Pappus, "deltaC_min_d")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Pappus, E_Pappus, "deltaC_max_d")
        self.assertEqual(lb, 3)
        lb = tdlib.lower_bound(V_Pappus, E_Pappus, "deltaC_least_c")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Pappus, E_Pappus, "LBN_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Pappus, E_Pappus, "LBNC_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Pappus, E_Pappus, "LBP_deltaC")
        self.assertEqual(lb, 4)
        lb = tdlib.lower_bound(V_Pappus, E_Pappus, "LBPC_deltaC")
        self.assertEqual(lb, 4)

    def test_exact_decomposition_cutset(self):
        V, E, tw = tdlib.exact_decomposition_cutset(V_P6, E_P6)
        self.assertEqual(V, [[4, 5], [3, 4], [2, 3], [1, 2], [0, 1]])
        self.assertEqual(E, [0, 1, 1, 2, 2, 3, 3, 4])
        self.assertEqual(tw, 1)

        V, E, tw = tdlib.exact_decomposition_cutset(V_K5, E_K5)
        self.assertEqual(V, [[0, 1, 2, 3, 4]])
        self.assertEqual(E, [])
        self.assertEqual(tw, 4)

        V, E, tw = tdlib.exact_decomposition_cutset(V_Petersen, E_Petersen)
        self.assertEqual(V, [[4, 6, 7, 9], [0, 3, 4, 6, 7], [0, 2, 3, 6, 7], \
                             [3, 5, 6, 8], [0, 3, 5, 6, 7], [0, 1, 2, 3, 6]])
        self.assertEqual(E, [0, 1, 1, 2, 1, 4, 2, 5, 3, 4])
        self.assertEqual(tw, 4)

        V, E, tw = tdlib.exact_decomposition_cutset(V_Wagner, E_Wagner)
        self.assertEqual(V, [[0, 3, 4, 5], [0, 1, 2, 3, 5], [0, 3, 6, 7], \
                             [0, 2, 3, 5, 6]])
        self.assertEqual(E, [0, 1, 2, 3, 3, 1])
        self.assertEqual(tw, 4)

    """
    def test_exact_decomposition_dynamic(self):
        V, E, tw = tdlib.exact_decomposition_dynamic(V_P6, E_P6)
        self.assertEqual(V, [[4, 5], [3, 4], [2, 3], [1, 2], [0, 1]])
        self.assertEqual(E, [0, 1, 1, 2, 2, 3, 3, 4])
        self.assertEqual(tw, 1)

        V, E, tw = tdlib.exact_decomposition_dynamic(V_K5, E_K5)
        self.assertEqual(V, [[0, 1, 2, 3, 4]])
        self.assertEqual(E, [])
        self.assertEqual(tw, 4)
    """

    def test_seperator_algorithm(self):
        V, E, lb = tdlib.seperator_algorithm(V_P6, E_P6)
        self.assertEqual(V, [[0, 1, 2, 3], [1, 3, 4, 5]])
        self.assertEqual(E, [1,0])
        self.assertEqual(lb, 3)

        V, E, lb = tdlib.seperator_algorithm(V_K5, E_K5)
        self.assertEqual(V, [[0, 1, 2, 3, 4]])
        self.assertEqual(E, [])
        self.assertEqual(lb, 4)

        V, E, lb = tdlib.seperator_algorithm(V_Pappus, E_Pappus)
        self.assertEqual(V, [[0, 1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6, 7], \
                             [2, 3, 4, 5, 6, 7, 8], [3, 4, 5, 6, 7, 8, 9], \
                             [4, 5, 6, 7, 8, 9, 10], [5, 6, 7, 8, 9, 10, 11], \
                             [6, 7, 8, 9, 10, 11, 12, 14], \
                             [6, 8, 9, 10, 11, 12, 13, 14], [6, 10, 13, 14, 17], \
                             [8, 9, 10, 11, 12, 13, 14, 15, 16]])
        self.assertEqual(E, [1, 0, 2, 1, 3, 2, 4, 3, 5, 4, 6, 5, 7, 6, 8, 7, 9, 7])
        self.assertEqual(lb, 8)

    def test_minDegree_ordering(self):
        O = tdlib.minDegree_ordering(V_Petersen, E_Petersen)
        self.assertEqual(O, [0, 2, 6, 3, 5, 1, 4, 7, 8, 9])
        O = tdlib.minDegree_ordering(V_Pappus, E_Pappus)
        self.assertEqual(O, [0, 2, 4, 7, 9, 11, 13, 15, 17, 1, 10, 3, \
                             5, 6, 8, 12, 14, 16])

    def test_fillIn_ordering(self):
        O = tdlib.fillIn_ordering(V_Petersen, E_Petersen)
        self.assertEqual(O, [0, 2, 6, 3, 5, 1, 4, 7, 8, 9])
        O = tdlib.fillIn_ordering(V_Pappus, E_Pappus)
        self.assertEqual(O, [0, 2, 4, 7, 9, 11, 13, 15, 17, 1, \
                             10, 3, 5, 6, 8, 12, 14, 16])

    """
    def test_treedec_to_ordering(self):
        V, E, lb = tdlib.seperator_algorithm(V_P6, E_P6)
        O = tdlib.treedec_to_ordering(V, E)
        self.assertEqual(O, [0, 2, 1, 3, 4, 5])
    """

    def test_ordering_to_treedec(self):
        V, E, lb = tdlib.ordering_to_treedec(V_P6, E_P6, [1,3,0,2,4,5])
        self.assertEqual(V, [[5], [4, 5], [2, 4], [0, 2], [2, 3, 4], [0, 1, 2]])
        self.assertEqual(E, [0, 1, 1, 2, 2, 3, 2, 4, 3, 5])

    def test_trivial_decomposition(self):
        V, E = tdlib.trivial_decomposition(["a", "b", "c"], [])
        self.assertEqual(V, [["a", "b", "c"]])

        V, E = tdlib.trivial_decomposition(V_P6, E_P6)
        self.assertEqual(V, [[0,1,2,3,4,5]])

        V, E = tdlib.trivial_decomposition(V_K5, E_K5)
        self.assertEqual(V, [[0,1,2,3,4]])

    """
    def test_MSVS(self):
        V, E = tdlib.trivial_decomposition(V_P6, E_P6)
        V, E, w = tdlib.MSVS(V_P6, E_P6, V, E)
        self.assertEqual(w, 1)
    """

    """
    def test_is_valid_treedecomposition(self):
        V, E = tdlib.trivial_decomposition(V_P6, E_P6)
        
        status = tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E, message=False)
        self.assertEqual(status, 0)

        status = tdlib.is_valid_treedecomposition(V_K5, E_K5, V, E, message=False)
        self.assertEqual(status, -2)

        V, E, tw = tdlib.exact_decomposition_cutset(V_Petersen, E_Petersen) 
        status = tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E, message=False)
        self.assertEqual(status, 0)

        V, E, tw = tdlib.exact_decomposition_cutset(V_Wagner, E_Wagner) 
        status = tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V, E, message=False)
        self.assertEqual(status, 0)
    """

    """
    def test_random_valid_treedecomposition(self):
        status_list = list()
        correct_status = list()
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.3)
                N, M, tw = tdlib.exact_decomposition_cutset(V, E)
                status = tdlib.is_valid_treedecomposition(V, E, N, M, message=True)
                status_list.append(status)
                correct_status.append(0)
                if(status < 0):
                    print("error (validate decompositions)! graph: " + str(V) + ", " + str(E))
    """

    def test_random_validate_width(self):
        status = True
        for n in range(0, 13):
            for i in range(0, 100):
                V, E = randomGNP(n, 0.3)
                Q, R, w2 = tdlib.PP_FI_TM(V, E)
                isleq = tdlib.exact_decomposition_cutset_decision(V, E, w2)
                if(not isleq):
                    print("error (validate width)! graph: " + str(V) + ", " + str(E))
                    status = False

        self.assertEqual(status, True)

    """
    def test_independent_set_with_treedecomposition(self):
        V_T1, E_T1, lb = tdlib.PP_MD(V_Petersen, E_Petersen)
        IS1 = tdlib.max_independent_set_with_treedecomposition(V_Petersen, E_Petersen, V_T1, E_T1)
        self.assertEqual(len(IS1), 4)

        V_T2, E_T2, lb = tdlib.PP_MD(V_Wagner, E_Wagner)
        IS2 = tdlib.max_independent_set_with_treedecomposition(V_Wagner, E_Wagner, V_T2, E_T2)
        self.assertEqual(len(IS2), 3)

        V_T3, E_T3, lb = tdlib.PP_MD(V_Pappus, E_Pappus)
        IS3 = tdlib.max_independent_set_with_treedecomposition(V_Pappus, E_Pappus, V_T3, E_T3)
        self.assertEqual(len(IS3), 9)

    def test_vertex_cover_with_treedecomposition(self):
        V_T1, E_T1, lb = tdlib.PP_MD(V_Petersen, E_Petersen)
        VC1 = tdlib.min_vertex_cover_with_treedecomposition(V_Petersen, E_Petersen, V_T1, E_T1)
        self.assertEqual(len(VC1), 6)

        V_T2, E_T2, lb = tdlib.PP_MD(V_Wagner, E_Wagner)
        VC2 = tdlib.min_vertex_cover_with_treedecomposition(V_Wagner, E_Wagner, V_T2, E_T2)
        self.assertEqual(len(VC2), 5)

        V_T3, E_T3, lb = tdlib.PP_MD(V_Pappus, E_Pappus)
        VC3 = tdlib.min_vertex_cover_with_treedecomposition(V_Pappus, E_Pappus, V_T3, E_T3)
        self.assertEqual(len(VC3), 9)

    def test_dominating_set_with_treedecomposition(self):
        V_T1, E_T1, lb = tdlib.PP_MD(V_Petersen, E_Petersen)
        DS1 = tdlib.min_dominating_set_with_treedecomposition(V_Petersen, E_Petersen, V_T1, E_T1)
        self.assertEqual(len(DS1), 3)

        V_T2, E_T2, lb = tdlib.PP_MD(V_Wagner, E_Wagner)
        DS2 = tdlib.min_dominating_set_with_treedecomposition(V_Wagner, E_Wagner, V_T2, E_T2)
        self.assertEqual(len(DS2), 3)

        V_T3, E_T3, lb = tdlib.PP_MD(V_Pappus, E_Pappus)
        DS3 = tdlib.min_dominating_set_with_treedecomposition(V_Pappus, E_Pappus, V_T3, E_T3)
        self.assertEqual(len(DS3), 5)
    """

if __name__ == '__main__':
    unittest.main()

