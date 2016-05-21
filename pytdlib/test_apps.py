import tdlib
import unittest

from graphs import *

class TestTdLib_app(unittest.TestCase):

    def test_max_clique_with_treedecomposition(self):
        V_T1, E_T1, lb = tdlib.PP_MD(V_Petersen, E_Petersen)
        C1 = tdlib.max_clique_with_treedecomposition(V_Petersen, E_Petersen, V_T1, E_T1)
        self.assertEqual(len(C1), 2)

        V_T2, E_T2, lb = tdlib.PP_MD(V_Wagner, E_Wagner)
        C2 = tdlib.max_clique_with_treedecomposition(V_Wagner, E_Wagner, V_T2, E_T2)
        self.assertEqual(len(C2), 2)

        V_T3, E_T3, lb = tdlib.PP_MD(V_Pappus, E_Pappus)
        C3 = tdlib.max_clique_with_treedecomposition(V_Pappus, E_Pappus, V_T3, E_T3)
        self.assertEqual(len(C3), 2)

        V_T4, E_T4, lb = tdlib.PP_MD(V_K5, E_K5)
        C4 = tdlib.max_clique_with_treedecomposition(V_K5, E_K5, V_T4, E_T4)
        self.assertEqual(len(C4), 5)

    def test_max_independent_set_with_treedecomposition(self):
        V_T1, E_T1, lb = tdlib.PP_MD(V_Petersen, E_Petersen)
        IS1 = tdlib.max_independent_set_with_treedecomposition(V_Petersen, E_Petersen, V_T1, E_T1)
        self.assertEqual(len(IS1), 4)

        V_T2, E_T2, lb = tdlib.PP_MD(V_Wagner, E_Wagner)
        IS2 = tdlib.max_independent_set_with_treedecomposition(V_Wagner, E_Wagner, V_T2, E_T2)
        self.assertEqual(len(IS2), 3)

        V_T3, E_T3, lb = tdlib.PP_MD(V_Pappus, E_Pappus)
        IS3 = tdlib.max_independent_set_with_treedecomposition(V_Pappus, E_Pappus, V_T3, E_T3)
        self.assertEqual(len(IS3), 9)

        V_T4, E_T4, lb = tdlib.PP_MD(V_K5, E_K5)
        IS4 = tdlib.max_independent_set_with_treedecomposition(V_K5, E_K5, V_T4, E_T4)
        self.assertEqual(len(IS4), 1)

    def test_min_vertex_cover_with_treedecomposition(self):
        V_T1, E_T1, lb = tdlib.PP_MD(V_Petersen, E_Petersen)
        VC1 = tdlib.min_vertex_cover_with_treedecomposition(V_Petersen, E_Petersen, V_T1, E_T1)
        self.assertEqual(len(VC1), 6)

        V_T2, E_T2, lb = tdlib.PP_MD(V_Wagner, E_Wagner)
        VC2 = tdlib.min_vertex_cover_with_treedecomposition(V_Wagner, E_Wagner, V_T2, E_T2)
        self.assertEqual(len(VC2), 5)

        V_T3, E_T3, lb = tdlib.PP_MD(V_Pappus, E_Pappus)
        VC3 = tdlib.min_vertex_cover_with_treedecomposition(V_Pappus, E_Pappus, V_T3, E_T3)
        self.assertEqual(len(VC3), 9)

        V_T4, E_T4, lb = tdlib.PP_MD(V_K5, E_K5)
        VC4 = tdlib.min_vertex_cover_with_treedecomposition(V_K5, E_K5, V_T4, E_T4)
        self.assertEqual(len(VC4), 4)

    def test_min_dominating_set_with_treedecomposition(self):
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
    def test_min_feedback_vertex_set_with_treedecomposition(self):
        V_T1, E_T1, lb = tdlib.PP_MD(V_Petersen, E_Petersen)
        FVS1 = tdlib.min_feedback_vertex_set_with_treedecomposition(V_Petersen, E_Petersen, V_T1, E_T1)
        self.assertEqual(len(FVS1), 3)

        V_T2, E_T2, lb = tdlib.PP_MD(V_Wagner, E_Wagner)
        FVS2 = tdlib.min_feedback_vertex_set_with_treedecomposition(V_Wagner, E_Wagner, V_T2, E_T2)
        self.assertEqual(len(FVS2), 3)

        V_T3, E_T3, lb = tdlib.PP_MD(V_Pappus, E_Pappus)
        FVS3 = tdlib.min_feedback_vertex_set_with_treedecomposition(V_Pappus, E_Pappus, V_T3, E_T3)
        self.assertEqual(len(FVS3), 5)
    """

    def test_min_coloring_with_treedecomposition(self):
        V_T1, E_T1, lb = tdlib.PP_MD(V_Petersen, E_Petersen)
        C = tdlib.min_coloring_with_treedecomposition(V_Petersen, E_Petersen, V_T1, E_T1)
        self.assertEqual(len(C), 3)

        V_T2, E_T2, lb = tdlib.PP_MD(V_Wagner, E_Wagner)
        C = tdlib.min_coloring_with_treedecomposition(V_Wagner, E_Wagner, V_T2, E_T2)
        self.assertEqual(len(C), 3)

        V_T3, E_T3, lb = tdlib.PP_MD(V_Pappus, E_Pappus)
        C = tdlib.min_coloring_with_treedecomposition(V_Pappus, E_Pappus, V_T3, E_T3)
        self.assertEqual(len(C), 2)

        V_T4, E_T4, lb = tdlib.PP_MD(V_K5, E_K5)
        C = tdlib.min_coloring_with_treedecomposition(V_K5, E_K5, V_T4, E_T4)
        self.assertEqual(len(C), 5)


if __name__ == '__main__':
    unittest.main()
