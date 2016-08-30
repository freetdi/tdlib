import tdlib
import unittest
import random

from graphs import *

class TestTdLib_app(unittest.TestCase):
    def test_max_clique_with_treedecomposition_0a(self):
        V, E = cornercases[0]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.max_clique_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 0)

    def test_max_clique_with_treedecomposition_0b(self):
        V, E = cornercases[1]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.max_clique_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 1)

    def test_max_clique_with_treedecomposition_0c(self):
        V, E = cornercases[2]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.max_clique_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 1)

    def test_max_clique_with_treedecomposition_0d(self):
        V, E = cornercases[3]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.max_clique_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 5)

    def test_max_clique_with_treedecomposition_1(self):
        V_T, E_T, w = tdlib.PP_MD(V_P6, E_P6)
        S = tdlib.max_clique_with_treedecomposition(V_P6, E_P6, V_T, E_T)
        self.assertEqual(len(S), 2)

    def test_max_clique_with_treedecomposition_2(self):
        V_T, E_T, w = tdlib.PP_MD(V_K5, E_K5)
        S = tdlib.max_clique_with_treedecomposition(V_K5, E_K5, V_T, E_T)
        self.assertEqual(len(S), 5)

    def test_max_clique_with_treedecomposition_3(self):
        V_T, E_T, w = tdlib.PP_MD(V_Petersen, E_Petersen)
        S = tdlib.max_clique_with_treedecomposition(V_Petersen, E_Petersen, V_T, E_T)
        self.assertEqual(len(S), 2)

    def test_max_clique_with_treedecomposition_4(self):
        V_T, E_T, w = tdlib.PP_MD(V_Petersen_double, E_Petersen_double)
        S = tdlib.max_clique_with_treedecomposition(V_Petersen_double, E_Petersen_double, V_T, E_T)
        self.assertEqual(len(S), 2)

    def test_max_clique_with_treedecomposition_5(self):
        V_T, E_T, w = tdlib.PP_MD(V_Wagner, E_Wagner)
        S = tdlib.max_clique_with_treedecomposition(V_Wagner, E_Wagner, V_T, E_T)
        self.assertEqual(len(S), 2)

    def test_max_clique_with_treedecomposition_6(self):
        V_T, E_T, w = tdlib.PP_MD(V_Pappus, E_Pappus)
        S = tdlib.max_clique_with_treedecomposition(V_Pappus, E_Pappus, V_T, E_T)
        self.assertEqual(len(S), 2)

    def test_max_clique_with_treedecomposition_7(self):
        V_T, E_T, w = tdlib.PP_MD(V_Grid_5_5, E_Grid_5_5)
        S = tdlib.max_clique_with_treedecomposition(V_Grid_5_5, E_Grid_5_5, V_T, E_T)
        self.assertEqual(len(S), 2)

    def test_max_clique_with_treedecomposition_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                V_T, E_T, w = tdlib.PP_MD(V, E)
                S = tdlib.max_clique_with_treedecomposition(V, E, V_T, E_T)


    def test_max_independent_set_with_treedecomposition_0a(self):
        V, E = cornercases[0]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.max_independent_set_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 0)

    def test_max_independent_set_with_treedecomposition_0b(self):
        V, E = cornercases[1]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.max_independent_set_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 1)

    def test_max_independent_set_with_treedecomposition_0c(self):
        V, E = cornercases[2]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.max_independent_set_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 5)

    def test_max_independent_set_with_treedecomposition_0d(self):
        V, E = cornercases[3]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.max_independent_set_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 1)

    def test_max_independent_set_with_treedecomposition_1(self):
        V_T, E_T, w = tdlib.PP_MD(V_P6, E_P6)
        S = tdlib.max_independent_set_with_treedecomposition(V_P6, E_P6, V_T, E_T)
        self.assertEqual(len(S), 3)

    def test_max_independent_set_with_treedecomposition_2(self):
        V_T, E_T, w = tdlib.PP_MD(V_K5, E_K5)
        S = tdlib.max_independent_set_with_treedecomposition(V_K5, E_K5, V_T, E_T)
        self.assertEqual(len(S), 1)

    def test_max_independent_set_with_treedecomposition_3(self):
        V_T, E_T, w = tdlib.PP_MD(V_Petersen, E_Petersen)
        S = tdlib.max_independent_set_with_treedecomposition(V_Petersen, E_Petersen, V_T, E_T)
        self.assertEqual(len(S), 4)

    def test_max_independent_set_with_treedecomposition_4(self):
        V_T, E_T, w = tdlib.PP_MD(V_Petersen_double, E_Petersen_double)
        S = tdlib.max_independent_set_with_treedecomposition(V_Petersen_double, E_Petersen_double, V_T, E_T)
        self.assertEqual(len(S), 8)

    def test_max_independent_set_with_treedecomposition_5(self):
        V_T, E_T, w = tdlib.PP_MD(V_Wagner, E_Wagner)
        S = tdlib.max_independent_set_with_treedecomposition(V_Wagner, E_Wagner, V_T, E_T)
        self.assertEqual(len(S), 3)

    def test_max_independent_set_with_treedecomposition_6(self):
        V_T, E_T, w = tdlib.PP_MD(V_Pappus, E_Pappus)
        S = tdlib.max_independent_set_with_treedecomposition(V_Pappus, E_Pappus, V_T, E_T)
        self.assertEqual(len(S), 9)

    def test_max_independent_set_with_treedecomposition_7(self):
        V_T, E_T, w = tdlib.PP_MD(V_Grid_5_5, E_Grid_5_5)
        S = tdlib.max_independent_set_with_treedecomposition(V_Grid_5_5, E_Grid_5_5, V_T, E_T)
        self.assertEqual(len(S), 13)

    def test_max_independent_set_with_treedecomposition_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                V_T, E_T, w = tdlib.PP_MD(V, E)
                S = tdlib.max_independent_set_with_treedecomposition(V, E, V_T, E_T)


    def test_min_vertex_cover_with_treedecomposition_0a(self):
        V, E = cornercases[0]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.min_vertex_cover_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 0)

    def test_min_vertex_cover_with_treedecomposition_0b(self):
        V, E = cornercases[1]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.min_vertex_cover_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 0)

    def test_min_vertex_cover_with_treedecomposition_0c(self):
        V, E = cornercases[2]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.min_vertex_cover_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 0)

    def test_min_vertex_cover_with_treedecomposition_0d(self):
        V, E = cornercases[3]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.min_vertex_cover_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 4)

    def test_min_vertex_cover_with_treedecomposition_1(self):
        V_T, E_T, w = tdlib.PP_MD(V_P6, E_P6)
        S = tdlib.min_vertex_cover_with_treedecomposition(V_P6, E_P6, V_T, E_T)
        self.assertEqual(len(S), 3)

    def test_min_vertex_cover_with_treedecomposition_2(self):
        V_T, E_T, w = tdlib.PP_MD(V_K5, E_K5)
        S = tdlib.min_vertex_cover_with_treedecomposition(V_K5, E_K5, V_T, E_T)
        self.assertEqual(len(S), 4)

    def test_min_vertex_cover_with_treedecomposition_3(self):
        V_T, E_T, w = tdlib.PP_MD(V_Petersen, E_Petersen)
        S = tdlib.min_vertex_cover_with_treedecomposition(V_Petersen, E_Petersen, V_T, E_T)
        self.assertEqual(len(S), 6)

    def test_min_vertex_cover_with_treedecomposition_4(self):
        V_T, E_T, w = tdlib.PP_MD(V_Petersen_double, E_Petersen_double)
        S = tdlib.min_vertex_cover_with_treedecomposition(V_Petersen_double, E_Petersen_double, V_T, E_T)
        self.assertEqual(len(S), 12)

    def test_min_vertex_cover_with_treedecomposition_5(self):
        V_T, E_T, w = tdlib.PP_MD(V_Wagner, E_Wagner)
        S = tdlib.min_vertex_cover_with_treedecomposition(V_Wagner, E_Wagner, V_T, E_T)
        self.assertEqual(len(S), 5)

    def test_min_vertex_cover_with_treedecomposition_6(self):
        V_T, E_T, w = tdlib.PP_MD(V_Pappus, E_Pappus)
        S = tdlib.min_vertex_cover_with_treedecomposition(V_Pappus, E_Pappus, V_T, E_T)
        self.assertEqual(len(S), 9)

    def test_min_vertex_cover_with_treedecomposition_7(self):
        V_T, E_T, w = tdlib.PP_MD(V_Grid_5_5, E_Grid_5_5)
        S = tdlib.min_vertex_cover_with_treedecomposition(V_Grid_5_5, E_Grid_5_5, V_T, E_T)
        self.assertEqual(len(S), 12)

    def test_min_vertex_cover_with_treedecomposition_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                V_T, E_T, w = tdlib.PP_MD(V, E)
                S = tdlib.min_vertex_cover_with_treedecomposition(V, E, V_T, E_T)


    def test_min_dominating_set_with_treedecomposition_0a(self):
        V, E = cornercases[0]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.min_dominating_set_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 0)

    def test_min_dominating_set_with_treedecomposition_0b(self):
        V, E = cornercases[1]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.min_dominating_set_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 1)

    def test_min_dominating_set_with_treedecomposition_0c(self):
        V, E = cornercases[2]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.min_dominating_set_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 5)

    def test_min_dominating_set_with_treedecomposition_0d(self):
        V, E = cornercases[3]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.min_dominating_set_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 1)

    def test_min_dominating_set_with_treedecomposition_1(self):
        V_T, E_T, w = tdlib.PP_MD(V_P6, E_P6)
        S = tdlib.min_dominating_set_with_treedecomposition(V_P6, E_P6, V_T, E_T)
        self.assertEqual(len(S), 2)

    def test_min_dominating_set_with_treedecomposition_2(self):
        V_T, E_T, w = tdlib.PP_MD(V_K5, E_K5)
        S = tdlib.min_dominating_set_with_treedecomposition(V_K5, E_K5, V_T, E_T)
        self.assertEqual(len(S), 1)

    def test_min_dominating_set_with_treedecomposition_3(self):
        V_T, E_T, w = tdlib.PP_MD(V_Petersen, E_Petersen)
        S = tdlib.min_dominating_set_with_treedecomposition(V_Petersen, E_Petersen, V_T, E_T)
        self.assertEqual(len(S), 3)

    def test_min_dominating_set_with_treedecomposition_4(self):
        V_T, E_T, w = tdlib.PP_MD(V_Petersen_double, E_Petersen_double)
        S = tdlib.min_dominating_set_with_treedecomposition(V_Petersen_double, E_Petersen_double, V_T, E_T)
        self.assertEqual(len(S), 6)

    def test_min_dominating_set_with_treedecomposition_5(self):
        V_T, E_T, w = tdlib.PP_MD(V_Wagner, E_Wagner)
        S = tdlib.min_dominating_set_with_treedecomposition(V_Wagner, E_Wagner, V_T, E_T)
        self.assertEqual(len(S), 3)

    def test_min_dominating_set_with_treedecomposition_6(self):
        V_T, E_T, w = tdlib.PP_MD(V_Pappus, E_Pappus)
        S = tdlib.min_dominating_set_with_treedecomposition(V_Pappus, E_Pappus, V_T, E_T)
        self.assertEqual(len(S), 5)

    def test_min_dominating_set_with_treedecomposition_7(self):
        V_T, E_T, w = tdlib.PP_MD(V_Grid_5_5, E_Grid_5_5)
        S = tdlib.min_dominating_set_with_treedecomposition(V_Grid_5_5, E_Grid_5_5, V_T, E_T)
        self.assertEqual(len(S), 7)

    def test_min_dominating_set_with_treedecomposition_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                V_T, E_T, w = tdlib.PP_MD(V, E)
                S = tdlib.min_dominating_set_with_treedecomposition(V, E, V_T, E_T)


    def test_min_coloring_with_treedecomposition_0a(self):
        V, E = cornercases[0]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.min_coloring_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 0)

    def test_min_coloring_with_treedecomposition_0b(self):
        V, E = cornercases[1]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.min_coloring_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 1)

    def test_min_coloring_with_treedecomposition_0c(self):
        V, E = cornercases[2]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.min_coloring_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 1)

    def test_min_coloring_with_treedecomposition_0d(self):
        V, E = cornercases[3]
        V_T, E_T, w = tdlib.PP_MD(V, E)
        S = tdlib.min_coloring_with_treedecomposition(V, E, V_T, E_T)
        self.assertEqual(len(S), 5)

    def test_min_coloring_with_treedecomposition_1(self):
        V_T, E_T, w = tdlib.PP_MD(V_P6, E_P6)
        S = tdlib.min_coloring_with_treedecomposition(V_P6, E_P6, V_T, E_T)
        self.assertEqual(len(S), 2)

    def test_min_coloring_with_treedecomposition_2(self):
        V_T, E_T, w = tdlib.PP_MD(V_K5, E_K5)
        S = tdlib.min_coloring_with_treedecomposition(V_K5, E_K5, V_T, E_T)
        self.assertEqual(len(S), 5)

    def test_min_coloring_with_treedecomposition_3(self):
        V_T, E_T, w = tdlib.PP_MD(V_Petersen, E_Petersen)
        S = tdlib.min_coloring_with_treedecomposition(V_Petersen, E_Petersen, V_T, E_T)
        self.assertEqual(len(S), 3)

    def test_min_coloring_with_treedecomposition_4(self):
        V_T, E_T, w = tdlib.PP_MD(V_Petersen_double, E_Petersen_double)
        S = tdlib.min_coloring_with_treedecomposition(V_Petersen_double, E_Petersen_double, V_T, E_T)
        self.assertEqual(len(S), 3)

    def test_min_coloring_with_treedecomposition_5(self):
        V_T, E_T, w = tdlib.PP_MD(V_Wagner, E_Wagner)
        S = tdlib.min_coloring_with_treedecomposition(V_Wagner, E_Wagner, V_T, E_T)
        self.assertEqual(len(S), 3)

    def test_min_coloring_with_treedecomposition_6(self):
        V_T, E_T, w = tdlib.PP_MD(V_Pappus, E_Pappus)
        S = tdlib.min_coloring_with_treedecomposition(V_Pappus, E_Pappus, V_T, E_T)
        self.assertEqual(len(S), 2)

    def test_min_coloring_with_treedecomposition_7(self):
        V_T, E_T, w = tdlib.PP_MD(V_Grid_5_5, E_Grid_5_5)
        S = tdlib.min_coloring_with_treedecomposition(V_Grid_5_5, E_Grid_5_5, V_T, E_T)
        self.assertEqual(len(S), 2)

    def test_min_coloring_with_treedecomposition_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                V_T, E_T, w = tdlib.PP_MD(V, E)
                S = tdlib.min_coloring_with_treedecomposition(V, E, V_T, E_T)


if __name__ == '__main__':
    unittest.main()
