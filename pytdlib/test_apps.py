import tdlib
import unittest
import random

from graphs import *

class TestTdLib_app(unittest.TestCase):
    def test_max_clique_with_treedecomposition_0a(self):
        V, E = cornercases[0]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_clique_with_treedecomposition(G, T)
        self.assertEqual(len(S), 0)

    def test_max_clique_with_treedecomposition_0b(self):
        V, E = cornercases[1]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_clique_with_treedecomposition(G, T)
        self.assertEqual(len(S), 1)

    def test_max_clique_with_treedecomposition_0c(self):
        V, E = cornercases[2]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_clique_with_treedecomposition(G, T)
        self.assertEqual(len(S), 1)

    def test_max_clique_with_treedecomposition_0d(self):
        V, E = cornercases[3]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_clique_with_treedecomposition(G, T)
        self.assertEqual(len(S), 5)

    def test_max_clique_with_treedecomposition_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_clique_with_treedecomposition(G, T)
        self.assertEqual(len(S), 2)

    def test_max_clique_with_treedecomposition_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_clique_with_treedecomposition(G, T)
        self.assertEqual(len(S), 5)

    def test_max_clique_with_treedecomposition_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_clique_with_treedecomposition(G, T)
        self.assertEqual(len(S), 2)

    def test_max_clique_with_treedecomposition_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_clique_with_treedecomposition(G, T)
        self.assertEqual(len(S), 2)

    def test_max_clique_with_treedecomposition_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_clique_with_treedecomposition(G, T)
        self.assertEqual(len(S), 2)

    def test_max_clique_with_treedecomposition_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_clique_with_treedecomposition(G, T)
        self.assertEqual(len(S), 2)

    def test_max_clique_with_treedecomposition_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_clique_with_treedecomposition(G, T)
        self.assertEqual(len(S), 2)

    def test_max_clique_with_treedecomposition_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                G = Graph(V, E)
                T, w = tdlib.PP_MD(G)
                S = tdlib.max_clique_with_treedecomposition(G, T)


    def test_max_independent_set_with_treedecomposition_0a(self):
        V, E = cornercases[0]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_independent_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 0)

    def test_max_independent_set_with_treedecomposition_0b(self):
        V, E = cornercases[1]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_independent_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 1)

    def test_max_independent_set_with_treedecomposition_0c(self):
        V, E = cornercases[2]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_independent_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 5)

    def test_max_independent_set_with_treedecomposition_0d(self):
        V, E = cornercases[3]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_independent_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 1)

    def test_max_independent_set_with_treedecomposition_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_independent_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 3)

    def test_max_independent_set_with_treedecomposition_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_independent_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 1)

    def test_max_independent_set_with_treedecomposition_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_independent_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 4)

    def test_max_independent_set_with_treedecomposition_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_independent_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 8)

    def test_max_independent_set_with_treedecomposition_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_independent_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 3)

    def test_max_independent_set_with_treedecomposition_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_independent_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 9)

    def test_max_independent_set_with_treedecomposition_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.PP_MD(G)
        S = tdlib.max_independent_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 13)

    def test_max_independent_set_with_treedecomposition_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                G = Graph(V, E)
                T, w = tdlib.PP_MD(G)
                S = tdlib.max_independent_set_with_treedecomposition(G, T)


    def test_min_vertex_cover_with_treedecomposition_0a(self):
        V, E = cornercases[0]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_vertex_cover_with_treedecomposition(G, T)
        self.assertEqual(len(S), 0)

    def test_min_vertex_cover_with_treedecomposition_0b(self):
        V, E = cornercases[1]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_vertex_cover_with_treedecomposition(G, T)
        self.assertEqual(len(S), 0)

    def test_min_vertex_cover_with_treedecomposition_0c(self):
        V, E = cornercases[2]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_vertex_cover_with_treedecomposition(G, T)
        self.assertEqual(len(S), 0)

    def test_min_vertex_cover_with_treedecomposition_0d(self):
        V, E = cornercases[3]
        G = Graph(V, E)
        T,  w = tdlib.PP_MD(G)
        S = tdlib.min_vertex_cover_with_treedecomposition(G, T)
        self.assertEqual(len(S), 4)

    def test_min_vertex_cover_with_treedecomposition_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_vertex_cover_with_treedecomposition(G, T)
        self.assertEqual(len(S), 3)

    def test_min_vertex_cover_with_treedecomposition_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_vertex_cover_with_treedecomposition(G, T)
        self.assertEqual(len(S), 4)

    def test_min_vertex_cover_with_treedecomposition_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_vertex_cover_with_treedecomposition(G, T)
        self.assertEqual(len(S), 6)

    def test_min_vertex_cover_with_treedecomposition_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_vertex_cover_with_treedecomposition(G, T)
        self.assertEqual(len(S), 12)

    def test_min_vertex_cover_with_treedecomposition_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_vertex_cover_with_treedecomposition(G, T)
        self.assertEqual(len(S), 5)

    def test_min_vertex_cover_with_treedecomposition_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_vertex_cover_with_treedecomposition(G, T)
        self.assertEqual(len(S), 9)

    def test_min_vertex_cover_with_treedecomposition_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_vertex_cover_with_treedecomposition(G, T)
        self.assertEqual(len(S), 12)

    def test_min_vertex_cover_with_treedecomposition_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                G = Graph(V, E)
                T, w = tdlib.PP_MD(G)
                S = tdlib.min_vertex_cover_with_treedecomposition(G, T)


    def test_min_dominating_set_with_treedecomposition_0a(self):
        V, E = cornercases[0]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_dominating_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 0)

    def test_min_dominating_set_with_treedecomposition_0b(self):
        V, E = cornercases[1]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_dominating_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 1)

    def test_min_dominating_set_with_treedecomposition_0c(self):
        V, E = cornercases[2]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_dominating_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 5)

    def test_min_dominating_set_with_treedecomposition_0d(self):
        V, E = cornercases[3]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_dominating_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 1)

    def test_min_dominating_set_with_treedecomposition_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_dominating_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 2)

    def test_min_dominating_set_with_treedecomposition_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_dominating_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 1)

    def test_min_dominating_set_with_treedecomposition_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_dominating_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 3)

    def test_min_dominating_set_with_treedecomposition_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_dominating_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 6)

    def test_min_dominating_set_with_treedecomposition_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_dominating_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 3)

    def test_min_dominating_set_with_treedecomposition_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_dominating_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 5)

    def test_min_dominating_set_with_treedecomposition_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_dominating_set_with_treedecomposition(G, T)
        self.assertEqual(len(S), 7)

    def test_min_dominating_set_with_treedecomposition_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                G = Graph(V, E)
                T, w = tdlib.PP_MD(G)
                S = tdlib.min_dominating_set_with_treedecomposition(G, T)


    def test_min_coloring_with_treedecomposition_0a(self):
        V, E = cornercases[0]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_coloring_with_treedecomposition(G, T)
        self.assertEqual(len(S), 0)

    def test_min_coloring_with_treedecomposition_0b(self):
        V, E = cornercases[1]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_coloring_with_treedecomposition(G, T)
        self.assertEqual(len(S), 1)

    def test_min_coloring_with_treedecomposition_0c(self):
        V, E = cornercases[2]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_coloring_with_treedecomposition(G, T)
        self.assertEqual(len(S), 1)

    def test_min_coloring_with_treedecomposition_0d(self):
        V, E = cornercases[3]
        G = Graph(V, E)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_coloring_with_treedecomposition(G, T)
        self.assertEqual(len(S), 5)

    def test_min_coloring_with_treedecomposition_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_coloring_with_treedecomposition(G, T)
        self.assertEqual(len(S), 2)

    def test_min_coloring_with_treedecomposition_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_coloring_with_treedecomposition(G, T)
        self.assertEqual(len(S), 5)

    def test_min_coloring_with_treedecomposition_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_coloring_with_treedecomposition(G, T)
        self.assertEqual(len(S), 3)

    def test_min_coloring_with_treedecomposition_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_coloring_with_treedecomposition(G, T)
        self.assertEqual(len(S), 3)

    def test_min_coloring_with_treedecomposition_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_coloring_with_treedecomposition(G, T)
        self.assertEqual(len(S), 3)

    def test_min_coloring_with_treedecomposition_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_coloring_with_treedecomposition(G, T)
        self.assertEqual(len(S), 2)

    def test_min_coloring_with_treedecomposition_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.PP_MD(G)
        S = tdlib.min_coloring_with_treedecomposition(G, T)
        self.assertEqual(len(S), 2)

    def test_min_coloring_with_treedecomposition_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                G = Graph(V, E)
                T, w = tdlib.PP_MD(G)
                S = tdlib.min_coloring_with_treedecomposition(G, T)

if __name__ == '__main__':
    unittest.main()
