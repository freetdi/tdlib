import tdlib
import unittest

from graphs import *

class TestTdLib_pre(unittest.TestCase):
    def test_preprocessing_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            G_, B, lb = tdlib.preprocessing(G)
            self.assertEqual(G_.vertices(), [])
            self.assertEqual(G_.edges(), [])

    def test_preprocessing_1(self):
        G = Graph(V_P6, E_P6)
        G_, B, lb = tdlib.preprocessing(G)
        self.assertEqual(G_.vertices(), [])
        self.assertEqual(G_.edges(), [])
        for i in range(len(B)):
            B[i].sort()
        B.sort()
        self.assertEqual(B, [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5]])
        self.assertEqual(lb, 1)

    def test_preprocessing_2(self):
        G = Graph(V_K5, E_K5)
        G_, B, lb = tdlib.preprocessing(G)
        self.assertEqual(G_.vertices(), [])
        self.assertEqual(G_.edges(), [])
        for i in range(len(B)):
            B[i].sort()
        B.sort()
        self.assertEqual(B, [[0, 1, 2, 3, 4], [1, 2, 3, 4], [2, 3, 4], [3, 4]])
        self.assertEqual(lb, 4)


    def test_preprocessing_3(self):
        G = Graph(V_Petersen, E_Petersen)
        G_, B, lb = tdlib.preprocessing(G)
        self.assertEqual(G_.vertices(), V_Petersen)
        self.assertEqual(G_.edges(), [0,1,0,4,0,5,1,2,1,6,2,3,2,7,3, \
                              4,3,8,4,9,5,7,5,8,6,8,6,9,7,9])

        self.assertEqual(len(B), 0)
        self.assertEqual(lb, 4)

    def test_preprocessing_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        G_, B, lb = tdlib.preprocessing(G)
        self.assertEqual(len(B), 0)

    def test_preprocessing_5(self):
        G = Graph(V_Wagner, E_Wagner)
        G_, B, lb = tdlib.preprocessing(G)
        self.assertEqual(len(B), 0)

    def test_preprocessing_6(self):
        G = Graph(V_Pappus, E_Pappus)
        G_, B, lb = tdlib.preprocessing(G)
        self.assertEqual(len(B), 0)

    def test_preprocessing_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        G_, B, lb = tdlib.preprocessing(G)
        self.assertEqual(lb, 4)

    def test_preprocessing_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                G = Graph(V, E)
                G_, B, lb = tdlib.preprocessing(G)


    def test_PP_MD_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            T, w = tdlib.PP_MD(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_PP_MD_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.PP_MD(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 1)

    def test_PP_MD_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.PP_MD(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_PP_MD_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.PP_MD(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w >= 4 and w <= 5, True) #c++11 issuse

    def test_PP_MD_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.PP_MD(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w >= 4 and w <= 5, True) #c++11 issuse

    def test_PP_MD_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.PP_MD(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_PP_MD_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.PP_MD(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 6)

    def test_PP_MD_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.PP_MD(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 5)

    def test_PP_MD_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                G = Graph(V, E)
                T, w = tdlib.PP_MD(G)
                self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)


    def test_PP_FI_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            T, w = tdlib.PP_FI(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_PP_FI_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.PP_FI(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 1)

    def test_PP_FI_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.PP_FI(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_PP_FI_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.PP_FI(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w >= 4 and w <= 5, True) #c++11 issuse

    def test_PP_FI_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.PP_FI(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w >= 4 and w <= 5, True) #c++11 issuse

    def test_PP_FI_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.PP_FI(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_PP_FI_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.PP_FI(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 6)

    def test_PP_FI_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.PP_FI(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 5)

    def test_PP_FI_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                G = Graph(V, E)
                T, w = tdlib.PP_FI(G)
                self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)


    def test_PP_FI_TM_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            T, w = tdlib.PP_FI_TM(G)
            self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

    def test_PP_FI_TM_1(self):
        G = Graph(V_P6, E_P6)
        T, w = tdlib.PP_FI_TM(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 1)

    def test_PP_FI_TM_2(self):
        G = Graph(V_K5, E_K5)
        T, w = tdlib.PP_FI_TM(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_PP_FI_TM_3(self):
        G = Graph(V_Petersen, E_Petersen)
        T, w = tdlib.PP_FI_TM(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w >= 4 and w <= 5, True) #c++11 issuse

    def test_PP_FI_TM_4(self):
        G = Graph(V_Petersen_double, E_Petersen_double)
        T, w = tdlib.PP_FI_TM(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w >= 4 and w <= 5, True) #c++11 issuse

    def test_PP_FI_TM_5(self):
        G = Graph(V_Wagner, E_Wagner)
        T, w = tdlib.PP_FI_TM(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 4)

    def test_PP_FI_TM_6(self):
        G = Graph(V_Pappus, E_Pappus)
        T, w = tdlib.PP_FI_TM(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 6)

    def test_PP_FI_TM_7(self):
        G = Graph(V_Grid_5_5, E_Grid_5_5)
        T, w = tdlib.PP_FI_TM(G)
        self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)
        self.assertEqual(w, 5)

    def test_PP_FI_TM_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                G = Graph(V, E)
                T, w = tdlib.PP_FI_TM(G)
                self.assertEqual(tdlib.is_valid_treedecomposition(G, T), True)

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
