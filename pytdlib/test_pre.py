import tdlib
import unittest

from graphs import *

class TestTdLib_pre(unittest.TestCase):
    """ segfault
    def test_preprocessing_0(self):
        for V, E in cornercases:
            V_, E_, B, lb = tdlib.preprocessing(V, E)
            self.assertEqual(V_, [])
            self.assertEqual(E_, [])
    """

    def test_preprocessing_1(self):
        V_, E_, B, lb = tdlib.preprocessing(V_P6, E_P6)
        self.assertEqual(V_, [])
        self.assertEqual(E_, [])
        for i in range(len(B)):
            B[i].sort()
        B.sort()
        self.assertEqual(B, [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5]])
        self.assertEqual(lb, 1)

    def test_preprocessing_2(self):
        V_, E_, B, lb = tdlib.preprocessing(V_K5, E_K5)
        self.assertEqual(V_, [])
        self.assertEqual(E_, [])
        for i in range(len(B)):
            B[i].sort()
        B.sort()
        self.assertEqual(B, [[0, 1, 2, 3, 4], [1, 2, 3, 4], [2, 3, 4], [3, 4]])
        self.assertEqual(lb, 4)


    def test_preprocessing_3(self):
        V_, E_, B, lb = tdlib.preprocessing(V_Petersen, E_Petersen)
        self.assertEqual(V_, V_Petersen)
        self.assertEqual(E_, [0,1,0,4,0,5,1,2,1,6,2,3,2,7,3, \
                              4,3,8,4,9,5,7,5,8,6,8,6,9,7,9])

        self.assertEqual(len(B), 0)
        self.assertEqual(lb, 4)

    def test_preprocessing_4(self):
        V_, E_, B, lb = tdlib.preprocessing(V_Petersen_double, E_Petersen_double)
        self.assertEqual(len(B), 0)

    def test_preprocessing_5(self):
        V_, E_, B, lb = tdlib.preprocessing(V_Wagner, E_Wagner)
        self.assertEqual(len(B), 0)

    def test_preprocessing_6(self):
        V_, E_, B, lb = tdlib.preprocessing(V_Pappus, E_Pappus)
        self.assertEqual(len(B), 0)

    def test_preprocessing_7(self):
        V_, E_, B, lb = tdlib.preprocessing(V_Grid_5_5, E_Grid_5_5)
        self.assertEqual(lb, 4)

    """ segfault
    def test_preprocessing_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                V_T, E_T, B, lb = tdlib.preprocessing(V, E)
    """

    def test_PP_MD_0(self):
        for V, E in cornercases:
            V_T, E_T, w = tdlib.PP_MD(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    def test_PP_MD_1(self):
        V, E, w = tdlib.PP_MD(V_P6, E_P6)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E), True)
        self.assertEqual(w, 1)

    def test_PP_MD_2(self):
        V, E, w = tdlib.PP_MD(V_K5, E_K5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V, E), True)
        self.assertEqual(w, 4)

    def test_PP_MD_3(self):
        V, E, w = tdlib.PP_MD(V_Petersen, E_Petersen)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E), True)
        self.assertEqual(w >= 4 and w <= 5, True) #c++11 issuse

    def test_PP_MD_4(self):
        V, E, w = tdlib.PP_MD(V_Petersen_double, E_Petersen_double)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen_double, E_Petersen_double, V, E), True)
        self.assertEqual(w >= 4 and w <= 5, True) #c++11 issuse

    def test_PP_MD_5(self):
        V, E, w = tdlib.PP_MD(V_Wagner, E_Wagner)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V, E), True)
        self.assertEqual(w, 4)

    def test_PP_MD_6(self):
        V, E, w = tdlib.PP_MD(V_Pappus, E_Pappus)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V, E), True)
        self.assertEqual(w, 6)

    def test_PP_MD_7(self):
        V, E, w = tdlib.PP_MD(V_Grid_5_5, E_Grid_5_5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Grid_5_5, E_Grid_5_5, V, E), True)
        self.assertEqual(w, 5)

    def test_PP_MD_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                V_T, E_T, w = tdlib.PP_MD(V, E)
                self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    def test_PP_FI_0(self):
        for V, E in cornercases:
            V_T, E_T, w = tdlib.PP_FI(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    def test_PP_FI_1(self):
        V, E, w = tdlib.PP_FI(V_P6, E_P6)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E), True)
        self.assertEqual(w, 1)

    def test_PP_FI_2(self):
        V, E, w = tdlib.PP_FI(V_K5, E_K5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V, E), True)
        self.assertEqual(w, 4)

    def test_PP_FI_3(self):
        V, E, w = tdlib.PP_FI(V_Petersen, E_Petersen)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E), True)
        self.assertEqual(w >= 4 and w <= 5, True) #c++11 issuse

    def test_PP_FI_4(self):
        V, E, w = tdlib.PP_FI(V_Petersen_double, E_Petersen_double)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen_double, E_Petersen_double, V, E), True)
        self.assertEqual(w >= 4 and w <= 5, True) #c++11 issuse

    def test_PP_FI_5(self):
        V, E, w = tdlib.PP_FI(V_Wagner, E_Wagner)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V, E), True)
        self.assertEqual(w, 4)

    def test_PP_FI_6(self):
        V, E, w = tdlib.PP_FI(V_Pappus, E_Pappus)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V, E), True)
        self.assertEqual(w, 6)

    def test_PP_FI_7(self):
        V, E, w = tdlib.PP_FI(V_Grid_5_5, E_Grid_5_5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Grid_5_5, E_Grid_5_5, V, E), True)
        self.assertEqual(w, 5)

    def test_PP_FI_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                V_T, E_T, w = tdlib.PP_FI(V, E)
                self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    def test_PP_FI_TM_0(self):
        for V, E in cornercases:
            V_T, E_T, w = tdlib.PP_FI_TM(V, E)
            self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

    def test_PP_FI_TM_1(self):
        V, E, w = tdlib.PP_FI_TM(V_P6, E_P6)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_P6, E_P6, V, E), True)
        self.assertEqual(w, 1)

    def test_PP_FI_TM_2(self):
        V, E, w = tdlib.PP_FI_TM(V_K5, E_K5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_K5, E_K5, V, E), True)
        self.assertEqual(w, 4)

    def test_PP_FI_TM_3(self):
        V, E, w = tdlib.PP_FI_TM(V_Petersen, E_Petersen)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen, E_Petersen, V, E), True)
        self.assertEqual(w >= 4 and w <= 5, True) #c++11 issuse

    def test_PP_FI_TM_4(self):
        V, E, w = tdlib.PP_FI_TM(V_Petersen_double, E_Petersen_double)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Petersen_double, E_Petersen_double, V, E), True)
        self.assertEqual(w >= 4 and w <= 5, True) #c++11 issuse

    def test_PP_FI_TM_5(self):
        V, E, w = tdlib.PP_FI_TM(V_Wagner, E_Wagner)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Wagner, E_Wagner, V, E), True)
        self.assertEqual(w, 4)

    def test_PP_FI_TM_6(self):
        V, E, w = tdlib.PP_FI_TM(V_Pappus, E_Pappus)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Pappus, E_Pappus, V, E), True)
        self.assertEqual(w, 6)

    def test_PP_FI_TM_7(self):
        V, E, w = tdlib.PP_FI_TM(V_Grid_5_5, E_Grid_5_5)
        self.assertEqual(tdlib.is_valid_treedecomposition(V_Grid_5_5, E_Grid_5_5, V, E), True)
        self.assertEqual(w, 5)

    def test_PP_FI_TM_8(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                V_T, E_T, w = tdlib.PP_FI_TM(V, E)
                self.assertEqual(tdlib.is_valid_treedecomposition(V, E, V_T, E_T), True)

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
