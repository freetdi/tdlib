import tdlib
import unittest

from graphs import *

class TestTdLib_pre(unittest.TestCase):

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


if __name__ == '__main__':
    unittest.main()
