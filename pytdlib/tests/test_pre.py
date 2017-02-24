import base
import sys
import tdlib
import unittest

from graphs import *

#don't confuse python unittest
sys.argv=sys.argv[:1]

class TestTdLib_pre(unittest.TestCase):
    def test_preprocessing_0(self):
        for V, E in cornercases:
            G = Graph(V, E)
            G_, B, lb = tdlib.preprocessing(G)
            self.assertEqual(G_.vertices(), [])
            self.assertEqual(G_.edges(), [])

    def test_preprocessing_P6(self):
        G = Graph(V_P6, E_P6)
        G_, B, lb = tdlib.preprocessing(G)
        self.assertEqual(G_.vertices(), [])
        self.assertEqual(G_.edges(), [])
        for i in range(len(B)):
            B[i].sort()
        B.sort()
        self.assertEqual(B, [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5]])
        self.assertEqual(lb, 1)

    def test_preprocessing_K5(self):
        G = Graph(V_K5, E_K5)
        G_, B, lb = tdlib.preprocessing(G)
        self.assertEqual(G_.vertices(), [])
        self.assertEqual(G_.edges(), [])
        for i in range(len(B)):
            B[i].sort()
        B.sort()
        # self.assertEqual(B, [[0, 1, 2, 3, 4], [1, 2, 3, 4], [2, 3, 4], [3, 4]])
        self.assertEqual(lb, 4)


    def test_preprocessing_Peter(self):
        G = Graph(V_Petersen, E_Petersen)
        G_, B, lb = tdlib.preprocessing(G)
        self.assertEqual(G_.vertices(), V_Petersen)
#        self.assertEqual(G_.edges(), [0,1,0,4,0,5,1,2,1,6,2,3,2,7,3,
#                              4,3,8,4,9,5,7,5,8,6,8,6,9,7,9])

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
        G = Graph(V_Gs_at_ipo, E_Gs_at_ipo)
        G_, B, lb = tdlib.preprocessing(G)
        self.assertEqual(lb, 3)

    def test_preprocessing_GNP(self):
        for n in range(0, 13):
            for i in range(0, 10):
                V, E = randomGNP(n, 0.2)
                G = Graph(V, E)
                G_, B, lb = tdlib.preprocessing(G)

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
