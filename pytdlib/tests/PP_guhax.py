import base
import tdlib
import unittest

from graphs import *

class TestTdLib_pre(unittest.TestCase):
    def test_dc_guhax(self):

        vertices=[1,2,3,4,5,6,7,8]
        edges=[(1,2), (1, 3), (1, 4), (1, 5), (2, 4), (2, 6), (2, 7), (3, 5),
                (3, 6), (3, 8), (4, 7), (4, 8), (5, 7),( 5, 8),( 6, 7),( 6, 8),
                (7, 8)]
        G = Graph(vertices, edges)
        lb = tdlib.lower_bound(G, "deltaC_least_c")
        self.assertEqual(lb, 4)
    def test_preprocessing_guhax(self):

        vertices=[1,2,3,4,5,6,7,8]
        edges=[(1,2), (1, 3), (1, 4), (1, 5), (2, 4), (2, 6), (2, 7), (3, 5),
                (3, 6), (3, 8), (4, 7), (4, 8), (5, 7),( 5, 8),( 6, 7),( 6, 8),
                (7, 8)]
        G = Graph(vertices, edges)
        G_, B, lb = tdlib.preprocessing(G)
        self.assertEqual(lb, 4)
        # no rules apply
        self.assertEqual(B, [])

if __name__ == '__main__':
    unittest.main()

# vim:ts=8:sw=4:et
