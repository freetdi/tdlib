import tdlib
import unittest

from graphs import *

class TestTdLib(unittest.TestCase):
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

if __name__ == '__main__':
    unittest.main()
