import base
import sys
import tdlib
import unittest
from Graph import Graph

V = range(10)
E_ = [ (1,2),
      (1,3),
      (1,4),
      (1,5),
      (1,6),
      (1,7),
      (1,8),
      (1,10),
      (2,3),
      (2,4),
      (2,5),
      (2,6),
      (2,7),
      (2,8),
      (2,9),
      (2,10),
      (3,4),
      (3,5),
      (3,6),
      (3,8),
      (3,9),
      (4,5),
      (4,6),
      (4,7),
      (5,6),
      (5,7),
      (5,8),
      (5,9),
      (5,10),
      (6,9),
      (6,10),
      (7,10),
      (8,9)
]

E = [ (a-1, b-1) for (a,b) in E_ ]
G = Graph(V, E)

what = [ tdlib.exact_decomposition_dynamic,
         tdlib.exact_decomposition_cutset,
         tdlib.exact_decomposition_ex17,
         tdlib.boost_minDegree_decomp,
         tdlib.boost_minDegree_decomp,
         tdlib.fillIn_decomp,
         tdlib.seperator_algorithm,
         tdlib.minDegree_decomp ]

for a in what:
	T, tw = a(G)
	print("decomp")
	print(T.vertices())
	print(T.edges())
	print("width", tw)
	print()
