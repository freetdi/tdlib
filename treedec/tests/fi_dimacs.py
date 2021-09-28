import base
import tdlib
import unittest
import sys

# from graphs import *
import Dimacs
from Graph import Graph
from treedec.greedy import fi

V = Dimacs.V_9
E = Dimacs.E_9
V = Dimacs.V_11
E = Dimacs.E_11

G = Graph(V, E, type=1)
T, w = tdlib.PP_FI(G)

assert(tdlib.is_valid_treedecomposition(G, T))

print("done", w)
