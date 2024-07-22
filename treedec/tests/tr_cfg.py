import base
import tdlib
import sys

# from graphs import *
import CFGs
import NewCFGs
from Graph import Graph
from treedec.greedy import fi
from treedec.exact import tr
from treedec import _treedec as td

PREFIX = "CFGs"
COUNT = 1817

PREFIX = "NewCFGs"
COUNT = 48
not_connected = { 3,6,10,12,17,20,25,27,35,38,40,41,46 }
too_big = {2,7,8,14,28,24,16,32,30,34,37} # vertices > 100
skip = too_big.union(not_connected)

for i in range(COUNT):
    if i in skip: continue
    G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)), type=1)
    g = G._graph
    print("tr",i, g, end=" ")

    a = tr(g)
    a.do_it()
    t = td._balvvu_treedec()
    a.get_treedec(t)
    assert(tdlib.is_valid_treedecomposition(G, t))

    print("...bagsize", t.get_bagsize(), "lb", a.lower_bound_bagsize())

