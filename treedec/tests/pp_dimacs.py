import base
import tdlib
import sys

# from graphs import *
import Dimacs
from Graph import Graph
from treedec.greedy import fi
from treedec.pp import pp
from treedec import _treedec as td

PREFIX = "Dimacs"
COUNT = 81

for i in range(COUNT):
	G = Graph(eval(PREFIX+".V_"+str(i)), eval(PREFIX+".E_"+str(i)), type=1)
	g = G._graph

	a = pp(g)
	a.do_it()
	t = td._balvvu_treedec()
	a.get_treedec(t)
	assert(tdlib.is_valid_treedecomposition(G, t))

	print("pp",i, "bagsize", t.get_bagsize(), "lb", a.lower_bound_bagsize())
