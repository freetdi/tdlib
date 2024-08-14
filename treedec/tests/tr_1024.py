import base
import tdlib
import sys
from treedec import _graph as gr

# from graphs import *
from treedec.exact import tr
from treedec import _treedec as td



E_1 = [(i,i+1) for i in range(1023)]

g = gr._balsvu(E_1, 1024)

try:

    print("tr",0, g, end=" ")
    a = tr(g)
    a.do_it()
    t = td._balvvu_treedec()
    a.get_treedec(t)
    #assert(tdlib.is_valid_treedecomposition(G, t))
    print("...bagsize", t.get_bagsize(), "lb", a.lower_bound_bagsize())
except Exception as e:
    print(e)

try:
    g.add_edge(1023,1024)
    print("tr",1, g, end=" ")
    a = tr(g)
    a.do_it()
    t = td._balvvu_treedec()
    a.get_treedec(t)
    #assert(tdlib.is_valid_treedecomposition(G, t))
    print("...bagsize", t.get_bagsize(), "lb", a.lower_bound_bagsize())
except Exception as e:
    print(e)

