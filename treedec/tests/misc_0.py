#!/usr/bin/env python3

from treedec._graph import _balvvu
from treedec._treedec import _balvvu_treedec
from treedec._exact import _ex17_balvvu
from treedec._misc import check_treedec

g = _balvvu([(0,1),(1,2),(2,3),(3,4),(4,0)], 5)
a = _ex17_balvvu(g)
a.do_it()
t = _balvvu_treedec()
a.get_treedec(t)

assert(check_treedec(g, t) == 0)
for i in t.edges():
	print(i)

t.add_edge(0,3)
print("error", check_treedec(g, t))
assert(check_treedec(g, t) == -4)

t.add_vertex()
print("error", check_treedec(g, t))
assert(check_treedec(g, t) == -1)

t = _balvvu_treedec()
print("error", check_treedec(g, t))
assert(check_treedec(g, t) == -5)

