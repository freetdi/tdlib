
from treedec import _graph as gr
from treedec import _treedec as td
from treedec.exact import *

EU = [(0,1),(1,2),(2,3),(3,0)]

g = gr._balsvu(EU, 4)
print(g)
a = ta(g)
a.do_it()
print(a, a.bagsize())
t = td._balsvu_treedec()
a.get_treedec(t)

j=0
for i in t.vertices():
	print(i)
	B = t.get_bag(j);

	for k in B:
		print("", k, end="")
	print()
	j += 1

g = gr._balvvu(EU, 4)
print(g)
a = ta(g)
a.do_it()
t = td._balvvu_treedec()
a.get_treedec(t)
print(a, "bagsize", a.bagsize(), t)
print(a.get_ordering())

for i in t.edges():
	print(i)

j=0
for i in t.vertices():
	print(i)
	B = t.get_bag(j);

	for k in B:
		print("", k, end="")
	print()
	j += 1
