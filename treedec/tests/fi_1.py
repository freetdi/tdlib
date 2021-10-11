
from treedec import _graph as gr
from treedec.greedy import fi

E = [(1,2),(2,1),(2,3),(3,2)]
EU = [(1,2),(2,3)]

g = gr._gsgvvu16(EU, 4)
print(g)
a = fi(g)
print(a)
a.do_it()
t = a.get_ordering()
print(t)
# t = a.get_tree_decomposition()

g = gr._gsgvvu32(EU, 4)
print(g)
a = fi(g)
print(a)
a.do_it()
t = a.get_ordering()
print(t)
# t = a.get_tree_decomposition()

g = gr._gsgvvu64(EU, 4)
print(g)
a = fi(g)
print(a)
a.do_it()
t = a.get_ordering()
print(t)
# t = a.get_tree_decomposition()
