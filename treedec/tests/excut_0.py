
from treedec import _graph as gr
from treedec.exact import *

E = [(1,2),(2,1),(2,3),(3,2)]
EU = [(0,1),(1,2),(2,3),(3,0)] # BUG: paths don't work

g = gr._balsvu(EU, 4)
print(g)
a = cutset(g)
a.do_it()
print(a, a.bagsize())
# t = a.get_treedec()

g = gr._balvvu(EU, 4)
print(g)
a = cutset(g)
a.do_it()
print(a, a.bagsize())
# t = a.get_treedec()

