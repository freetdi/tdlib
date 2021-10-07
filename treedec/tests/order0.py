#!/usr/bin/env python3

import treedec.order as tdo
import treedec._graph as gr
from treedec import _treedec as td
from treedec.exact import *

EU = [(0,1),(1,2),(2,3),(3,0)]

g = gr._balvvu(EU, 4)
print(g)
a = ta(g)
a.do_it()
t = td._balvvu_treedec()
a.get_treedec(t)
print(a, "bagsize", a.bagsize(), t)

o = tdo.treedec_to_ordering(g)
o.set_treedec(t)

ord = o.get_ordering()
print(ord)

