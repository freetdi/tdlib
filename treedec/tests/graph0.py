
from treedec import graph, Graph, subgraph

g = graph(17)

v = g.vertices();

iv = iter(v)
a = next(iv);
assert(a==0)
a = next(iv);
assert(a==1)

c=0
for i in v:
	c += 1
	print(i)

assert(c==17)

G = Graph(["Hello","World", "3.1","4"], [("3.1","4"), ("Hello", "4")])
v = G.vertices();

iv = iter(v)
a = next(iv);
assert(a=='Hello')
a = next(iv);
assert(a=='World')

c = 0
for i in v:
	c += 1
	print(i)

assert(c==4)

e = G.edges()
c = 0
for i in e:
	c += 1
	print(i.source())
	print(i.target())
	print(i) #incomplete.

assert(c==2)

print(g.backend_typename())
assert(g.num_vertices()==17)
assert(g.num_edges()==0)

h = subgraph([0,1,2],[(1,2),(0,2)])
assert(h.num_vertices()==3)
assert(h.num_edges()==2)

g = h #?!
assert(g.num_vertices()==3)
assert(g.num_edges()==2)

#does not work
# import copy
# g = copy.copy(h)
# assert(g.num_vertices()==3)
# assert(g.num_edges()==2)
