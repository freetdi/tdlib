from treedec._graph import _balvvu, _balsvd, _balsvu, _balvvd

g = _balvvd([(1,2), (3,4)], 5)

print(g)

e = g.edges();
print("edges")
# print(e) # incomplete()
for i in e:
	print(i)

v = g.vertices();
for i in v:
	print(i)

iv = iter(v)
a = next(iv);
assert(a==0)
a = next(iv);
assert(a==1)

c=0
for i in v:
	assert(i==c)
	c += 1

assert(c==5)

h = _balsvu([(1,2),(0,2)], 3)
c = 0
for i in h.edges():
	print(i)
	c += 1
assert(c==2)
assert(h.num_edges() == 2)
print(h.add_edge(1, 2))
assert(h.num_edges() == 2)
print(h.add_edge(1, 0))
assert(h.num_edges() == 3)

h = _balvvu([(1,2),(0,2)], 8)
assert(h.num_vertices()==8)
assert(h.num_edges()==2)

c = 0
A = h.edges()

for i in h.edges():
	c += 1
	print("edg", c, i)
assert(c==2)

h.add_edge(1,2)

g = _balvvd([(1,2),(2,3)], 4)
assert(g.num_vertices()==4)
assert(g.num_edges()==2)

g = _balvvd([(1,2),(2,3),(2,1),(3,2)], 3)
for i in g.edges():
	print(i)
assert(g.num_vertices()==4)
assert(g.num_edges()==4)
