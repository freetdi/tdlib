
from treedec import graph_bal, pp
from treedec import graph_balu

g = graph_bal([(1,2),(2,3)], 4)
assert(g.num_vertices()==4)
assert(g.num_edges()==2)

if 0:
	g = graph_bal([(1,2),(2,3)], 3)
	for i in g.edges():
		print(i)
	assert(g.num_vertices()==4)

	pp.preprocessing(g)


g = graph_balu([(1,2),(2,3)], 3)
for i in g.edges():
	print(i)
assert(g.num_vertices()==4)

pp.preprocessing(g)
