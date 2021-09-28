# from graphs import *
# move to tdlib?
from treedec import graph, subgraph, preprocessing

V_Wagner = range(8)
E_Wagner = [(0, 1),(0, 4), (0, 7), (1, 2), (1, 5), (2, 3), (2, 6), (3, 4), \
            (3, 7), (4, 5), (5, 6), (6, 7)]

V_P6 = range(6)
E_P6 = [(0,1),(1,2),(2,3),(3,4),(4,5)]

V_Wagner = list(V_Wagner) # incomplete
V_P6 = list(V_P6) # incomplete

gg = graph(len(V_Wagner), E_Wagner)
assert(gg.num_vertices()==len(V_Wagner))

g = subgraph(V_Wagner, E_Wagner)
p = subgraph(V_P6, E_P6)

assert(g.num_vertices()==len(V_Wagner))

print(g.num_edges())
print(len(E_Wagner))

assert(g.num_edges()==len(E_Wagner))
preprocessing(g)

# no edges removed
assert(g.num_edges()==len(E_Wagner))

assert(p.num_edges()==5)
preprocessing(p)

# fully preprocessable
assert(p.num_edges()==0)

p = graph(len(V_P6), E_P6)
a, b, c = preprocessing(p)
