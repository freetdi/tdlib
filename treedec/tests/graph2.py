#import tdlib
from treedec import graph, subgraph, preprocessing, exception_invalid
from treedec import Graph

# missing test
# g = graph.graph([0,2,2],[(1,2),(0,2)])

x = True
try:
	print("trying bogus edg1")
	g = subgraph([0,1,3],[(1,2),(0,2)])
	print("should not get here1")
	x = False
except exception_invalid as e:
	print("caught...")
	print(e)

try:
	print("trying bogus edg")
	g = graph(3,[(3,0)])
	print("should not get here2")
	x = False
except exception_invalid as e:
	print(e)

assert(x);

G = Graph([0,1, 3,4], [(0,1), (3,4)])
G_, B, lb = preprocessing(G)
print(B, lb)
assert(lb==1)

print("testing 3")

G = Graph(["0","1", "3","4"], [("0","1"), ("3","4")])
G_, B, lb = preprocessing(G)
print(B, lb)
assert(lb==1)

print("=========================")
try:
	G = Graph(["0","1", "3","4"], [("0","1"), ("3","4")])
except exception_invalid as e:
	print(e)
	assert(0)

assert(G.num_vertices()==4)
assert(G.num_edges()==2)
G_, B, lb = preprocessing(G)
B.sort()
print(B, lb)
assert(lb==1)
for b in B:
	b.sort()
assert(B==[['0', '1'], ['3', '4']])

print("=========================")
try:
	# this does not work (never did? or is it python3?)
	G = Graph(["0",1, "3",[1,2,3,4]], [("0",1), ("3",[1,2,3,4])])
	assert(0)
except TypeError as e:
	print(e)

if 0:
	assert(G.num_edges()==2)
	assert(G.num_vertices()==4)
	G_, B, lb = preprocessing(G)
	B.sort()
	print(B, lb)
	assert(lb==1)
	for b in B:
		b.sort()
	print(B)
	assert(B==[[1, '0'], [[1,2,3,4], '3']])


G = subgraph([0,1, 3,4], [(0,1), (3,4)])
G_, B, lb = preprocessing(G)
B.sort()
for b in B:
	b.sort()
print(B, lb)
assert(B==[[0,1], [3,4]])
assert(lb==1)

print("testing vvd")

G = subgraph([0,1, 3,4], [(0,1), (3,4)], "badjl_vvd")
G_, B, lb = preprocessing(G)
B.sort()
for b in B:
	b.sort()
print(B, lb)
assert(B==[[0,1], [3,4]])
assert(lb==1)

print("ipo")
V_Gs_at_ipo = [1,2,3,4,5,6,7,8]
E_Gs_at_ipo = [(1,2),(1,3),(1,4),(2,6),(2,7),(3,6),(3,8),(4,7),(4,8),(5,6),(5,7),(5,8)]

V_GsF__dn_ = [1,2,3,4,5,6,7]
E_GsF__dn_ = [(1,2),(1,5),(1,7),(2,6),(2,7),(3,4),(3,6),(3,7),(4,5),(4,7),(5,6)]
E_GsF__dn_0 = [(0,1),(0,4),(0,6),(1,5),(1,6),(2,3),(2,5),(2,6),(3,4),(3,6),(4,5)]

# ouch. V_Gs_at_ipo starts at 1.
G = subgraph(V_Gs_at_ipo, E_Gs_at_ipo, "")
G_, B, lb = preprocessing(G)
assert(lb == 3)

G = Graph(V_GsF__dn_, E_GsF__dn_, "badjl_vvd")
assert(G.num_edges()==11)
G_, B, lb = preprocessing(G)
assert(lb== 3)

G = graph(len(V_GsF__dn_), E_GsF__dn_0, "badjl_vvd")
assert(G.num_vertices()==7)
assert(G.num_edges()==11)
G_, B, lb = preprocessing(G)
assert(lb== 3)
