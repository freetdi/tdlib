from treedec._graph import _gsgvvu32

g = _gsgvvu32([(1,2), (3,4)], 5)

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
