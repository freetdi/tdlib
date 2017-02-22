
def is_permutation(x):
	c=[0]*len(x)
	for i in x:
		c[i]=1;
	for i in c:
		if(not i):
			return False
	return True;
