# (c) 2021 Felix Salfelder
# GPLv3+

from . import _graph as G
from . import _exact as ex

def cutset(g):
	if isinstance(g, G._balvvu):
		return ex._cutset_balvvu(g)
	elif isinstance(g, G._balsvu):
		return ex._cutset_balsvu(g)
#	elif isinstance(g, G._balvvd):
#		return ex._cutset_balvvd(g)
#	elif isinstance(g, G._balsvd):
#		return ex._cutset_balsvd(g)
	else:
		raise ValueError("cutset: can't handle " + str(g))

def ta(g):
	if isinstance(g, G._balvvu):
		return ex._ta_balvvu(g)
	elif isinstance(g, G._balsvu):
		return ex._ta_balsvu(g)
#	elif isinstance(g, G._balvvd):
#		return ex._ta_balvvd(g)
#	elif isinstance(g, G._balsvd):
#		return ex._ta_balsvd(g)
	else:
		raise ValueError("ta: can't handle " + str(g))

def ex17(g):
	if isinstance(g, G._balvvu):
		return ex._ex17_balvvu(g)
	elif isinstance(g, G._balsvu):
		return ex._ex17_balsvu(g)
#	elif isinstance(g, G._balvvd):
#		return ex._ex17_balvvd(g)
#	elif isinstance(g, G._balsvd):
#		return ex._ex17_balsvd(g)
	else:
		raise ValueError("ex17: can't handle " + str(g))


def tr(g):
	def size_tag(s):
		upper_bounds = [32,64,128,192,256,512,1024]
		for upper_bound  in upper_bounds:
			if s<= upper_bound:
				return upper_bound
		raise NotImplementedError(f"Handling for more than 1024 vertices is not implemented. Received: {s}.")

	size = g.num_vertices()
	if isinstance(g, G._balvvu):
		alg = "_tr_%i_balvvu"%(size_tag(size))
		return ex.__getattribute__(alg)(g)
	elif isinstance(g, G._balsvu):
		alg = "_tr_%i_balsvu"%(size_tag(size))
		return ex.__getattribute__(alg)(g)
#	elif isinstance(g, G._balvvd):
#		return ex._tr_balvvd(g)
#	elif isinstance(g, G._balsvd):
#		return ex._tr_balsvd(g)
	else:
		raise ValueError("tr: can't handle " + str(g))
