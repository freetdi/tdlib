# (c) 2021 Felix Salfelder
# GPLv3+

from . import _graph as G
from . import _lb as lb

def deltaC_least_c(g):
	if isinstance(g, G._balvvu):
		return lb._dclc_balvvu(g)
	elif isinstance(g, G._balsvu):
		return lb._dclc_balsvu(g)
#	elif isinstance(g, G._balvvd):
#		return ex._cutset_balvvd(g)
#	elif isinstance(g, G._balsvd):
#		return ex._cutset_balsvd(g)
	else:
		raise ValueError("deltaC_least_c: can't handle " + str(g))


# TODO: add the others. (macro?)
