# (c) 2021 Felix Salfelder
# GPLv3+

from . import _graph as G
from . import _order as o

def treedec_to_ordering(g):
	if isinstance(g, G._balvvd):
		return o._toeo_balvvd(g)
	elif isinstance(g, G._balvvu):
		return o._toeo_balvvu(g)
	elif isinstance(g, G._balsvd):
		return o._toeo_balsvd(g)
	elif isinstance(g, G._balsvu):
		return o._toeo_balsvu(g)
	elif isinstance(g, G._gsgvvu16):
		return o._toeo_gsgvvu16(g)
	elif isinstance(g, G._gsgvvu32):
		return o._toeo_gsgvvu32(g)
	elif isinstance(g, G._gsgvvu64):
		return o._toeo_gsgvvu64(g)
	else:
		raise ValueError("toeo: can't handle " + str(g))
