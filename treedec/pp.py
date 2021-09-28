# (c) 2021 Felix Salfelder
# GPLv3+

from . import _graph as G
from . import _pp as c

def pp(g):
	if isinstance(g, G._balvvd):
		return c._pp_balvvd(g)
	elif isinstance(g, G._balvvu):
		return c._pp_balvvu(g)
	elif isinstance(g, G._balsvd):
		return c._pp_balsvd(g)
	elif isinstance(g, G._balsvu):
		return c._pp_balsvu(g)
	elif isinstance(g, G._gsgvvu16):
		return c._pp_gsgvvu16(g)
	elif isinstance(g, G._gsgvvu32):
		return c._pp_gsgvvu32(g)
	elif isinstance(g, G._gsgvvu64):
		return c._pp_gsgvvu64(g)
	else:
		raise ValueError("pp: can't handle " + str(g))
