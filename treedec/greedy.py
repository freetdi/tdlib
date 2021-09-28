# (c) 2021 Felix Salfelder
# GPLv3+

from . import _graph as G
from . import _greedy as tg

def fi(g):
	if isinstance(g, G._balvvd):
		return tg._fi_balvvd(g)
	elif isinstance(g, G._balvvu):
		return tg._fi_balvvu(g)
	elif isinstance(g, G._balsvd):
		return tg._fi_balsvd(g)
	elif isinstance(g, G._balsvu):
		return tg._fi_balsvu(g)
	else:
		raise ValueError("fi: can't handle " + str(g))

def bmd(g):
	if isinstance(g, G._balvvd):
		return tg._bmd_balvvd(g)
	elif isinstance(g, G._balvvu):
		return tg._bmd_balvvu(g)
	elif isinstance(g, G._balsvd):
		return tg._bmd_balsvd(g)
	elif isinstance(g, G._balsvu):
		return tg._bmd_balsvu(g)
	elif isinstance(g, G._gsgvvu16):
		return tg._bmd_gsgvvu16(g)
	elif isinstance(g, G._gsgvvu32):
		return tg._bmd_gsgvvu32(g)
	elif isinstance(g, G._gsgvvu64):
		return tg._bmd_gsgvvu64(g)
	else:
		raise ValueError("bmd: can't handle " + str(g))
