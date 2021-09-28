# (c) 2021 Felix Salfelder
# GPLv3+

from . import _graph as G
from . import _comb as c

def ppfi(g):
	if isinstance(g, G._balvvd):
		return c._ppfi_balvvd(g)
	elif isinstance(g, G._balvvu):
		return c._ppfi_balvvu(g)
	elif isinstance(g, G._balsvd):
		return c._ppfi_balsvd(g)
	elif isinstance(g, G._balsvu):
		return c._ppfi_balsvu(g)
	else:
		raise ValueError("ppfi: can't handle " + str(g))

def ppfitm(g):
	if isinstance(g, G._balvvd):
		return c._ppfitm_balvvd(g)
	elif isinstance(g, G._balvvu):
		return c._ppfitm_balvvu(g)
	elif isinstance(g, G._balsvd):
		return c._ppfitm_balsvd(g)
	elif isinstance(g, G._balsvu):
		return c._ppfitm_balsvu(g)
	else:
		raise ValueError("ppfitm: can't handle " + str(g))
