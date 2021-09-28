from treedec import lb, exact, comb, greedy, pp
from treedec._misc import check_treedec
from treedec._graph import _balsvu
from treedec._treedec import _balsvu_treedec, _balvvu_treedec

# if gala... try?
from treedec._treedec import _gsgvvu16_treedec, _gsgvvu32_treedec, _gsgvvu64_treedec

def boost_minDegree_decomp(g):
	a = greedy.bmd(g._graph)
	a.do_it()

	t = _gsgvvu16_treedec()
	t = _gsgvvu64_treedec()
	t = _gsgvvu32_treedec()
	t = _balvvu_treedec() # set
	a.get_treedec(t)
	return t, a.bagsize()-1

def exact_decomposition_ex17(g):
	a = exact.ex17(g._graph)
	a.do_it()

	t = _balsvu_treedec()
	a.get_treedec(t)
	return t, a.bagsize()-1

def lower_bound(g, what):
	if(what == "deltaC_least_c"):
		a = lb.deltaC_least_c(g._graph)
	else:
		print("lower_bound incomplete", what)

	a.do_it()
	return a.lower_bound_bagsize()-1

def is_valid_treedecomposition(G, T):
	err = check_treedec(G._graph, T);
	if err:
		print("check", err)
		return False
	else:
		return True



def fillIn_ordering(g):
	incomplete()

def fillIn_decomp(g):
	a = greedy.fi(g._graph)
	a.do_it()

	t = _balsvu_treedec()
	a.get_treedec(t)
	print("treedec numverices", t.num_vertices())
	return t, a.bagsize()-1

def PP_FI_TM(g):
	a = comb.ppfitm(g._graph)
	a.do_it()

	t = _balsvu_treedec()
	a.get_treedec(t)
	return t, a.bagsize()-1

def preprocessing(g):
	a = pp(g._graph)
	a.do_it()

	t = _balsvu_treedec()
	a.get_treedec(t)
	return t, None, a.bagsize()-1

def PP_FI(g):
	a = comb.ppfi(g._graph)
	a.do_it()

	t = _balsvu_treedec()
	a.get_treedec(t)
	return t, a.bagsize()-1

def trivial_decomposition(g):
	incomplete()
def PP_MD(g):
	incomplete()
def seperator_algorithm(g):
	incomplete()
def minDegree_decomp(g):
	incomplete()

# TODO: need more of these to run legacy tests.
