
#include "boost_compat.h"
#include <boost/python.hpp>
#include "lower_bounds.hpp"
#include "graph_py.hpp"
#include "algo_py.hpp"

using treedec::lb::impl::deltaC_least_c;
using treedec::algo::default_config;

typedef deltaC_least_c<_balsvu, default_config> _dclc_balsvu;
typedef deltaC_least_c<_balsvd, default_config> _dclc_balsvd;
typedef deltaC_least_c<_balvvu, default_config> _dclc_balvvu;
typedef deltaC_least_c<_balvvd, default_config> _dclc_balvvd;

namespace py = boost::python;

// BUG: copy from exact. move to algo_py
#define xCOMMON_ALGO_IFACE(A, G) \
	py::class_< A ## G, boost::noncopyable>(# A # G, py::no_init) \
			 .def("__init__", py::make_constructor(algo_factory<A ## G, G>)) \
			 .def("__repr__", &make_string<A ## G, # A # G > ) \
			 .def("do_it", static_cast<void ( A ## G:: *)()>(&A ## G ::do_it) )

		    //.def("get_treedec", & A ## G ::get_tree_decomposition< G ## _treedec>)
//        impl::deltaC_least_c<G_t> deltaC_least_c(G);
//        deltaC_least_c.do_it();
//        return (int)deltaC_least_c.lower_bound_bagsize()-1;
			 // .def("bagsize", &A ## G::bagsize)

BOOST_PYTHON_MODULE(_lb)
{ itested();

//	to_python_converter<std::pair<int, bool>, pair2tuple<int, bool> >(); // hack
//	to_python_converter<std::pair<int, int>, pair2tuple<int, int> >();
//	to_python_converter<std::pair<long, long>, pair2tuple<long, long> >();
//	to_python_converter<std::pair<unsigned, unsigned>, pair2tuple<unsigned, unsigned> >();
//	to_python_converter<std::pair<unsigned long, unsigned long>, pair2tuple<unsigned long, unsigned long> >();

	COMMON_ALGO_IFACE(_dclc, _balvvu);
	// COMMON_ALGO_IFACE(_cutset, _balvvd); does not compile.
	COMMON_ALGO_IFACE(_dclc, _balsvu);
	// COMMON_ALGO_IFACE(_cutset, _balsvd); does not compile.

}
