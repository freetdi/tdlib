// Felix Salfelder, 2021
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option) any
// later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.


#include "boost_compat.h"
#include "graph_py.hpp"
#include "boost_graph.hpp"
#include <vector>
#include <set>
#include "algo_py.hpp"
#include "elimination_orderings.hpp"

using treedec::to_elimination_ordering;

typedef to_elimination_ordering<_balsvu> _toeo_balsvu;
typedef to_elimination_ordering<_balsvd> _toeo_balsvd;
typedef to_elimination_ordering<_balvvu> _toeo_balvvu;
typedef to_elimination_ordering<_balvvd> _toeo_balvvd;

#ifdef HAVE_GALA_GRAPH_H
typedef to_elimination_ordering<_gsgvvu32> _toeo_gsgvvu32;
#endif

BOOST_PYTHON_MODULE(_order)
{ itested();
	namespace py=boost::python;
	using namespace boost::python;

//	to_python_converter<std::pair<int, bool>, pair2tuple<int, bool> >(); // hack
//	to_python_converter<std::pair<int, int>, pair2tuple<int, int> >();
//	to_python_converter<std::pair<long, long>, pair2tuple<long, long> >();
//	to_python_converter<std::pair<unsigned, unsigned>, pair2tuple<unsigned, unsigned> >();
//	to_python_converter<std::pair<unsigned long, unsigned long>, pair2tuple<unsigned long, unsigned long> >();

//	class_<std::pair<py::object, py::object> >("edge")
//		    .def("source", &::detail::get_first)
//			 .def("target", &::detail::get_second)
//			 ;
	COMMON_ALGO_IFACE(_toeo, _balvvu)
		    .def("set_treedec", static_cast<void ( _toeo_balvvu:: * )(_balvvu_treedec const*)> \
					     ( &_toeo_balvvu::set_tree_decomposition<_balvvu_treedec>) );
	COMMON_ALGO_IFACE(_toeo, _balvvd);
	COMMON_ALGO_IFACE(_toeo, _balsvu);
	COMMON_ALGO_IFACE(_toeo, _balsvd);

#ifdef HAVE_GALA_GRAPH_H
	COMMON_ALGO_IFACE(_toeo, _gsgvvu32)
		    .def("set_treedec", static_cast<void ( _toeo_gsgvvu32:: * )(_gsgvvu32_treedec const*)> \
					     ( &_toeo_gsgvvu32::set_tree_decomposition<_gsgvvu32_treedec>) );
#endif

}

// ..?
// template size_t wrap::graph<_balsvu>::preprocessing(pp_data_type&);
// template size_t wrap::graph<_balsvd>::preprocessing(pp_data_type&);
// template size_t wrap::graph<_balvvu>::preprocessing(pp_data_type&);
// template size_t wrap::graph<_balvvd>::preprocessing(pp_data_type&);
