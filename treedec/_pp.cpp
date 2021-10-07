// Felix Salfelder, 2017
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2, or (at your option) any
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
#include "preprocessing.hpp"
#include "algo_py.hpp"

using treedec::impl::preprocessing;
//using treedec::algo::default_config;
using treedec::impl::draft::pp_cfg;

typedef preprocessing<_balsvu, pp_cfg> _pp_balsvu;
typedef preprocessing<_balsvd, pp_cfg> _pp_balsvd;
typedef preprocessing<_balvvu, pp_cfg> _pp_balvvu;
typedef preprocessing<_balvvd, pp_cfg> _pp_balvvd;

#ifdef HAVE_GALA_GRAPH_H
typedef preprocessing<_gsgvvu32, pp_cfg> _pp_gsgvvu32;
#endif

namespace detail{

template<class tdt>
struct pp : boost::static_visitor<size_t> {
	pp(tdt& b) : _b(b) { }
	template<class G>
	size_t operator()(G& t) const {
	//	std::vector<std::vector<int> > c_bags; // BUG
	//	incomplete. DO NOT USE
	//	partly implemented in graph_py...
	//	partly implemented in preprocessing.hpp
	//
	//	TODO:
	//	- use pp class directly
	//	- produce something useful *here*
#if 1
		int c_lb=-1;
		treedec::preprocessing(t, _b, c_lb);
		return c_lb+1;
#else
		impl::preprocessing<G> a(t);
		a.do_it();
#endif

	}
	tdt& _b;
};

}

// YUCK.
typedef std::vector< boost::tuple< unsigned,
		  std::set<unsigned> > > pp_data_type;


namespace treedec{
namespace frontend{

#if 0
template<class A>
template<class tdt>
size_t graph<A>::preprocessing(tdt& td_bags)
{

#if 0 // does not work. carful.
	incomplete();
//	put proper implementation here;
	return 0;
}

template<class A>
template<>
size_t anygraph_<A>::preprocessing(pp_data_type& td_bags)
{ untested();

#endif
	//	incomplete. DO NOT USE
	//	partly implemented in graph_py...
	//	partly implemented in preprocessing.hpp
	auto X=::detail::pp<tdt>(td_bags);
	return boost::apply_visitor(X, _g);

//	int c_lb=-1;
///	std::vector<std::vector<int> > c_bags; // BUG
//	treedec::preprocessing(_g, c_bags, c_lb);
//
//	return X._lb_bs;
}

template
size_t graph<std::vector<unsigned> >::preprocessing(pp_data_type&);
template
size_t graph<boost::typed_identity_property_map<unsigned long> >::preprocessing(pp_data_type&);
template
size_t Graph::preprocessing(pp_data_type&);

#endif
} // frontend

} //treedec

namespace wrap{


template<class graph>
//boost::python::object
int preprocessing(graph& g)
{
	 pp_data_type b;
   // size_t lb_bs = g.preprocessing(b); // pass pythonstuff?
	// boost::python::list a, ret;
#if 1
		int c_lb=-1;
		treedec::preprocessing(g, b, c_lb);
		return c_lb+1;
#else
		impl::preprocessing<G> a(g);
		a.do_it();

	// tdt& _b;

	 /// YUCK
	 // pass something resonable to PP
	 // WRONG PLACE. has nothing to do with python.
	 for(auto const& i: b){
		 boost::python::list c;
		 unsigned j = boost::get<0>(i);
		 c.append(g.maphack(j)); // this is the eliminated node.
		                         // bad/FIXME: the order is lost.

		 auto const& bagi=boost::get<1>(i);
		 for(auto const& j : bagi){
			 c.append(g.maphack(j));
		 }

		 a.append(c);
	 }

	 ret.append(g); // FIXME (remove)
	 ret.append(a);
	 ret.append(int(lb_bs)-1); // yuck
	 return ret;
#endif
}
} // wrap


BOOST_PYTHON_MODULE(_pp)
{ itested();
	boost::python::numpy::initialize();
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
	COMMON_ALGO_IFACE(_pp, _balvvu);
	COMMON_ALGO_IFACE(_pp, _balvvd);
	COMMON_ALGO_IFACE(_pp, _balsvu);
	COMMON_ALGO_IFACE(_pp, _balsvd);

#ifdef HAVE_GALA_GRAPH_H
	COMMON_ALGO_IFACE(_pp, _gsgvvu32)
		    .def("get_treedec", static_cast<void ( _pp_gsgvvu32:: * )(_gsgvvu32_treedec&)> \
					     ( &_pp_gsgvvu32::get_tree_decomposition< _gsgvvu32_treedec>) );
#endif


	// OLD
	def("preprocessing", wrap::preprocessing<_balvvu>);
	def("preprocessing", wrap::preprocessing<_balvvd>);
	def("preprocessing", wrap::preprocessing<_balsvu>);
	def("preprocessing", wrap::preprocessing<_balsvd>);

}

// template int wrap::preprocessing<graph_balu>(graph_balu&);
// template int wrap::preprocessing<graph_bal>(graph_bal&);

// BUG
template size_t wrap::graph<_balsvu>::preprocessing(pp_data_type&);
template size_t wrap::graph<_balsvd>::preprocessing(pp_data_type&);
template size_t wrap::graph<_balvvu>::preprocessing(pp_data_type&);
template size_t wrap::graph<_balvvd>::preprocessing(pp_data_type&);

