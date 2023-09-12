// Felix Salfelder, 2017, 2021
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
//

#define BOOST_BIND_GLOBAL_PLACEHOLDERS // for now?
#define BOOST_ALLOW_DEPRECATED_HEADERS // boost/python.hpp?
#include "graph_traits.hpp"
#include "treedec_traits.hpp"
#include <set>
#include <typeinfo>
#include <vector>
#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include "trace.hpp"
#include "graph_iter.hpp"
#include "graph_py.hpp"
#include "boost_graph.hpp"
#include "util_py.hpp"
#include "treedec.hpp"

namespace py=boost::python;

template<class T1, class T2>
struct pair2tuple {
  static PyObject* convert(const std::pair<T1, T2>& pair) {
    return py::incref(py::make_tuple(pair.first, pair.second).ptr());
  }
};

template<class T>
struct set2list {
  static PyObject* convert(const std::set<T>&) {
    return py::incref(py::make_tuple(0, 0).ptr());
    return nullptr;
  }
};
template<class T>
struct set2listp {
  static PyObject* convert(const std::set<T>*&) {
    return py::incref(py::make_tuple(0, 0).ptr());
    return nullptr;
  }
};

namespace detail{

#if 0
unsigned stringtotype(const std::string& type)
{
	unsigned typenumber=graph::b_default;
	if(type==""){
	}else if(type=="badjl_vvd"){
		typenumber = graph::b_badjl_vvd;
	}else{ untested();
		throw exception_invalid("type");
	}
	return typenumber;
}

template<class what>
static boost::shared_ptr<what> graph_factory(
		boost::python::list const& x, boost::python::list const& edgl, std::string const& type)
{ untested();
	typedef boost::python::stl_input_iterator<typename what::vertex_descriptor> iterator_type;
	iterator_type i(x);
	iterator_type e;
	unsigned typenumber=stringtotype(type);

	typedef typename what::vertex_descriptor vd;
	typedef ::detail::pyedge_input_iterator<std::pair<vd, vd> > pi;

	trace1("fact", typenumber);
	return boost::shared_ptr<what>(new what(pi(edgl), pi(), i, e, typenumber));
}
#endif

struct below_checker{
	below_checker(unsigned max=-1u):_m(max){}
	bool operator()(unsigned what) const{return what<_m;}

	unsigned _m;
};

#if 0
template<class what>
static boost::shared_ptr<what> graph_factory_n(
		 unsigned max, boost::python::list const& y, std::string const& type)
{
	trace1("fact_n", max);
	unsigned typenumber=stringtotype(type);

	typedef ::detail::pyedge_input_iterator<std::pair<int, int>, below_checker > pi;

	below_checker C(max);

	trace1("fact", typenumber);
	return boost::shared_ptr<what>(
			new what(pi(y,C), pi(), max, typenumber));
}
#endif

template<class what>
static boost::shared_ptr<what> graph_factory_n0(
		unsigned x, boost::python::list const& y)
{
	return graph_factory_n<what>(x, y, "");
}

template<class what>
static boost::shared_ptr<what> init_boost_graph(size_t n)
{ untested();
	return boost::shared_ptr<what>(new what(n));
}

template<class what>
static boost::shared_ptr<what> init_boost_graph0( boost::python::list const& E, size_t n)
{
//	typedef boost::python::stl_input_iterator<typename what::vertex_descriptor> iterator_type;
//	iterator_type i(x);
//	iterator_type e;

	typedef typename what::vertex_descriptor vd;
	typedef ::detail::pyedge_input_iterator<std::pair<vd, vd> > pi;

	what* ng = new what(pi(E), pi(), n);
	return boost::shared_ptr<what>(ng);
}


#if 0
// want: const reference.
// not possible in python, maybe
template<class P>
typename P::value_type get_first(P const& p){
	return p.first;
}
template<class P>
typename P::value_type get_second(P const& p){
	return p.second;
}
#endif

namespace py=boost::python;
py::object get_first(std::pair<py::object, py::object>& p)
{
	return p.first;
}
namespace py=boost::python;
py::object get_second(std::pair<py::object, py::object>& p)
{
	return p.second;
}

} // detail

template<class T, StringLiteral n>
static std::string treedec_repr(T const& g)
{
//	return typeid(g).name();
	auto v = boost::num_vertices(g);
	return std::string(n.value) + " with " + std::to_string(v) + " bags";
}

template<class T>
//typename treedec::treedec_traits<T>::bag_type const& treedec_get_bag(T& t, unsigned i)
auto const& treedec_get_bag(T& t, unsigned i)
{
	return boost::get(treedec::bag_t(), t, i);
}

template<class T>
//typename treedec::treedec_traits<T>::bag_type const& treedec_get_bag(T& t, unsigned i)
size_t treedec_bagsize(T& t)
{
	return treedec::get_bagsize(t);
}

#define COMMON_TREEDEC_IFACE(gid) \
	py::class_< gid, boost::noncopyable >(#gid, py::init<>()) \
		.def("__init__", py::make_constructor(::detail::init_boost_graph< gid >)) \
		.def("__init__", py::make_constructor(::detail::init_boost_graph0< gid >)) \
		.def(py::init<size_t>()) \
		.def("num_vertices", &boost_num_vertices< gid >) \
		.def("num_edges",    &boost_num_edges< gid >) \
		.def("edges",        &boost_edges< gid >) \
		.def("vertices",     &boost_vertices< gid >) \
		.def("add_vertex",   &boost_add_vertex< gid >) \
		.def("add_edge",     &boost_add_edge< gid >) \
		.def("__repr__",     &treedec_repr<gid, #gid> ) \
		.def("get_bagsize",  &treedec_bagsize<gid>) \
		.def("get_bag",      &treedec_get_bag<gid>, py::return_value_policy<py::copy_const_reference>() )
                                                                          // BUG ^^

template<class G>
size_t boost_num_vertices(G const& g)
{
	return boost::num_vertices(g);
}
template<class G>
size_t boost_num_edges(G const& g)
{
	return boost::num_edges(g);
}

template<class G>
struct edg2pair{
	edg2pair(G const& g) : _g(g)
	{
	}
	typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
    typedef std::pair<vertex_descriptor, vertex_descriptor> edge_descriptor;
	template<class E>
	// std::pair<vertex_descriptor, vertex_descriptor>
	std::pair<unsigned, unsigned>
	operator()(E edg) const {
		auto s=boost::source(edg, _g);
		auto t=boost::target(edg, _g);
		return std::make_pair(s, t);
	}
	//idmap_type const& _m;
	G const& _g;
};

template<class edge_descriptor>
	using edge_iterator_ = any_iterator<const edge_descriptor,
			  std::input_iterator_tag,
			  const edge_descriptor>;

template<class G>
pair_< edge_iterator_< pair_<typename boost::graph_traits<G>::vertex_descriptor> > >
boost_edges(G const& g)
{
	typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
	typedef std::pair<vertex_descriptor, vertex_descriptor> edge_descriptor;
	typedef typename boost::graph_traits<G>::edge_iterator geit_t;
	typedef any_iterator<const edge_descriptor,
			  std::input_iterator_tag,
			  const edge_descriptor> edge_iterator;

	auto x = boost::edges(g);
	geit_t b(x.first);
	geit_t e(x.second);

	// translation
	auto bt=make_transform_iterator(b, edg2pair<G>(g));
	auto et=make_transform_iterator(e, edg2pair<G>(g));
	return std::make_pair(edge_iterator(bt), edge_iterator(et));
	// return std::make_pair(edge_iterator(), edge_iterator());
}

template<class T>
std::string typestring(){
	return "incomplete";
}
template<>
std::string typestring<unsigned>(){
	return "unsigned";
}


template<class T>
struct bag_wrap {
    typedef typename T::value_type V;
    typedef typename T::const_iterator it;
    static V& get(T const& x, int i) {
		 PyErr_SetString( PyExc_IndexError, "TODO");
		 py::throw_error_already_set();

    }

	 static std::pair<it, it> range(T const& s){
		 auto b = s.begin();
		 auto e = s.end();

		 return std::make_pair(b,e);
	 }

//	 static auto pyrange(T const& s){
//		 auto b = s.begin();
//		 auto e = s.end();
//
//		 return new py::range(b,e);
//	 }

	 static std::string repr(T const& t) {
		 auto s = typestring<typename T::value_type>();
		 return "set with " + std::to_string(t.size()) + " elements of type " + s;
	 }
};

BOOST_PYTHON_MODULE(_treedec)
{ itested();
	py::to_python_converter<std::pair<int, bool>, pair2tuple<int, bool> >(); // needed?

//	py::to_python_converter<std::set<const unsigned int>, set2list<const unsigned int> >();
	py::to_python_converter<std::pair<int, int>, pair2tuple<int, int> >();
	py::to_python_converter<std::pair<long, long>, pair2tuple<long, long> >();
	py::to_python_converter<std::pair<unsigned, unsigned>, pair2tuple<unsigned, unsigned> >();
	py::to_python_converter<std::pair<unsigned long, unsigned long>, pair2tuple<unsigned long, unsigned long> >();

// 	std::pair<boost::detail::undirected_edge_iter<std::_List_iterator<boost::list_edge<unsigned long, boost::no_property> >

//    to_python_converter<Vertex, cached_vertex_object<Graph> >();
//    to_python_converter<Edge, cached_edge_object<Edge> >();

//	py::to_python_converter<std::set<unsigned int>, set2list<unsigned int> >();
	// py::to_python_converter<const std::set<unsigned int>, set2list<unsigned int> >();
	// https://wiki.python.org/moin/boost.python/StlContainers ?
	py::class_<std::set<unsigned> >("set_unsigned")
#if 0
		doesn't compile on mac w/ clang++14, boost 1.82
		.def("__iter__",     py::range(&std::set<unsigned>::begin, &std::set<unsigned>::end) )
#else
		// supposedly the same (?) c.f.
		// https://wiki.python.org/moin/boost.python/iterator
		.def("__iter__",     py::iterator<std::set<unsigned>>())
#endif
	//	.def("__iter__",     &set_wrap<std::set<unsigned>>::range )
		.def("__repr__",     &bag_wrap<std::set<unsigned>>::repr );

#if 1
	py::class_<std::vector<uint16_t> >("vec_uint16")
		.def("__iter__",     py::iterator< std::vector< uint16_t>>())
	//	.def("__iter__",     &bag_wrap<std::set<unsigned>>::range )
		.def("__repr__",     &bag_wrap<std::vector<uint16_t>>::repr );

	py::class_<std::vector<uint32_t> >("vec_uint32") \
		.def("__iter__",     py::iterator< std::vector< uint32_t>>())
	//	.def("__iter__",     &bag_wrap<std::set<unsigned>>::range )
		.def("__repr__",     &bag_wrap<std::vector<int32_t>>::repr );

	py::class_<std::vector<uint64_t> >("vec_uint64") \
		.def("__iter__",     py::iterator< std::vector< uint64_t>>())
	//	.def("__iter__",     &bag_wrap<std::set<uint64_t>>::range )
		.def("__repr__",     &bag_wrap<std::vector<uint64_t>>::repr );
#endif

	py::class_<std::pair<py::object, py::object> >("edge")
		    .def("source", &::detail::get_first)
			 .def("target", &::detail::get_second)
			 ;

	// these all seem to use std::set<unsigned> as bags atm.
 	COMMON_TREEDEC_IFACE(_balvvd_treedec);
 	COMMON_TREEDEC_IFACE(_balsvd_treedec);
 	COMMON_TREEDEC_IFACE(_balvvu_treedec);
 	COMMON_TREEDEC_IFACE(_balsvu_treedec);

#ifdef HAVE_GALA_GRAPH_H
 	COMMON_TREEDEC_IFACE(_gsgvvu16_treedec);
 	COMMON_TREEDEC_IFACE(_gsgvvu32_treedec);
 	COMMON_TREEDEC_IFACE(_gsgvvu64_treedec);
#endif
}

//template class wrap::graph<graph_bal, boost::identity_property_map>;
//template class wrap::graph<graph_balu, boost::identity_property_map>;
// template class tf::graph<boost::identity_property_map>;
