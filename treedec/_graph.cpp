// Felix Salfelder, 2017, 2021
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
//

#include "boost_compat.h"

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


// unneeded. probably
// typedef wrap::graph<graph_bal, boost::identity_property_map > graph_bal_wrap;
// typedef wrap::graph<graph_balu, boost::identity_property_map > graph_balu_wrap;

using namespace boost::python; // BUG
namespace py = boost::python;

typedef std::vector< boost::tuple< unsigned,
		  std::set<unsigned> > > pp_data_type;

template<class T1, class T2>
struct pair2tuple {
  static PyObject* convert(const std::pair<T1, T2>& pair) {
    return incref(py::make_tuple(pair.first, pair.second).ptr());
  }
};

template<class T>
struct addedgeret {
  static PyObject* convert( std::pair< std::pair<T, T>, bool> const& r) {
    auto p = py::make_tuple(r.first.first, r.first.second);

    return incref(py::make_tuple(p, r.second).ptr());
  }
};

#if 0
struct makeunsignedlong {
  static PyObject* convert(unsigned long const & x) {
    return incref(new PyObject()); // py::make_integer(x));
  }
};
struct makeunsigned {
  static PyObject* convert(unsigned const & x) {
    return incref(new PyObject()); // py::make_integer(x));
  }
};
struct makelong {
  static PyObject* convert(long const & x) {
    return incref(new PyObject()); // py::make_integer(x));
  }
};
#endif

// BUG? detail..
static PyObject *pythonExceptionType = NULL;
static PyObject* createExceptionClass(const char* name, PyObject* baseTypeObj = PyExc_Exception)
{
    using std::string;
    namespace bp = boost::python;

    const string scopeName = bp::extract<string>(bp::scope().attr("__name__"));
    const string qualifiedName0 = scopeName + "." + name;
    PyObject* typeObj = PyErr_NewException(qualifiedName0.c_str(), baseTypeObj, 0);
    if (!typeObj) bp::throw_error_already_set();
    bp::scope().attr(name) = bp::handle<>(bp::borrowed(typeObj));
    return typeObj;
}

namespace detail{

template<class what>
static boost::shared_ptr<what> graph_factory0(
		boost::python::list const& x, boost::python::list const& y)
{
	return graph_factory<what>(x, y, "");
}

// PyObject* ex_invalid_pytype=NULL;
//
// void trans_ex(exception_invalid const &e)
// { itested();
// 	assert(ex_invalid_pytype!=NULL);
// 	boost::python::object pythonExceptionInstance(e);
// 	PyErr_SetObject(ex_invalid_pytype, pythonExceptionInstance.ptr());
// 	// PyErr_SetString(PyExc_RuntimeError, e.what());
// }

static void translate_exception_invalid(exception_invalid const &e)
{
    using namespace boost;
    python::object exc_t(python::handle<>(python::borrowed(pythonExceptionType)));
    exc_t.attr("cause") = python::object(e); // add the wrapped exception to the Python exception
//    exc_t.attr("what") = python::object(e.what()); // for convenience
    PyErr_SetString(pythonExceptionType, e.what()); // the string is used by print(exception) in python
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

	typedef typename boost::graph_traits<what>::vertex_descriptor vd;
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

#define COMMON_GRAPH_PY_IFACE(gid) \
		.def(init<size_t>()) \
		.def(init<std::string const&>()) \
		.def(init<size_t, std::string const&>()) \
		.def("backend_typename", & gid ::backend_typename) \
		.def("num_vertices",     & gid ::num_vertices) \
		.def("num_edges",        & gid ::num_edges) \
		.def("add_edge",         & gid ::add_edge) \
		.def("vertices",         & gid ::vertices) \
		.def("edges",            & gid ::edges) \
		; \
	typedef std::pair<gid ::vertex_iterator, gid ::vertex_iterator> gid ## _vip; \
	class_< gid ## _vip >( #gid "_vertex_range", no_init ) \
		 .def("__iter__" , py::range(& gid ## _vip::first, & gid ## _vip::second)) \
		; \
	typedef std::pair<gid ::edge_iterator, gid ::edge_iterator> gid ## _eip; \
	class_< gid ## _eip >( #gid "_edge_range", no_init ) \
		 .def("__iter__" , py::range(& gid ## _eip::first, & gid ## _eip::second)) \
		 ;

template<class T, StringLiteral n>
static std::string make_string(T const& g)
{
//	return typeid(g).name();
	auto v = boost::num_vertices(g);
	auto e = boost::num_edges(g);
	return std::string(n.value) + " on " + std::to_string(v) + " vertices "
		+ "with " + std::to_string(e) + " edges";
}

#define COMMON_BOOST_GRAPH_PY_IFACE(gid) \
	class_< gid, boost::noncopyable >(#gid, init<>()) \
		.def("__init__", make_constructor(::detail::init_boost_graph< gid >)) \
		.def("__init__", make_constructor(::detail::init_boost_graph0< gid >)) \
		.def(init<size_t>()) \
		.def("num_vertices", &boost_num_vertices< gid >) \
		.def("num_edges",    &boost_num_edges< gid >) \
		.def("vertices",     &boost_vertices< gid >) \
		.def("edges",        &boost_edges< gid >) \
		.def("add_vertex",   &boost_add_vertex< gid >) \
		.def("add_edge",     &boost_add_edge< gid >) \
		.def("__repr__",     &make_string<gid, #gid> ) \
		; \
	typedef pair_<typename boost::graph_traits<gid>::vertex_iterator > gid ## _vip; \
	class_< gid ## _vip >( #gid "_vertex_range", no_init ) \
		 .def("__iter__" , py::range(& gid ## _vip::first, & gid ## _vip::second)) \
		; \
	typedef pair_<typename boost::graph_traits<gid>::edge_iterator > gid ## _eip; \
	class_< gid ## _eip >( #gid "_edge_range", no_init ) \
		 .def("__iter__" , py::range(& gid ## _eip::first, & gid ## _eip::second)) \
		 ;

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
	std::pair<vertex_descriptor, vertex_descriptor>
	// std::pair<unsigned, unsigned>
	operator()(E edg) const {
		auto s=boost::source(edg, _g);
		auto t=boost::target(edg, _g);
		return std::make_pair(s, t);
	}
	//idmap_type const& _m;
	G const& _g;
};


class my_vd{
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
}


BOOST_PYTHON_MODULE(_graph)
{
	namespace py = boost::python;

	//py::to_python_converter<my_vd, makemyvd>();
#if 0
	py::to_python_converter<const unsigned, makeunsigned>();
	py::to_python_converter<unsigned, makeunsigned>();

	py::to_python_converter<const long, makelong>();

	py::to_python_converter<unsigned long, makeunsignedlong>();
	py::to_python_converter<const unsigned long, makeunsignedlong>();
	py::to_python_converter<long, makelong>();
#endif

#if 1 // TMP. add_edge return type
	py::to_python_converter<std::pair<int, bool>, pair2tuple<int, bool> >(); // hack
#endif
	py::to_python_converter<std::pair<int, int>, pair2tuple<int, int> >();
	py::to_python_converter<std::pair<long, long>, pair2tuple<long, long> >();
	py::to_python_converter<std::pair<unsigned, unsigned>, pair2tuple<unsigned, unsigned> >();
	py::to_python_converter<std::pair<unsigned long, unsigned long>, pair2tuple<unsigned long, unsigned long> >();

	py::to_python_converter<std::pair<std::pair< unsigned long, unsigned long >, bool>, addedgeret<unsigned long> >();

// 	std::pair<boost::detail::undirected_edge_iter<std::_List_iterator<boost::list_edge<unsigned long, boost::no_property> >

//    to_python_converter<Vertex, cached_vertex_object<Graph> >();
//    to_python_converter<Edge, cached_edge_object<Edge> >();

	class_<std::pair<py::object, py::object> >("edge")
		    .def("source", &::detail::get_first)
			 .def("target", &::detail::get_second)
			 ;

//	class_<graph_bal_wrap>("graph_bal_wrap", init<>())
//		//.def(init< boost::python::list const&, boost::python::list const&>())
//		.def("__init__", make_constructor(::detail::graph_factory_n<graph_bal_wrap>))
//		.def("__init__", make_constructor(::detail::graph_factory_n0<graph_bal_wrap>))
//		COMMON_GRAPH_PY_IFACE(graph_bal_wrap)
//		;

//	class_<graph_balu_wrap>("graph_balu_wrap", init<>())
//		//.def(init< boost::python::list const&, boost::python::list const&>())
//		.def("__init__", make_constructor(::detail::graph_factory_n<graph_balu_wrap>))
//		.def("__init__", make_constructor(::detail::graph_factory_n0<graph_balu_wrap>))
//		COMMON_GRAPH_PY_IFACE(graph_balu_wrap)
//		;
//
//
#if 1
	typedef std::pair<wrap::graph<_balvvu>::edge_iterator,
	                  wrap::graph<_balvvu>::edge_iterator> graph_eip_1;
	class_< graph_eip_1 >( "graph_edge_range", no_init ) \
		 .def("__iter__" , range(& graph_eip_1::first, & graph_eip_1::second)) \
		;

	typedef std::pair<wrap::graph<_gsgvvu16, std::vector<uint16_t> >::edge_iterator,
	                  wrap::graph<_gsgvvu16, std::vector<uint16_t> >::edge_iterator> graph_eip_16;
	class_< graph_eip_16 >( "graph_edge_range", no_init ) \
		 .def("__iter__" , range(& graph_eip_16::first, & graph_eip_16::second)) \
		;

	typedef std::pair<wrap::graph<_gsgvvu32, std::vector<uint32_t>>::edge_iterator,
	                  wrap::graph<_gsgvvu32, std::vector<uint32_t>>::edge_iterator> graph_eip_32;
	class_< graph_eip_32 >( "graph_edge_range", no_init ) \
		 .def("__iter__" , range(& graph_eip_32::first, & graph_eip_32::second)) \
		;
#endif

	typedef std::pair<wrap::graph<_balvvu>::raw_vertex_iterator,
	                  wrap::graph<_balvvu>::raw_vertex_iterator> graph_vip_1;
	class_< graph_vip_1 >( "graph_vertex_range", no_init ) \
		 .def("__iter__" , range(& graph_vip_1::first, & graph_vip_1::second)) \
		;

#if 0 // ?

	typedef std::pair<graph ::raw_vertex_iterator, graph ::raw_vertex_iterator> graph_vip; \
	class_< graph_vip >( "graph_vertex_range", no_init ) \
		 .def("__iter__" , range(& graph_vip::first, & graph_vip::second)) \
		; \
	typedef std::pair<graph ::edge_iterator, graph ::edge_iterator> graph_eip; \
	class_< graph_eip >( "graph_edge_range", no_init ) \
		 .def("__iter__" , range(& graph_eip::first, & graph_eip::second)) \
		 ;
#endif

// class_<graph>("graph", init<>())
// 	//.def(init< boost::python::list const&, boost::python::list const&>())
// 	.def("__init__", make_constructor(::detail::graph_factory_n<graph>))
// 	.def("__init__", make_constructor(::detail::graph_factory_n0<graph>))
//	   COMMON_GRAPH_PY_IFACE(graph)
// 	;
// 
// 	class_<subgraph>("subgraph", init<>())
// 		//.def(init< boost::python::list const&, boost::python::list const&>())
// 		.def("__init__", make_constructor(::detail::graph_factory0<subgraph>))
// 		.def("__init__", make_constructor(::detail::graph_factory<subgraph>))
// 		COMMON_GRAPH_PY_IFACE(subgraph)
// 		;
// 
// 	class_<Graph>("Graph", init<>())
// 		//.def(init< boost::python::list const&, boost::python::list const&>())
// 		.def("__init__", make_constructor(::detail::graph_factory0<Graph>))
// 		.def("__init__", make_constructor(::detail::graph_factory<Graph>))
// 		COMMON_GRAPH_PY_IFACE(Graph)
// 		;
//

	typedef pair_<typename boost::graph_traits<_gsgvvu32>::edge_iterator > test_eip;
	class_< test_eip >( "test_edge_range", no_init )
		 .def("__iter__" , py::range(& test_eip::first, & test_eip::second)) ;

#ifdef HAVE_GALA_GRAPH_H
	COMMON_BOOST_GRAPH_PY_IFACE(_gsgvvu16)
	COMMON_BOOST_GRAPH_PY_IFACE(_gsgvvu32)
	COMMON_BOOST_GRAPH_PY_IFACE(_gsgvvu64)
#endif

	COMMON_BOOST_GRAPH_PY_IFACE(_balvvd)
	COMMON_BOOST_GRAPH_PY_IFACE(_balsvd)
	COMMON_BOOST_GRAPH_PY_IFACE(_balvvu)
	COMMON_BOOST_GRAPH_PY_IFACE(_balsvu)

	class_<exception_invalid> cpp_exception_invalid_class("exception_invalid", init<std::string>());
	cpp_exception_invalid_class.def("what", &exception_invalid::what);

	pythonExceptionType = createExceptionClass("exception_invalid");
	boost::python::register_exception_translator<exception_invalid>(&::detail::translate_exception_invalid);
}

//template class wrap::graph<graph_bal, boost::identity_property_map>;
//template class wrap::graph<graph_balu, boost::identity_property_map>;
//template class tf::graph<boost::identity_property_map>;
