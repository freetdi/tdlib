// Felix Salfelder 2017
//
// (c) 2017 Felix Salfelder
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
//
#ifndef TDI_DIRECTED_VIEW_H
#define TDI_DIRECTED_VIEW_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
// #include <boost/graph/copy.hpp>
#include "treedec_traits.hpp"

namespace treedec{

namespace draft{

namespace detail{

template<class G, class X=void>
struct dwt{
	typedef typename graph_traits<G>::directed_type type;

	static std::string dbg(){ return "wrap directed\n"; }

	static size_t init(G const& g){
		return boost::num_vertices(g);
	}

	template<class GG, class H>
	static void copy(GG const& g, H& h){
		assert(boost::is_undirected(g));
		auto p=boost::edges(g);
		for(; p.first!=p.second; ++p.first){
			auto V=boost::source(*p.first, g);
			auto W=boost::target(*p.first, g);
			boost::add_edge(V, W, h);
			boost::add_edge(W, V, h);
		}
		trace2("", boost::num_edges(g), boost::num_edges(h));
		assert(boost::num_edges(g)*2 == boost::num_edges(h));
	}
};

template<class G>
struct dwt<G,
	typename std::enable_if< std::is_convertible<
	typename boost::graph_traits<G>::directed_category,
  	boost::directed_tag>::value && !std::is_convertible<
   typename boost::graph_traits<G>::directed_category,
   boost::bidirectional_tag>::value, void >::type >
{
	typedef typename graph_traits<G>::directed_type& type;
	static std::string dbg(){ return "dummy wrapper\n"; }

	static G& init(G& g){
		return g;
	}

	template<class GG, class H>
	static void copy(GG const&, H&){
	}
};

template<class G>
struct dwt<G,
	typename std::enable_if< std::is_convertible<
	typename boost::graph_traits<G>::directed_category,
  	boost::bidirectional_tag>::value, void >::type >
{
	typedef typename graph_traits<G>::directed_type type;
	static std::string dbg(){ return "bidir wrapper\n"; }

	static size_t init(G&){
		return 0;
	}

	// check: do we need this copy?
	template<class GG, class H>
	static void copy(GG const& g, H& h){
		assert(!boost::num_vertices(h));
		assert(boost::is_directed(h));
		auto p=boost::edges(g);
		for(; p.first!=p.second; ++p.first){
			auto e=*p.first;
			auto V=boost::source(e, g);
			auto W=boost::target(e, g);
			boost::add_edge(V, W, h);
			boost::add_edge(W, V, h);
			trace2("bidir cp", W, V);
		}
		trace2("bidir cp", boost::num_edges(g), boost::num_edges(h));
		assert(2*boost::num_edges(g) == boost::num_edges(h));
	}
};

}

template<class G>
class directed_view{

public:
	typedef detail::dwt<G> wrapper_help;
	typedef typename detail::dwt<G>::type backend_type;
	typedef typename std::remove_reference<backend_type>::type wrapped_type;

	typedef typename boost::graph_traits<wrapped_type> wrapped_traits;
	typedef typename wrapped_traits::vertex_descriptor vertex_descriptor;
	typedef typename wrapped_traits::edge_descriptor edge_descriptor;
	typedef typename wrapped_traits::out_edge_iterator out_edge_iterator;
	typedef typename wrapped_traits::adjacency_iterator adjacency_iterator;
	typedef typename wrapped_traits::vertex_iterator vertex_iterator;
	typedef typename wrapped_traits::edge_iterator edge_iterator;
	typedef typename wrapped_traits::vertices_size_type vertices_size_type;
	typedef typename wrapped_traits::degree_size_type degree_size_type;
	typedef typename wrapped_traits::edges_size_type edges_size_type;

	typedef typename wrapped_traits::edge_parallel_category edge_parallel_category;

	typedef typename boost::directed_tag directed_category;

//	typedef typename wrapped_traits::traversal_category traversal_category;
	typedef typename boost::adjacency_graph_tag traversal_category;

private:
	directed_view(){ unreachable(); }
	directed_view(const directed_view& ) { unreachable();}
public:
	directed_view(G& g, bool commit=false)
	 : _g(wrapper_help::init(g)),
	   _commit(commit)
	{
		assert(boost::is_directed(_g));

		// no, only copies one edge per edge
		//boost::copy_graph(g, _g);
		//
		trace1("wrapping", wrapper_help::dbg());

		wrapper_help::copy(g, _g);
	}

	~directed_view(){
		if(_commit){ untested();
			incomplete();
		}else{
		}
	}

public: // dangerous?
	wrapped_type* operator->(){ return &_g; }
	wrapped_type& operator*(){ return _g; }
	wrapped_type const& operator*() const{ return _g; }
public: // unnecessary
	std::pair< typename boost::graph_traits<G>::vertex_iterator,
		typename boost::graph_traits<G>::vertex_iterator>
	vertices() const{
		return boost::vertices(_g);
	}
	std::pair<adjacency_iterator, adjacency_iterator>
	adjacent_vertices(vertex_descriptor v) const{
		return boost::adjacent_vertices(v, _g);
	}
	typename boost::graph_traits<G>::vertices_size_type
	degree(typename boost::graph_traits<G>::vertex_descriptor v) const{ untested();
      wrapped_type* g = const_cast<wrapped_type*>(&_g);
		return boost::out_degree(v, *g);
	}
private:
	backend_type _g;
	bool _commit;
}; // directed_view

} // draft

#if 1
// default directed view
// there are better ways
template<class G>
struct directed_view_select{
    typedef treedec::draft::directed_view<G> type;
};
#endif

} // treedec

namespace boost{

template<class G>
struct graph_traits<treedec::draft::directed_view<G> > {
	typedef treedec::draft::directed_view<G> B;
	typedef typename treedec::draft::directed_view<G>::wrapped_type WW;
	typedef typename std::remove_reference<WW>::type W;

	typedef typename B::vertices_size_type vertices_size_type;
	typedef typename B::degree_size_type degree_size_type;
	typedef typename B::edges_size_type edges_size_type;

	typedef typename B::vertex_descriptor vertex_descriptor;
	typedef typename B::edge_descriptor edge_descriptor;

	typedef typename B::vertex_iterator vertex_iterator;
	typedef typename B::edge_iterator edge_iterator;
	typedef typename B::out_edge_iterator out_edge_iterator; //?
	typedef typename B::directed_category directed_category;
	typedef typename B::adjacency_iterator adjacency_iterator;

	// hmm.
	// typedef typename boost::disallow_parallel_edge_tag edge_parallel_category;
	typedef typename B::edge_parallel_category edge_parallel_category;

	typedef typename B::traversal_category traversal_category;

	static vertex_descriptor null_vertex(){
		return 0; // graph_traits<W>::null_vertex();
	}
};

template<class G, class T>
struct property_map<treedec::draft::directed_view<G>, T> {
	typedef typename property_map<
		typename treedec::draft::directed_view<G>::wrapped_type, T>::const_type const_type;
	typedef typename property_map<
		typename treedec::draft::directed_view<G>::wrapped_type, T>::type type;
};

template<class G>
degree_property_map<treedec::draft::directed_view<G> >
get(vertex_degree_t t, treedec::draft::directed_view<G> const& g)
{
	return get(t, *g);
}

template<class G>
degree_property_map<treedec::draft::directed_view<G> >
get(vertex_degree_t t, treedec::draft::directed_view<G>& g)
{
	return get(t, *g);
}

template<class T, class G>
typename property_map<typename treedec::draft::directed_view<G>::wrapped_type,
                      T>::type
get(T t, treedec::draft::directed_view<G>& g)
{
	return get(t, *g);
}

template<class T, class G>
typename property_map<typename treedec::draft::directed_view<G>::wrapped_type,
                      T>::const_type
get(T t, treedec::draft::directed_view<G> const& g)
{
	return get(t, *g);
}

template<class T, class G, class V>
typename property_map<typename treedec::draft::directed_view<G>::wrapped_type,
                      T>::value_type
get(T t, treedec::draft::directed_view<G> const& g, V v)
{
	return get(t, *g, v);
}

template<class T, class G, class V, class W>
void put(T t, treedec::draft::directed_view<G>& g, V v, W w)
{
	return put(t, *g, v, w);
}

template<class G>
typename property_map<typename treedec::draft::directed_view<G>::wrapped_type,
                      boost::edge_all_t>::const_type
get(edge_all_t, treedec::draft::directed_view<G> const& g)
{
	return get(edge_all, *g);
}

template<class G>
typename treedec::draft::directed_view<G>::vertices_size_type
num_vertices(treedec::draft::directed_view<G> const& g)
{
	return num_vertices(*g);
}

template<class G>
typename graph_traits<G>::edges_size_type
num_edges(treedec::draft::directed_view<G> const& g)
{
	return num_edges(*g);
}

template<class VD, class G>
typename graph_traits<G>::vertices_size_type
out_degree(VD v, treedec::draft::directed_view<G> const& g)
{
	return out_degree(v, *g);
}

template<class VD, class G>
typename graph_traits<G>::vertices_size_type
	// BUG: out_degree!
degree(VD v, treedec::draft::directed_view<G> const& g)
{
	return out_degree(v, *g);
}

template<class VD, class G>
void clear_vertex(VD v, treedec::draft::directed_view<G>& g)
{
	return clear_vertex(v, *g);
}

template<class VD, class G>
std::pair<typename treedec::draft::directed_view<G>::edge_descriptor, bool >
add_edge(VD v, VD w, treedec::draft::directed_view<G>& g)
{
#ifndef NDEBUG
	// hmm directed_view disallows multiedges.
	if(edge(v, w, *g).second){
		incomplete();
		// that will not work with vectors
		return std::make_pair(edge(v, w, *g).first, true);
	}
#endif
	return add_edge(v, w, *g);
}

template<class G>
void remove_edge(
		typename treedec::draft::directed_view<G>::vertex_descriptor v,
		typename treedec::draft::directed_view<G>::vertex_descriptor w,
		treedec::draft::directed_view<G>& g)
{
	return remove_edge(v, w, *g);
}

template<class G, class P>
void remove_out_edge_if(
		typename treedec::draft::directed_view<G>::vertex_descriptor v,
		P const& p,
		treedec::draft::directed_view<G>& g)
{
	return remove_out_edge_if(v, p, *g);
}

template<class G>
std::pair<typename treedec::draft::directed_view<G>::edge_descriptor, bool >
edge(
	typename treedec::draft::directed_view<G>::vertex_descriptor v,
	typename treedec::draft::directed_view<G>::vertex_descriptor w,
	treedec::draft::directed_view<G>& g)
{
	return edge(v, w, *g);
}

template<class G>
typename treedec::draft::directed_view<G>::vertex_descriptor
source(typename treedec::draft::directed_view<G>::edge_descriptor w,
		treedec::draft::directed_view<G> const& g)
{
	return source(w, *g);
}

template<class G>
typename treedec::draft::directed_view<G>::vertex_descriptor
target(typename treedec::draft::directed_view<G>::edge_descriptor w,
		treedec::draft::directed_view<G> const& g)
{
	return target(w, *g);
}

template<class G>
std::pair<
	typename treedec::draft::directed_view<G>::vertex_iterator,
	typename treedec::draft::directed_view<G>::vertex_iterator>
vertices(treedec::draft::directed_view<G> const& g)
{
	return vertices(*g);
}

template<class G>
std::pair<
	typename treedec::draft::directed_view<G>::edge_iterator,
	typename treedec::draft::directed_view<G>::edge_iterator>
edges(treedec::draft::directed_view<G> const& g)
{
	return edges(*g);
}

template<class G>
std::pair<
	typename treedec::draft::directed_view<G>::out_edge_iterator,
	typename treedec::draft::directed_view<G>::out_edge_iterator>
out_edges(
		typename treedec::draft::directed_view<G>::vertex_descriptor v,
		treedec::draft::directed_view<G> const& g)
{
	return out_edges(v, *g);
}

template<class G>
std::pair<
	typename treedec::draft::directed_view<G>::adjacency_iterator,
	typename treedec::draft::directed_view<G>::adjacency_iterator>
adjacent_vertices(
		typename treedec::draft::directed_view<G>::vertex_descriptor v,
		treedec::draft::directed_view<G> const& g)
{
	return g.adjacent_vertices(v);
}

    template<class G>
    struct vertex_bundle_type<treedec::draft::directed_view<G> > {
      typedef typename vertex_bundle_type<G>::type type;
    };

} // boost

#endif
