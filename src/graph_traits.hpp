// Felix Salfelder, 2016
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
//   traits for tdlib graphs.

#ifndef TREEDEC_GRAPH_TRAITS_HPP
#define TREEDEC_GRAPH_TRAITS_HPP

// BUG, some old code depends on this (and shouldn't)
#define TD_GRAPH_TRAITS_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <set>
#include "trace.hpp"
#include "container.hpp" // HACK

namespace treedec{

// BUG
#if defined(TREEDEC_DEFS_NETWORK_FLOW) || defined(TD_DEFS_NETWORK_FLOW)
#warning deprecated macro use
#else
#define TREEDEC_DEFS_NETWORK_FLOW

// old BUG (already spread out)
#define TD_DEFS_NETWORK_FLOW

struct bagsize_t{
	unsigned dummy;
};

} // treedec

namespace boost{

// this is forbidden?
template<class G>
unsigned& get(treedec::bagsize_t, G&){
	unreachable();
	static unsigned udummy;
	return udummy;
}

} // boost

namespace treedec{


struct Vertex_NF{
    bool visited;
    int predecessor;
};

struct Edge_NF{
    bool path; //true if a path uses the edge
};

#endif

// ouch. this is actually bag_t
// how to fix that now?
// "bag" is also used as shorthand for bag access...

// BUG: move to treedec_traits

struct bag_t{ //
  typedef boost::vertex_property_tag kind;
  std::set<unsigned int> bag; // yikes. old way.
  bag_t& operator=(std::set<unsigned int> const& b){
	  bag = b;
	  return *this;
  }
  bag_t& operator=(std::vector<unsigned int> const& b){ untested();
	  bag.clear();
	 	for( auto& i : b){ untested();
			bag.insert(i);
		}
	  return *this;
  }
  bag_t& operator=(boost::property<bag_t, std::vector<unsigned int> > const& b)
  {
	  bag.clear();
	 	for( auto& i : b.m_value){ untested();
			bag.insert(i);
		}
	  return *this;
  }
  bag_t& operator=(boost::property<treedec::bag_t, std::set<unsigned int> > const& b)
  { untested();
	  bag = b.m_value;
	  return *this;
  }

  operator std::set<unsigned>() const{
	  return bag;
  }
};

// dont define twice (in old code)
#define TREEDEC_STRUCT_BAG

}// treedec

// KLUGE: put it here...
// (and cross fingers)
//using bag = treedec::bag_t;

namespace treedec{


template<class G_t>
struct graph_traits_base : public boost::graph_traits<G_t> {};
// kludge for balu
// TODO: use graph_traits. see below
// namespace detail{ not yet
template<class G>
struct treedec_chooser{ //
    typedef unsigned value_type;
    typedef std::set<unsigned> bag_type;
#ifdef AVOID_BUNDLES_PROPERTY_BAGS
	 typedef boost::adjacency_list<boost::vecS, boost::vecS,
	                           boost::undirectedS,
										boost::property<bag_t, bag_type> > type;
#else
	 typedef boost::adjacency_list<boost::vecS, boost::vecS,
	                           boost::undirectedS, bag_t> type;
#endif
};
// } not yet
//
namespace detail{

template<class G_t>
struct default_directed_select{
   typedef void type; // bug?
};

// maybe just all? be conservative
template<class X, class Y, class Z>
struct default_directed_select< boost::adjacency_list<X, Y, boost::directedS, Z> >{
   typedef boost::adjacency_list<X, Y, boost::directedS, Z> type;
};

template<class X, class Y, class Z>
struct default_directed_select< boost::adjacency_list<X, Y, boost::bidirectionalS, Z> >{
   typedef boost::adjacency_list<X, Y, boost::directedS, Z> type;
};

template<class X, class Y, class Z>
struct default_directed_select< boost::adjacency_list<X, Y, boost::undirectedS, Z> >{
   typedef boost::adjacency_list<X, Y, boost::directedS, Z> type;
};

} // detail

// this makes some sense...
template<class G_t>
struct graph_traits : public graph_traits_base<G_t> { //
    typedef typename treedec_chooser<G_t>::type treedec_type;
    typedef typename std::set<unsigned> outedge_set_type;
    typedef typename boost::adjacency_list<boost::vecS, boost::vecS,
                boost::bidirectionalS, Vertex_NF, Edge_NF> directed_overlay;
    typedef typename boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> immutable_type;
    typedef typename detail::default_directed_select<G_t>::type directed_type;
};

// obsolete. use graph_traits directly.
template<class G>
struct outedge_set{
    typedef typename graph_traits<G>::outedge_set_type type;
};

// test if v is a valid vertex_descriptor of g
template<typename G>
inline bool is_valid(typename boost::graph_traits<G>::vertex_iterator const &,
		const G&)
{
    return true;
}
template<typename G>
inline bool is_valid(typename boost::graph_traits<G>::vertex_descriptor const & ,
		const G& )
{
	// bug: missing gala override
   // return v < boost::num_vertices(g);
	return true;
}


namespace detail{
// working around balu and bag
template <bool, typename T = void>
struct vdstuff {
    typedef unsigned type;
    typedef std::set<unsigned> bag_type;
};

// temporary hack don't touch.
// specialize teedec_traits (below) in case you must.
template <typename T>
struct vdstuff<false, T> {
	// T is not bag_t. so probably a boost::property
    typedef typename T::value_type::value_type type;
    typedef typename T::value_type bag_type;
};
} //detail

template<class T>
struct treedec_traits{
// TODO should be this (does not work, why?)
//    typedef typename boost::graph_traits<T>::vertex_property_type vertex_property_type;
    typedef typename T::vertex_property_type vertex_property_type;
    typedef typename detail::vdstuff<
       boost::is_same<vertex_property_type, bag_t >::value,
         vertex_property_type >::type vd_type;

    typedef typename detail::vdstuff<
       boost::is_same<vertex_property_type, treedec::bag_t >::value,
          vertex_property_type >::bag_type bag_type;

};

} // treedec

// return "id" where the vertex_descriptor might make more sense.
// (transitional interface)
template<typename G>
inline unsigned get_vd(const G&, const typename boost::graph_traits<G>::vertex_descriptor& v )
{
    // works with "TREEDEC_graph_t" (augmented adj_list)
    //return g[v].id;
    return v;
}

//Return the internal vertex position.
//To be used as a narrower alternative to vertex_descriptor.
//Positions are in {0, 1, ..., num_vertices-1}, where applicable.
//(One you use the vertex descriptor in boost graphs with vertex container 'vecS').
// this position must be stable under copy and assignment operations.

namespace treedec{

// chooose deg implementation for graph backend.
// to be accessed through graph_traits
// obsolete. now in CFG, if necessary at all.
template<class G_t>
struct deg_chooser;

template<class G>
void check(G const&)
{
}

namespace detail{
  template<class G>
  class shared_adj_iter;
}

template<typename vertex_descriptor>
struct vertex_callback{
    virtual ~vertex_callback(){};
    virtual void operator()(vertex_descriptor)=0;
};

template<typename G_t>
struct edge_callback{
    typedef typename boost::graph_traits<G_t>::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    virtual ~edge_callback(){};
    virtual void operator()(vertex_descriptor, vertex_descriptor)=0;
    void operator()(edge_descriptor)
    { incomplete();
    }
};

template<typename G_t>
struct graph_callback{ // fixme: union of the above?
    typedef typename boost::graph_traits<G_t>::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    virtual ~graph_callback(){};
    virtual void operator()(vertex_descriptor)=0;
    virtual void operator()(vertex_descriptor, vertex_descriptor)=0;
};

template<class G>
inline std::pair<typename boost::graph_traits<G>::edge_descriptor, bool>
add_edge(typename boost::graph_traits<G>::vertex_descriptor x,
		   typename boost::graph_traits<G>::vertex_descriptor y, G& g);

template<class G>
inline typename boost::graph_traits<G>::edges_size_type num_edges(G const& g);

template<class G>
struct graph_helper{
    // stub. incomplete.
	static_assert(sizeof(G)==0, "need specialization");

	template<class S>
	static void close_neighbourhood(S&, G const&){
		static_assert(sizeof(S)==0, "need specialization");
	}
	template<class S>
	static void open_neighbourhood(S&, G const&){
		static_assert(sizeof(S)==0, "need specialization");
	}
};

} // treedec


//  template <class PropertyMap, class Reference, class K, class V>
//  inline void
//  put(const put_get_helper<Reference, PropertyMap>& pa, K k, const V& v)
//  {
//    static_cast<const PropertyMap&>(pa)[k] = v;
//  }
//

namespace boost {

// TODO: move to treedec_traits?
namespace bagstuff {

	template<class T, class X=void>
	struct gtob{
		typedef decltype( boost::vertex_bundle_type<T>::type::bag ) type;
	};
	template<class T>
	struct gtob<T, decltype( boost::vertex_bundle_type<T>::type::bag )>{
		typedef decltype( boost::vertex_bundle_type<T>::type::bag ) type;
	};

	template<class G>
	struct const_treebagpmap
		: public put_get_helper<typename gtob<G>::type, const_treebagpmap<G> >
	{ //
		typedef typename gtob<G>::type B;
		const_treebagpmap(G const& g) : _g(g){}

		std::set<unsigned>& operator[](unsigned v) const{
			auto& g=const_cast<G&>(_g); // huh?
			return g[v].bag;
		}
		G const&_g;
	};
	template<class G>
	struct treebagpmap : public put_get_helper<typename gtob<G>::type, treebagpmap<G> >
	{
		typedef typename gtob<G>::type B;
		treebagpmap(G& g) : _g(g){
			untested();
		}
		B & operator[](unsigned v) const{
			auto& g=const_cast<G&>(_g); // huh?
			return g[v].bag;
		}

		B& operator[](unsigned v){
			return _g[v].bag;
		}
		G&_g;
	};

} // bagstuff


	template<class C, class T>
	struct bag_sfinae { typedef T type; };

	template<class G>
	struct property_map<G,
		typename bag_sfinae<
			decltype( boost::vertex_bundle_type<G>::type::bag ),
		   typename std::enable_if<
			   1||std::is_same< typename boost::vertex_bundle_type<G>::type,
			                  treedec::bag_t>::value, vertex_all_t
							>::type
		 >::type
		 >
	{
		typedef bagstuff::treebagpmap<G>  type;
		typedef bagstuff::const_treebagpmap<G> const_type;
	};


} // boost

#endif
