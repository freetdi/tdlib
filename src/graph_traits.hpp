// Felix Salfelder, 2016
//
// (c) 2016 Goethe-Universität Frankfurt
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
//   traits for tdlib graphs and treedecs.
//
//

#ifndef TD_GRAPH_TRAITS_HPP
#define TD_GRAPH_TRAITS_HPP

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <set>

namespace treedec{
#ifndef TD_DEFS_NETWORK_FLOW
#define TD_DEFS_NETWORK_FLOW

struct Vertex_NF{
    bool visited;
    int predecessor;
};

struct Edge_NF{
    bool path; //true if a path uses the edge
};

#endif

struct bag_t{ //
    std::set<unsigned int> bag;
};

// dont define twice (in old code)
#define TD_STRUCT_BAG

}// treedec

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
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, bag_t> type;
};
// } not yet

// this makes some sense...
template<class G_t>
struct graph_traits : public graph_traits_base<G_t> { //
    typedef typename treedec_chooser<G_t>::type treedec_type;
    typedef typename std::set<unsigned> outedge_set_type;
    typedef typename boost::adjacency_list<boost::vecS, boost::vecS,
                boost::bidirectionalS, Vertex_NF, Edge_NF> directed_overlay;
    typedef typename boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> immutable_type;
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
inline bool is_valid(typename boost::graph_traits<G>::vertex_descriptor const &,
		const G&)
{
    return true;
}


namespace detail{ //
// working around balu and bag
template <bool, typename T = void>
struct vdstuff { //
    typedef unsigned type;
    typedef std::set<unsigned> bag_type;
};

// temporary hack don't touch.
// specialize teedec_traits (below) in case you must.
template <typename T>
struct vdstuff<false, T> { //
    typedef typename T::value_type type;
    typedef T bag_type;
};
} //detail

template<class T>
struct treedec_traits{ //
// TODO should be this (does not work, why?)
//    typedef typename boost::graph_traits<T>::vertex_property_type vertex_property_type;
    typedef typename T::vertex_property_type vertex_property_type;
    typedef typename detail::vdstuff<
       boost::is_same<vertex_property_type, bag_t >::value,
         vertex_property_type >::type vd_type;

    typedef typename detail::vdstuff<
       boost::is_same<vertex_property_type, bag_t >::value,
         vertex_property_type >::bag_type bag_type;

};

}

// return "id" where the vertex_descriptor might make more sense.
// (transitional interface)
template<typename G>
inline unsigned get_vd(const G&, const typename boost::graph_traits<G>::vertex_descriptor& v )
{
    // works with "TD_graph_t" (augmented adj_list)
    //return g[v].id;
    return v;
}

//Return the internal vertex position.
//To be used as a narrower alternative to vertex_descriptor.
//Positions are in {0, 1, ..., num_vertices-1}, where applicable.
//(One you use the vertex descriptor in boost graphs with vertex container 'vecS').
// this position must be stable under copy and assignment operations.
// // BUG: namespace
template<typename G_t>
inline unsigned
   get_pos(typename boost::graph_traits<G_t>::vertex_descriptor v, const G_t& G);

namespace treedec{

// chooose deg implementation for graph backend.
// to be accessed through graph_traits
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
struct vertex_callback{ //
    virtual ~vertex_callback(){};
    virtual void operator()(vertex_descriptor)=0;
};

template<typename G_t>
struct edge_callback{ //
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

} // treedec

#endif
