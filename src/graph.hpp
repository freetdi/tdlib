// Felix Salfelder, 2015 - 2016
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

/*
 * Graph routines related to tree decompositions.
 *
 * Implemented for generic BGL graphs. Possibly (much) faster for particular
 * graphs.
 *
 * template<typename nIter_t, typename G_t>
 *    void make_clique(nIter_t nIter, G_t &G)
 *
 * template<typename B_t, typename nIter_t, typename G_t>
 *    void copy_neighbourhood(B_t &B, nIter_t nIter, G_t &G)
 *
 * template<class G>
 * void detach_neighborhood(typename boost::graph_traits<G>::vertex_descriptor& c,
 *    G& g, typename noboost::outedge_set<G>::type& bag)
 */

#ifndef TD_NOBOOST_H
#define TD_NOBOOST_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "degree.hpp"
#include "std.hpp"
#include "trace.hpp"

// OUCH
//#include "fill.hpp"

#ifndef TD_STRUCT_BAG
#define TD_STRUCT_BAG
struct bag{
    std::set<unsigned int> bag;
    typedef std::set<unsigned int> bagtype;
};
#endif

namespace noboost{
#if 0 // later, need c++11
    template<typename G>
    using vertex_iterator = typename boost::graph_traits<G>::vertex_iterator;
    template<typename G>
    using vertex_descriptor = typename boost::graph_traits<G>::vertex_descriptor;
    template<typename G>
    using adjacency_iterator = typename boost::graph_traits<G>::adjacency_iterator;
#define vertex_iterator_G vertex_iterator<G>
#define vertex_descriptor_G typename vertex_descriptor<G>
#define adjacency_iterator_G typename adjacency_iterator<G>
#else
#define vertex_iterator_G typename boost::graph_traits<G>::vertex_iterator
#define vertex_descriptor_G typename boost::graph_traits<G>::vertex_descriptor
#define adjacency_iterator_G typename boost::graph_traits<G>::adjacency_iterator
#endif

template<class G>
void check(G const&)
{
}

template<typename G>
void remove_vertex(vertex_iterator_G u, G &g)
{
    remove_vertex(*u, g);
}

template<typename vertex_descriptor>
struct vertex_callback{
    virtual ~vertex_callback(){};
    virtual void operator()(vertex_descriptor)=0;
};

template<typename G_t>
struct edge_callback{
    typedef typename boost::graph_traits<G_t>::edge_descriptor edge_descriptor;
    virtual ~edge_callback(){};
    virtual void operator()(edge_descriptor)=0;
};

}

namespace treedec{
template<typename G_t>
struct graph_callback{ // fixme: union of the above?
    typedef typename boost::graph_traits<G_t>::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    virtual ~graph_callback(){};
    virtual void operator()(edge_descriptor)=0;
    virtual void operator()(vertex_descriptor)=0;
};
} //treedec

namespace noboost{

//Vertex v will remain as isolated node.
//Calls cb on neighbors if degree drops by one,
//before dg will drop.
// FIXME: pass edge or edge iterator.
template<typename G>
void contract_edge(vertex_descriptor_G v,
                   vertex_descriptor_G target,
                   G &g,
                   bool /*erase*/=false,
                   vertex_callback<vertex_descriptor_G>* cb=NULL)
{ untested();
    adjacency_iterator_G I, E;
    for(boost::tie(I, E)=boost::adjacent_vertices(v, g); I!=E; ++I){
        assert(boost::edge(v, *I, g).second);
        if(*I != target){
            bool added=boost::add_edge(target, *I, g).second;
            if(added){
                //rebasing edge from I-v to I-target.
            }else if(cb){
                //did not add, degree will drop by one.
                (*cb)(*I);
            }
        }else{
        }
    }

    boost::clear_vertex(v, g);
    assert(!boost::edge(v, target, g).second);
}

//Vertex v will remain as isolated node, unless erase.
//If erase, v will be deleted (not collapsed).
template<typename G>
inline void contract_edge(vertex_iterator_G v,
                   vertex_descriptor_G into,
                   G &g,
                   bool erase=false,
                   vertex_callback<vertex_descriptor_G>* cb=NULL)
{ untested();
    contract_edge(*v, into, g, erase, cb);
    if(erase){
        noboost::remove_vertex(v, g);
    }
}

// turn vertex range into clique.
// call cb on newly created edges and incident vertices.
template<typename B, typename E, typename G_t>
size_t /*hmm*/ make_clique(B nIt1, E nEnd, G_t &G, treedec::graph_callback<G_t>* cb=NULL)
{
    typedef typename treedec::graph_callback<G_t> CB;
    size_t counter=0;
    B nIt2;
    for(; nIt1!=nEnd; ++nIt1){
#if 1
        if(cb){
            (*cb)(*nIt1); // hmm
        }
#endif
        nIt2 = nIt1;
        nIt2++;
        for(; nIt2 != nEnd; nIt2++){

            BOOST_AUTO(ep, boost::add_edge(*nIt1, *nIt2, G));
            if(ep.second){ untested();
               ++counter;

               if(cb){
#if 0 // doesn't work. removing center...?
            (*cb)(*nIt1); // hmm
            (*cb)(*nIt2); // hmm
#endif
                   (*cb)(ep.first);
               }
            }
        }
    }
    return counter;
}

// convenience wrapper (boost-style iterator pair)
template<typename nIter_t, typename G_t>
void make_clique(nIter_t nIter, G_t &G, treedec::graph_callback<G_t>* cb=NULL)
{ itested();
    typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nEnd;
    boost::tie(nIt1, nEnd) = nIter;
    make_clique(nIt1, nEnd, G, cb);
}

// FIXME: wrong name, unused?
template<typename B_t, typename nIter_t, typename G_t>
// void copy_vertex_range
void fetch_neighbourhood(B_t &B, nIter_t nIter, G_t &G)
{ untested();
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = nIter; nIt != nEnd; nIt++){
        B.insert(*nIt);
    }
}

template <typename G_t>
inline typename boost::graph_traits<G_t>::vertex_descriptor
   get_min_degree_vertex(const G_t &G, bool ignore_isolated_vertices=false)
{
    unsigned int min_degree = UINT_MAX;

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    boost::tie(vIt, vEnd) = boost::vertices(G);
    typename boost::graph_traits<G_t>::vertex_descriptor min_vertex = *vIt++;
    for(; vIt != vEnd; vIt++){
        unsigned int degree = boost::degree(*vIt, G);
        if(degree <= min_degree){
            if(ignore_isolated_vertices && degree == 0){ continue; }
            min_degree = degree;
            min_vertex = *vIt;
        }
    }

    return min_vertex;
}

template <typename G_t>
inline typename boost::graph_traits<G_t>::vertex_descriptor
   get_least_common_vertex(const typename boost::graph_traits<G_t>::vertex_descriptor &min_vertex, const G_t &G)
{
    typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
    boost::tie(nIt1, nEnd) = boost::adjacent_vertices(min_vertex, G);
    typename boost::graph_traits<G_t>::vertex_descriptor w = *nIt1;

    unsigned int min_common = UINT_MAX;

    for(; nIt1 != nEnd; nIt1++){
        unsigned int cnt_common = 0;
        nIt2 = boost::adjacent_vertices(min_vertex, G).first;
        for(; nIt2 != nEnd; nIt2++){
            if(boost::edge(*nIt1, *nIt2, G).second){
                cnt_common++;
            }
        }
        if(cnt_common < min_common){
            w = *nIt1;
            min_common = cnt_common;
        }
    }

    return w;
}

template<typename G_t>
unsigned int eliminate_vertex(typename boost::graph_traits<G_t>::vertex_descriptor v, G_t &G){
    noboost::make_clique(boost::adjacent_vertices(v, G), G);
    unsigned int deg = boost::degree(v, G);
    boost::clear_vertex(v, G);
    return deg;
}

template <typename G_t>
inline void make_degree_sequence(const G_t &G,
          std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &degree_sequence)
{
    unsigned int max_degree = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::degree(*vIt, G);
        max_degree = (degree>max_degree)? degree : max_degree;
    }

    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > buckets(max_degree+1);
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::degree(*vIt, G);
        if(degree > 0){
            buckets[degree].push_back(*vIt);
        }
    }
    for(unsigned int i = 1; i <= max_degree; i++){
        for(unsigned int j = 0; j < buckets[i].size(); j++){
            degree_sequence.push_back(buckets[i][j]);
        }
    }
}

//Return the internal vertex position.
//To be used as a narrower alternative to vertex_descriptor.
//Positions are in {0, 1, ..., num_vertices-1}, where applicable.
//(One you use the vertex descriptor in boost graphs with vertex container 'vecS').
// this position must be stable under copy and assignment operations.
template<typename G_t>
inline unsigned int get_pos(const typename boost::graph_traits<G_t>::vertex_descriptor v, const G_t& G){
    return boost::get(boost::get(boost::vertex_index, G), v);
}

// return "id" where the vertex_descriptor might make more sense.
// (transitional interface)
template<typename G>
inline unsigned get_vd(const G& g, const vertex_descriptor_G& v )
{
    // works with "TD_graph_t" (augmented adj_list)
    //return g[v].id;
    return v;
}

// test if v is a valid vertex_descriptor of g
template<typename G>
inline bool is_valid(const vertex_descriptor_G& v, const G& g)
{
    // base case: don't know. lets say yes (better in assertions)
    return true;
}

template<class G>
struct outedge_set{
    typedef std::set<unsigned> type;
//	typedef std::set type;
};

// kludge for balu
template<class G>
struct treedec_chooser{
    typedef unsigned value_type;
    typedef std::set<unsigned> bag_type;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, bag> type;
};

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
    typedef typename T::value_type type;
    typedef T bag_type;
};
} //detail

// FIXME: not part of noboost
template<class T>
struct treedec_traits{
    typedef typename detail::vdstuff<
       boost::is_same< typename T::vertex_property_type, bag >::value,
        typename T::vertex_property_type >::type vd_type;

    typedef typename detail::vdstuff<
       boost::is_same< typename T::vertex_property_type, bag >::value,
        typename T::vertex_property_type >::bag_type bag_type;

};

namespace detail{
template<class B, class T, class V>
struct tmpbaghack{
    static typename treedec_traits<T>::bag_type& get_bag(T& t, V& v)
    {
        return t[v];
    }
    static typename treedec_traits<T>::bag_type const& get_bag(T const& t, V const& v)
    {
        return t[v];
    }
};

template<class T_t, class V>
struct tmpbaghack<bag, T_t, V>{
    static typename treedec_traits<T_t>::bag_type& get_bag(T_t& t, V& v)
    {
        return t[v].bag;
    }
    static typename treedec_traits<T_t>::bag_type const& get_bag(T_t const& t, V const& v)
    {
        return t[v].bag;
    }
};
} //detail

template<typename T_t>
inline typename treedec_traits<T_t>::bag_type& bag(
	const typename boost::graph_traits<T_t>::vertex_descriptor& v,
        T_t& T)
{
    typedef    typename T_t::vertex_property_type b; //>::bag_type b;
    return detail::tmpbaghack<b,T_t,const typename boost::graph_traits<T_t>::vertex_descriptor&>::get_bag(T, v);
}

template<typename T_t>
inline typename treedec_traits<T_t>::bag_type const& bag(
        const typename boost::graph_traits<T_t>::vertex_descriptor& v,
        T_t const& T)
{
    typedef    typename T_t::vertex_property_type b; //>::bag_type b;
    return detail::tmpbaghack<b,T_t,const typename boost::graph_traits<T_t>::vertex_descriptor&>::get_bag(T, v);
}

template<class G_t>
struct deg_chooser{
    typedef typename misc::DEGS<G_t> type;
    typedef type degs_type; // transition? don't use.
};

} // namespace noboost


namespace treedec{


// put neighbors of c into bag
// isolate c in g
// return bag
// optionally: pass pointer to bag for storage.
template<class G>
inline void detach_neighborhood(
        typename boost::graph_traits<G>::vertex_descriptor& c,
        G& g, typename noboost::outedge_set<G>::type& bag);


// count number of edges missing in 1-neighborhood of v
template <typename G_t>
inline size_t count_missing_edges(
        const typename boost::graph_traits<G_t>::vertex_descriptor v, G_t const &G)
{ itested();
    size_t missing_edges = 0;

    typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
    for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(v, G); nIt1 != nEnd; nIt1++){
        nIt2 = nIt1;
        nIt2++;
        for(; nIt2 != nEnd; nIt2++){
            if(!boost::edge(*nIt1, *nIt2, G).second){
                ++missing_edges;
            }
        }
    }
    return missing_edges;
}
} // treedec

namespace treedec{

    // TODO: cleanup
using noboost::outedge_set;

// collect neighbors of c into bag
// remove c from graph (i.e. isolate).
// turn subgraph induced by bag into clique.
//
// optionally:
// - call back on each newly created edge, and on each vertex incident to such
//   edge.
// - provide memory through bag pointer.
//
template<class G>
size_t /*G::edge_index_type?*/ make_clique_and_detach(
        typename boost::graph_traits<G>::vertex_descriptor c,
        G& g,
        typename outedge_set<G>::type& bag,
        treedec::graph_callback<G>* cb=NULL)
{
    detach_neighborhood(c, g, bag);
    return noboost::make_clique(bag.begin(), bag.end(), g, cb);
}
}//treedec


namespace treedec{
namespace detail{

    // iterate over edges adjacent to both v and s
    // implementation: iterate outedges(v), skip non-outedges(s).
    template<class G>
    class shared_adj_iter : public boost::graph_traits<G>::adjacency_iterator {
    public:
        typedef typename boost::graph_traits<G>::adjacency_iterator adjacency_iterator;
        typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
        shared_adj_iter(adjacency_iterator v, adjacency_iterator ve,
                        vertex_descriptor s, G const& g)
            : adjacency_iterator(v), _ve(ve), _s(s), _g(g)
        {
            skip();
        }
        shared_adj_iter(const shared_adj_iter& p)
            : adjacency_iterator(p), _ve(p._ve), _s(p._s), _g(p._g)
        {
        }

        shared_adj_iter& operator++(){ untested();
            assert(_ve!=adjacency_iterator(*this));
            assert(adjacency_iterator(*this)!=_ve);
            adjacency_iterator::operator++();
            skip();

            return *this;
        }
    private:
        void skip()
        {
            while(true){
                if(typename boost::graph_traits<G>::adjacency_iterator(*this)==_ve){ untested();
                    return;
                }else if(!boost::edge(**this, _s, _g).second){ untested();
                    adjacency_iterator::operator++();
                }else{ untested();
                    return;
                }
            }
        }
        adjacency_iterator _ve;
        vertex_descriptor _s;
        G const& _g;
    };

}// detail

template<class G>
std::pair<detail::shared_adj_iter<G>, detail::shared_adj_iter<G> >
    common_out_edges(typename boost::graph_traits<G>::vertex_descriptor v,
                     typename boost::graph_traits<G>::vertex_descriptor w,
                     const G& g)
{ itested();
    BOOST_AUTO(p, boost::adjacent_vertices(v, g));
    typedef typename detail::shared_adj_iter<G> Iter;
    return std::make_pair(Iter(p.first, p.second, w, g),
                          Iter(p.second, p.second, w, g));
}

} // treedec

#endif //TD_NOBOOST_H

// vim:ts=8:sw=4:et
