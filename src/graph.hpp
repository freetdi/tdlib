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
 *    G& g, typename outedge_set<G>::type& bag)
 */

#ifndef TD_GRAPH_H
#define TD_GRAPH_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "degree.hpp"
#include "std.hpp"
#include "trace.hpp"

// OUCH
//#include "fill.hpp"

#ifndef TD_STRUCT_BAG
#define TD_STRUCT_BAG
struct bag{ //
    std::set<unsigned int> bag;
    typedef std::set<unsigned int> bagtype;
};
#endif

namespace treedec{ //
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
{untested();
    remove_vertex(*u, g);
}

template<typename vertex_descriptor>
struct vertex_callback{ //
    virtual ~vertex_callback(){};
    virtual void operator()(vertex_descriptor)=0;
};

template<typename G_t>
struct edge_callback{ //
    typedef typename boost::graph_traits<G_t>::edge_descriptor edge_descriptor;
    virtual ~edge_callback(){};
    virtual void operator()(edge_descriptor)=0;
};

template<typename G_t>
struct graph_callback{ // fixme: union of the above?
    typedef typename boost::graph_traits<G_t>::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    virtual ~graph_callback(){};
    virtual void operator()(edge_descriptor)=0;
    virtual void operator()(vertex_descriptor)=0;
};

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
{
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
    if(erase){untested();
        remove_vertex(v, g);
    }
}

} // noboost

namespace treedec{
// turn vertex range into clique.
// call cb on newly created edges and incident vertices.
// returns the number of newly created edges.
template<typename B, typename E, typename G_t>
size_t /*hmm*/ make_clique(B nIt1, E nEnd, G_t &G, treedec::graph_callback<G_t>* cb=NULL)
{itested();
    size_t counter=0;
    B nIt2;
    for(; nIt1!=nEnd; ++nIt1){
#if 1
        if(cb){untested();
            (*cb)(*nIt1); // hmm
        }
#endif
        nIt2 = nIt1;
        nIt2++;
        for(; nIt2 != nEnd; nIt2++){

            BOOST_AUTO(ep, boost::add_edge(*nIt1, *nIt2, G));
            if(ep.second){
               ++counter;

               if(cb){untested();
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

template<typename nIter_t, typename G_t>
void make_clique(nIter_t nIter, G_t &G, treedec::graph_callback<G_t>* cb=NULL)
{ itested();
    typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nEnd;
    boost::tie(nIt1, nEnd) = nIter;
    make_clique(nIt1, nEnd, G, cb);
}

// insert the neighbors of v in G into B
template<typename B_t, typename V_t, typename G_t>
void insert_neighbours(B_t &B, V_t v, G_t const &G)
{ itested();
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G);
    for(; nIt!=nEnd; ++nIt){itested();
        B.insert(*nIt);
    }
}

template<typename B_t, typename V_t, typename G_t>
void insert_neighbours(B_t &B, V_t v, V_t w, G_t const &G)
{
    insert_neighbours(B, v, G);
    insert_neighbours(B, w, G);
}

template<typename B_t, typename V_t, typename G_t>
void insert_neighbours(B_t &B, V_t v, V_t w, V_t x, G_t const &G)
{ untested();
    insert_neighbours(B, v, w, G);
    insert_neighbours(B, x, G);
}

// insert the neighbors of v in G into B
// B starts empty.
// equivalent to B = outedge_set(v, G) where applicable
template<typename B_t, typename V_t, typename G_t>
void assign_neighbours(B_t &B, V_t v, G_t const &G)
{
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G);
    for(; nIt!=nEnd; ++nIt){itested();
        B.insert(B.end(), *nIt);
    }
}

// insert neighbors of both vertices, after clearing B
template<typename B_t, typename V_t, typename G_t>
void assign_neighbours(B_t &B, V_t v, V_t w, V_t x, G_t const &G)
{
    B.clear();
    insert_neighbours(B, v, w, G);
    insert_neighbours(B, x, G);
}

#if 0 // unneeded wrappers
//transitional wrapper. don't use.
template<typename B, typename E, typename G_t>
size_t /*hmm*/ make_clique(B nIt1, E nEnd, G_t &G, treedec::graph_callback<G_t>* cb=NULL)
{
    return treedec::make_clique(nIt1, nEnd, G, cb);
}

// convenience wrapper (boost-style iterator pair)
// transitional. do not use
template<typename nIter_t, typename G_t>
void make_clique(nIter_t nIter, G_t &G, treedec::graph_callback<G_t>* cb=NULL)
{ itested();
    typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nEnd;
    boost::tie(nIt1, nEnd) = nIter;
    treedec::make_clique(nIt1, nEnd, G, cb);
}
#endif

#if 0
// FIXME: wrong name, used in preprocessing.
// FIXME: does not use G...
template<typename B_t, typename nIter_t, typename G_t>
// void paste_vertex_range
void fetch_neighbourhood(B_t &B, nIter_t nIter, G_t &G)
{ itested();
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = nIter; nIt != nEnd; nIt++){itested();
        // FIXME: hint?
        B.insert(*nIt);
    }
}
#endif

// FIXME: is this required?!
template <typename G_t>
inline typename boost::graph_traits<G_t>::vertex_descriptor
   get_min_degree_vertex(const G_t &G, bool ignore_isolated_vertices=false)
{untested();
    unsigned int min_degree = UINT_MAX;

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    boost::tie(vIt, vEnd) = boost::vertices(G);
    typename boost::graph_traits<G_t>::vertex_descriptor min_vertex = *vIt++;
    for(; vIt != vEnd; vIt++){untested();
        unsigned int degree = boost::degree(*vIt, G);
        if(degree <= min_degree){untested();
            if(ignore_isolated_vertices && degree == 0){ continue; }
            min_degree = degree;
            min_vertex = *vIt;
        }
    }

    return min_vertex;
}

template <typename G_t>
inline typename boost::graph_traits<G_t>::vertex_descriptor
   get_least_common_vertex(const typename boost::graph_traits<G_t>::vertex_descriptor &min_vertex,
           const G_t &G);

#if 0 // unused (hopefully)
// FIXME: it's make_clique.
template<typename G_t>
unsigned int eliminate_vertex(typename boost::graph_traits<G_t>::vertex_descriptor v, G_t &G)
{ unreachable(); // bogus function.
    noboost::make_clique(boost::adjacent_vertices(v, G), G);
    unsigned int deg = boost::degree(v, G);
    boost::clear_vertex(v, G);
    return deg; // does not make sense. the caller already knows the degree.
}
#endif

// FIXME: what does this do?
template <typename G_t>
inline void make_degree_sequence(const G_t &G,
          std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &degree_sequence)
{untested();
    unsigned int max_degree = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){untested();
        unsigned int degree = boost::degree(*vIt, G);
        max_degree = (degree>max_degree)? degree : max_degree;
    }

    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > buckets(max_degree+1);
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){untested();
        unsigned int degree = boost::degree(*vIt, G);
        if(degree > 0){untested();
            buckets[degree].push_back(*vIt);
        }
    }
    for(unsigned int i = 1; i <= max_degree; i++){untested();
        for(unsigned int j = 0; j < buckets[i].size(); j++){untested();
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
inline unsigned int get_pos(const typename boost::graph_traits<G_t>::vertex_descriptor v, const G_t& G){itested();
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
{untested();
    // base case: don't know. lets say yes (better in assertions)
    return true;
}

template<class G>
struct outedge_set{ //
    typedef std::set<unsigned> type;
//	typedef std::set type;
};

// kludge for balu
// TODO: move to graph_traits?!
template<class G>
struct treedec_chooser{ //
    typedef unsigned value_type;
    typedef std::set<unsigned> bag_type;
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, bag> type;
};

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

// FIXME: not part of noboost
template<class T>
struct treedec_traits{ //
    typedef typename detail::vdstuff<
       boost::is_same< typename T::vertex_property_type, bag >::value,
        typename T::vertex_property_type >::type vd_type;

    typedef typename detail::vdstuff<
       boost::is_same< typename T::vertex_property_type, bag >::value,
        typename T::vertex_property_type >::bag_type bag_type;

};

namespace detail{ //
template<class B, class T, class V>
struct tmpbaghack{ //
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
struct tmpbaghack<bag, T_t, V>{ //
    static typename treedec_traits<T_t>::bag_type& get_bag(T_t& t, V& v)
    {itested();
        return t[v].bag;
    }
    static typename treedec_traits<T_t>::bag_type const& get_bag(T_t const& t, V const& v)
    {itested();
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

// FIXME: move to graph_traits
template<class G_t>
struct deg_chooser{ //
    typedef typename misc::DEGS<G_t> type;
    typedef type degs_type; // transition? don't use.
};

// put neighbors of c into bag
// isolate c in g
// return bag
// optionally: pass pointer to bag for storage.
template<class G>
inline void detach_neighborhood(
        typename boost::graph_traits<G>::vertex_descriptor& c,
        G& g, typename outedge_set<G>::type& bag);


// count number of edges missing in 1-neighborhood of v
template <typename G_t>
inline size_t count_missing_edges(
        const typename boost::graph_traits<G_t>::vertex_descriptor v, G_t const &G)
{ itested();
    size_t missing_edges = 0;

    typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
    for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(v, G); nIt1 != nEnd; nIt1++){untested();
        nIt2 = nIt1;
        nIt2++;
        for(; nIt2 != nEnd; nIt2++){untested();
            if(!boost::edge(*nIt1, *nIt2, G).second){untested();
                ++missing_edges;
            }
        }
    }
    return missing_edges;
}

// collect neighbors of c into bag
// remove c from graph (i.e. isolate).
// turn subgraph induced by bag into clique.
//
// optionally:
// - call back on each newly created edge, and on each vertex incident to such
//   edge.
// - provide memory through bag pointer.
//
template<class G, class B>
size_t /*G::edge_index_type?*/ make_clique_and_detach(
        typename boost::graph_traits<G>::vertex_descriptor c,
        G& g, B& bag,
        treedec::graph_callback<G>* cb=NULL)
{
    detach_neighborhood(c, g, bag);
    return make_clique(bag.begin(), bag.end(), g, cb);
}

namespace detail{ //

    // iterate over edges adjacent to both v and s
    // implementation: iterate outedges(v), skip non-outedges(s).
    template<class G>
    class shared_adj_iter : public boost::graph_traits<G>::adjacency_iterator{ //
    public:
        typedef typename boost::graph_traits<G>::adjacency_iterator adjacency_iterator;
        typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
    public: // construct
        shared_adj_iter(adjacency_iterator v, adjacency_iterator ve,
                        vertex_descriptor s, G const& g)
            : adjacency_iterator(v), _ve(ve),
              _s(s), _g(g)
        {
            skip();
        }
        shared_adj_iter(vertex_descriptor v,
                        vertex_descriptor s, G const& g)
            : adjacency_iterator(boost::adjacent_vertices(v, g).first),
              _ve(boost::adjacent_vertices(v, g).second),
              _s(s), _g(g)
        {untested();
            skip();
        }
        shared_adj_iter(const shared_adj_iter& p)
            : adjacency_iterator(p), _ve(p._ve), _s(p._s), _g(p._g)
        {
        }
    public: //ops
        shared_adj_iter& operator++(){
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
                if(typename boost::graph_traits<G>::adjacency_iterator(*this)==_ve){
                    return;
                }else if(!boost::edge(**this, _s, _g).second){
                    adjacency_iterator::operator++();
                }else{
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
    inline common_out_edges(typename boost::graph_traits<G>::vertex_descriptor v,
                     typename boost::graph_traits<G>::vertex_descriptor w,
                     const G& g)
{ itested();
    typedef typename detail::shared_adj_iter<G> Iter;
    BOOST_AUTO(p, boost::adjacent_vertices(v, g));
    return std::make_pair(Iter(p.first, p.second, w, g),
                          Iter(p.second, p.second, w, g));
}

// check if two vertices have the same neighborhood.
// return false if not.
template<class G_t>
bool check_twins(typename boost::graph_traits<G_t>::vertex_descriptor v,
                 typename boost::graph_traits<G_t>::vertex_descriptor w,
                 const G_t& G)
{ itested();
    typename outedge_set<G_t>::type N1, N2;
    assign_neighbours(N1, v, G);
    assign_neighbours(N2, w, G);

    return(N1==N2);
}

} // treedec


// transition
namespace noboost{
    using treedec::bag;
    using treedec::check;
    using treedec::deg_chooser;
    using treedec::contract_edge;
    using treedec::edge_callback;
//    using treedec::eliminate_vertex;
//    using treedec::fetch_neighbourhood; // obsolete?
    using treedec::get_min_degree_vertex;
    using treedec::get_pos;
    using treedec::get_vd;
    using treedec::is_valid;
    using treedec::make_clique;
    using treedec::make_degree_sequence;
    using treedec::outedge_set;
    using treedec::remove_vertex;
    using treedec::treedec_traits;
    using treedec::treedec_chooser;
    using treedec::vertex_callback;
} // noboost

#endif // guard
// vim:ts=8:sw=4:et
