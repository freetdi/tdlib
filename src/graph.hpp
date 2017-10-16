// Felix Salfelder, 2015 - 2016
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

#include <boost/graph/adjacency_list.hpp>

#include "error.hpp"
#include "graph_traits.hpp"
#include "treedec_traits.hpp"
#include "platform.hpp"
#include "random_generators.hpp"
#include "trace.hpp"
#include "error.hpp"
#include "bits/bool.hpp"

// not yet
// #include <gala/boost.h>
// #include "graph_gala.hpp"

#include "degree.hpp"

#include "generic_elimination_search_overlay.hpp"
#include "induced_subgraph.hpp"

namespace treedec{

#define vertex_iterator_G typename boost::graph_traits<G>::vertex_iterator
#define vertex_descriptor_G typename boost::graph_traits<G>::vertex_descriptor
#define adjacency_iterator_G typename boost::graph_traits<G>::adjacency_iterator


// obsolete forward? better don't use
#if 0
template<typename G_t>
inline typename boost::graph_traits<G_t>::vertices_size_type
   get_pos(typename boost::graph_traits<G_t>::vertex_descriptor v, G_t const& G)
{
    return boost::get(boost::vertex_index, G, v);
}
#endif

//Vertex v will remain as isolated node.
//Calls cb on neighbors if degree drops by one,
//before dg will drop.
// TODO: pass edge or edge iterator.
template<typename G>
void contract_edge(vertex_descriptor_G v,
                   vertex_descriptor_G target,
                   G &g,
                   vertex_callback<vertex_descriptor_G>* cb=NULL)
{
    adjacency_iterator_G I, E;
    for(boost::tie(I, E)=boost::adjacent_vertices(v, g); I!=E; ++I){
        assert(boost::edge(v, *I, g).second);
        if(*I != target){
            bool added=treedec::add_edge(target, *I, g).second;
            if(added){
                //rebasing edge from I-v to I-target.
            }
            else if(cb){
                //did not add, degree will drop by one.
                (*cb)(*I);
            }
        }
    }

    boost::clear_vertex(v, g);
    assert(!boost::edge(v, target, g).second);
}


// turn vertex range into clique.
// call cb on newly created edges and incident vertices.
// (no longer?!) returns the number of newly created edges.
template<typename B, typename E, typename G_t>
void make_clique(B nIt1, E nEnd, G_t &G, typename treedec::graph_callback<G_t>* cb=NULL)
{
    typedef typename boost::graph_traits<G_t>::edge_descriptor edge_descriptor;
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
            std::pair<edge_descriptor, bool> ep=boost::add_edge(*nIt1, *nIt2, G);
            if(ep.second){
               if(cb){
                   (*cb)(*nIt1, *nIt2);
               }
            }
            if(!boost::is_directed(G)){
            }
            else if( boost::edge(*nIt2, *nIt1, G).second ){
                // change later.
            }
            else{
                boost::add_edge(*nIt2, *nIt1, G);
            }
        }
    }
}

template<typename nIter_t, typename G_t>
void make_clique(nIter_t nIter, G_t &G, treedec::graph_callback<G_t>* cb=NULL)
{
    typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nEnd;
    boost::tie(nIt1, nEnd) = nIter;
    make_clique(nIt1, nEnd, G, cb);
}

// insert the neighbors of v in G into B
template<typename B_t, typename V_t, typename G_t>
void insert_neighbours(B_t &B, V_t v, G_t const &G)
{
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G);
    for(; nIt!=nEnd; ++nIt){
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
{
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
    for(; nIt!=nEnd; ++nIt){
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

// find a minimum degree vertex using linear search.
template <typename G_t>
inline typename boost::graph_traits<G_t>::vertex_descriptor
   get_min_degree_vertex(const G_t &G, bool ignore_isolated_vertices=false)
{
    unsigned int min_degree = UINT_MAX;

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    boost::tie(vIt, vEnd) = boost::vertices(G);
    typename boost::graph_traits<G_t>::vertex_descriptor min_vertex = *vIt++;
    for(; vIt != vEnd; vIt++){
        unsigned int degree = boost::out_degree(*vIt, G);
        if(degree <= min_degree){
            if(ignore_isolated_vertices && degree == 0){ continue; }
            min_degree = degree;
            min_vertex = *vIt;
        }
    }

    return min_vertex;
}

template <typename G_t, class M>
inline typename boost::graph_traits<G_t>::vertex_descriptor
   get_least_common_vertex(const typename boost::graph_traits<G_t>::vertex_descriptor min_vertex,
           M& marker,
           const G_t &G);

//Copy vertices of G into degree_sequence, ordered by degree, starting with lowest.
template <typename G_t>
inline void make_degree_sequence(const G_t &G,
          std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &degree_sequence)
{
    unsigned int max_degree = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::out_degree(*vIt, G);
        max_degree = (degree>max_degree)? degree : max_degree;
    }

    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > buckets(max_degree+1);
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::out_degree(*vIt, G);
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

// obsolete?
template<class G>
struct outedge_set;

template <typename G_t>
std::pair<typename boost::graph_traits<typename graph_traits<G_t>::directed_overlay>::vertex_descriptor,
          typename boost::graph_traits<typename graph_traits<G_t>::directed_overlay>::vertex_descriptor>
    make_digraph_with_source_and_sink(G_t const &G, std::vector<BOOL> const &disabled,
                 unsigned num_dis,
                 typename graph_traits<G_t>::directed_overlay& dg,
                 std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &idxMap,
                 typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> const &X,
                 typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> const &Y);


template<class G_t>
struct deg_chooser { //
    typedef typename misc::DEGS<G_t> type;
    typedef type degs_type; // transition? don't use.
};

// put neighbors of c into bag
// isolate c in g
// return bag
// optionally: pass pointer to bag for storage.

//TODO: not here.
template<class G_t, class B_t>
inline void detach_neighborhood(
        typename boost::graph_traits<G_t>::vertex_descriptor c,
        G_t& g, typename std::set<B_t> &bag)
{
    assert(boost::is_undirected(g));

    typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
    // inefficient.

    for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(c, g); nIt1 != nEnd; nIt1++)
    {
        bag.insert(get_vd(g, *nIt1));
    }
    boost::clear_vertex(c, g);
}

//TODO: not here.
template<class G_t, class B_t>
inline void detach_neighborhood(
        typename boost::graph_traits<G_t>::vertex_descriptor c,
        G_t& g, typename std::vector<B_t> &bag)
{
    assert(boost::is_undirected(g));

    typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
    // inefficient.

    unsigned i = 0;
    for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(c, g); nIt1 != nEnd; nIt1++)
    {
        bag[i++] = get_vd(g, *nIt1);
    }
    boost::clear_vertex(c, g);
}


// count number of edges missing in 1-neighborhood of v
// OBSOLETE
// OBSOLETE
// OBSOLETE
template <typename G_t>
inline size_t count_missing_edges(
        const typename boost::graph_traits<G_t>::vertex_descriptor v, G_t const &G)
{
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
void make_clique_and_detach(
        typename boost::graph_traits<G>::vertex_descriptor c,
        G& g, B& bag,
        treedec::graph_callback<G>* cb=NULL)
{
    detach_neighborhood(c, g, bag);
    make_clique(bag.begin(), bag.end(), g, cb);
}

namespace detail{ //

    // iterate over edges adjacent to both v and s
    // implementation: iterate outedges(v), skip non-outedges(s).
    template<class G>
    class shared_adj_iter : public boost::graph_traits<G>::adjacency_iterator{ //
    public:
        typedef typename boost::graph_traits<G>::adjacency_iterator adjacency_iterator;
        typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
        typedef typename boost::graph_traits<G>::adjacency_iterator baseclass;
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
                if(baseclass(*this)==_ve){
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

// inefficient.
// possibly obsolete.
template<class G>
std::pair<detail::shared_adj_iter<G>, detail::shared_adj_iter<G> >
inline common_out_edges(typename boost::graph_traits<G>::vertex_descriptor v,
        typename boost::graph_traits<G>::vertex_descriptor w,
        const G& g)
{
    typedef typename detail::shared_adj_iter<G> Iter;

    auto p=boost::adjacent_vertices(v, g);
    return std::make_pair(Iter(p.first, p.second, w, g),
                          Iter(p.second, p.second, w, g));
}

// check if two vertices have the same neighborhood.
// return false if not.
template<class G_t>
bool check_twins(typename boost::graph_traits<G_t>::vertex_descriptor v,
                 typename boost::graph_traits<G_t>::vertex_descriptor w,
                 const G_t& G)
{
    typename graph_traits<G_t>::outedge_set_type N1, N2;
    assign_neighbours(N1, v, G);
    assign_neighbours(N2, w, G);

    return(N1==N2);
}

// create ig.
// edge(v, w, ig) == edge(ordering[v], ordering[w], g)
// TODO: inverse map?
// TODO: partial/induced graph?
template<class G>
void immutable_clone(G const &g, typename graph_traits<G>::immutable_type& ig,
   std::vector<typename boost::graph_traits<G>::vertex_descriptor>* ordering=NULL)
{
    typedef typename graph_traits<G>::immutable_type immutable_type;
//    typedef typename graph_traits<G>::vertex_iterator vertex_iterator;
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G>::edge_iterator edge_iterator;

    if(!ordering){ incomplete();
    }
    assert(ordering->size()==boost::num_vertices(g)); // for now.


    if(boost::num_vertices(g) != boost::num_vertices(ig)){ untested();
        // drop a new one... (for now?)
        // inefficient.
        ig = MOVE(immutable_type(boost::num_vertices(g)));
    }else if(boost::num_edges(ig)){ untested();
        for(unsigned i=0; i<boost::num_vertices(g); ++i){
            boost::clear_vertex(i, ig);
        }
    }

    std::vector<unsigned> inverse_ordering(boost::num_vertices(g));
    for(unsigned o=0; o<boost::num_vertices(g); ++o){
        assert((*ordering)[o] < inverse_ordering.size());
        inverse_ordering[(*ordering)[o]] = o;
    }

    for(unsigned o=0; o<boost::num_vertices(g); ++o){
        assert((*ordering)[inverse_ordering[o]] == o);
    }

    edge_iterator e, eend;
    boost::tie(e,eend) = boost::edges(g);
    for(;e!=eend; ++e){
        vertex_descriptor s = boost::source(*e, g);
        vertex_descriptor t = boost::target(*e, g);
        assert(inverse_ordering[s] < boost::num_vertices(ig));
        assert(inverse_ordering[t] < boost::num_vertices(ig));
        assert((s==t) == (inverse_ordering[s]==inverse_ordering[t]));
        boost::add_edge(inverse_ordering[s], inverse_ordering[t], ig);
    }
}

// check if {*vd1, *vd2} is a subset of a bag adjacent to t.
template<class VD_t, class T_t>
class is_in_neighbour_bd{ //
public:
    is_in_neighbour_bd(T_t const& T,
        typename boost::graph_traits<T_t>::vertex_descriptor t)
       : _T(T), _t(t)
    {
    }
public:
    bool operator() (VD_t vd1, VD_t vd2)
    {
        assert(vd1!=vd2);
        if(vd1<vd2){
        }

        typedef typename boost::graph_traits<T_t>::adjacency_iterator bag_iterator;
        bag_iterator nIt, nEnd;
        boost::tie(nIt, nEnd) = boost::adjacent_vertices(_t, _T);
        for(; nIt != nEnd; nIt++){
            BOOST_AUTO(const& ibag, bag(*nIt, _T));

            // BUG, does not work on vectors.
            if(ibag.find(vd1)==ibag.end()){
            }else if(ibag.find(vd2)==ibag.end()){
            }else{
                return true;
            }
        }
        // a = vd1;
        // b = vd2;
        return false;
    }
private: //data
    T_t const &_T;
    typename boost::graph_traits<T_t>::vertex_descriptor _t;
public: // HACK
    long unsigned a, b;
};

// clone subgraph induced by bag into ig.
// store map V(H) -> X \subset V(G) in vdMap
// add more edges for MSVS, TODO: implement properly (how?)
template<class G_t, class T_t, class IG_t, class M_t>
inline void induced_subgraph_with_extra_edges
(G_t const &G, IG_t& ig, // typename graph_traits<G>::immutable_type& ig,
       T_t const& T, typename boost::graph_traits<T_t>::vertex_descriptor bd,
    //   URGHS. no default types without c++11.
     M_t* vdMap /*=NULL*/);


// insert reverse edges.
// dont know how to express "orientation" yet
template<class G>
void make_symmetric(G&g, bool force_oriented=false)
{ untested();
    // assert(g.is_directed);
    auto EE=boost::edges(g);
    if ( /*g.oriented ||*/ force_oriented) { untested();
        for(;EE.first!=EE.second; ++EE.first){
            auto s=boost::source(*EE.first, g);
            auto t=boost::target(*EE.first, g);
            boost::add_edge(t,s,g);
        }
    }else{ untested();
        for(;EE.first!=EE.second; ++EE.first){
            auto s=boost::source(*EE.first, g);
            auto t=boost::target(*EE.first, g);
            if(!boost::edge(s,t,g).second){ untested();
                unreachable(); //?
                boost::add_edge(s,t,g);
            }else if(!boost::edge(t,s,g).second){ untested();
                boost::add_edge(t,s,g);
            }
        }
    }
}

namespace detail{

template<class G, class X=void>
struct edge_helper{
};

template<class G>
struct edge_helper<G, typename std::enable_if< std::is_convertible<
                    typename boost::graph_traits<G>::directed_category*, boost::undirected_tag*
                    >::value, void>::type > {
    typedef typename boost::graph_traits<G>::edges_size_type size_type;
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;

    static size_type num(G const& g){
        return boost::num_edges(g);
    }
    static std::pair<typename boost::graph_traits<G>::edge_descriptor, bool>
        add(vertex_descriptor x, vertex_descriptor y, G& g){

	return boost::add_edge(x, y, g);
    }
};

template<class G>
struct edge_helper<G, typename std::enable_if< std::is_convertible<
                    typename boost::graph_traits<G>::directed_category*, boost::directed_tag*
                    >::value, void>::type > {
    typedef typename boost::graph_traits<G>::edges_size_type size_type;
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;

    typedef typename std::enable_if< std::is_same< boost::directed_tag,
                                typename boost::graph_traits<G>::directed_category >::value, G >::type T;

    static size_type num(G const& g){
        assert(1 ^ boost::num_edges(g)) ;
        return boost::num_edges(g)/2;
    }
    static std::pair<typename boost::graph_traits<G>::edge_descriptor, bool>
        add(vertex_descriptor x, vertex_descriptor y, G& g){

	boost::add_edge(y, x, g);
	return boost::add_edge(x, y, g);
    }
};

}

template<class G>
inline typename boost::graph_traits<G>::edges_size_type
  num_edges(G const& g)
{
    return detail::edge_helper<G>::num(g);
}

template<class G>
inline std::pair<typename boost::graph_traits<G>::edge_descriptor, bool>
add_edge(typename boost::graph_traits<G>::vertex_descriptor x,
		   typename boost::graph_traits<G>::vertex_descriptor y, G& g)
{
    return detail::edge_helper<G>::add(x, y, g);
}

template<class S, class G>
void open_neighbourhood(S& s, G const& g)
{
    graph_helper<G>::open_neighbourhood(s, g);
}

template<class S, class G>
void close_neighbourhood(S& s, G const& g)
{
    graph_helper<G>::close_neighbourhood(s, g);
}

template<class S, class G>
void saturate(S& s, G const& g)
{
    graph_helper<G>::saturate(s, g);
}

} // treedec

#endif // guard
// vim:ts=8:sw=4:et
