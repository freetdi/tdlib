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
 * Implements functions on generic BGL graphs that can be implemented (much)
 * faster for particular graph classes.
 */

#ifndef TD_NOBOOST_H
#define TD_NOBOOST_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "TD_degree.hpp"

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

template<typename G>
struct vertex_callback{
    virtual ~vertex_callback(){};
    virtual void operator()(vertex_descriptor_G)=0;
};

//Vertex v will remain as isolated node.
//Calls cb on neighbors if degree drops by one,
//before dg will drop.
template<typename G>
void contract_edge(vertex_descriptor_G v,
                   vertex_descriptor_G target,
                   G &g,
                   bool /*erase*/=false,
                   vertex_callback<G>* cb=NULL)
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
                   vertex_callback<G>* cb=NULL)
{
    contract_edge(*v, into, g, erase, cb);
    if(erase){
        noboost::remove_vertex(v, g);
    }
}

template<typename nIter_t, typename G_t>
void make_clique(nIter_t nIter, G_t &G){
    typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
    for(boost::tie(nIt1, nEnd) = nIter; nIt1 != nEnd; nIt1++){
        nIt2 = nIt1;
        nIt2++;
        for(; nIt2 != nEnd; nIt2++){
            boost::add_edge(*nIt1, *nIt2, G);
        }
    }
}

template<typename B_t, typename nIter_t, typename G_t>
void fetch_neighbourhood(B_t &B, nIter_t nIter, G_t &G){
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
inline unsigned get_pos(const typename boost::graph_traits<G_t>::vertex_descriptor v, const G_t& G){
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

// TRANSITIONAL. (remove later)
template<typename T_t>
inline typename treedec_traits<T_t>::bag_type& bag_(T_t& T,
	 const typename boost::graph_traits<T_t>::vertex_descriptor& v)
{
    return bag(v,T);
}


template<typename T_t>
inline typename treedec_traits<T_t>::bag_type const& bag_(T_t const& T,
       const typename boost::graph_traits<T_t>::vertex_descriptor& v)
{
    return bag(v,T);
}
// end TRANSITIONAL



} // namespace noboost


#endif //TD_NOBOOST_H

// vim:ts=8:sw=4:et
