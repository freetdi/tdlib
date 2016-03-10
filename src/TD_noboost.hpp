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
                   bool /*erase*/=true,
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
                   bool erase=true,
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


template<typename G>
inline unsigned get_id(const G& g, const vertex_descriptor_G& v )
{
    // works with "TD_graph_t" (augmented adj_list)
    return g[v].id;
}

// return the internal vertex position.
// to be used as a narrower alternative to vertex_descriptor.
// positions are in {0, 1, ..., num_vertices-1}, where applicable
template<typename G>
inline unsigned get_pos(const G& /*g*/, const vertex_descriptor_G& /*v*/ )
{
    assert(false); // incomplete.
    return 0;
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

// is this really part of noboost?
template<class T>
struct treedec_traits{
    typedef unsigned vd_type;
    typedef std::set<vd_type> bag_type;
};

template<typename T>
inline typename treedec_traits<T>::bag_type& bag(T& t,
	const typename boost::graph_traits<T>::vertex_descriptor& v)
{
    return t[v].bag;
}

} // namespace noboost

namespace detail{

    // DRAFT. no useful interface
template<typename G>
struct degree_mod : public noboost::vertex_callback<G>{
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
    degree_mod(misc::DEGS<G>* d, G* g) : _degs(d), _g(g){}

    // reinsert with degree-1
    void operator()(vertex_descriptor v){ untested();
        unsigned deg=boost::degree(v,*_g);
        (void)deg;
        assert(deg);
        bool done=(*_degs)[deg-1].insert(v).second;
        (void)done;
        assert(done);
    }
    //private: not yet.
//    BUG:: hardcoded type
        misc::DEGS<G>* _degs;
    private:
        degree_mod(const degree_mod&){}
        G* _g;
};

} //namespace detail

#endif //TD_NOBOOST_H

// vim:ts=8:sw=4:et
