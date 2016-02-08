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

/* implements functions on generic BGL graphs that can be implemented (much)
 * faster for particular graph classes.
 */

#ifndef TD_NOBOOST_H
#define TD_NOBOOST_H

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

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

    // vertex v will remain as isolated node
    // calls cb on neighbors if degree drops by one,
    //          before dg will drop.
    template<typename G>
    void contract_edge(vertex_descriptor_G v,
                       vertex_descriptor_G target,
                       G &g,
                       bool erase=true,
                       vertex_callback<G>* cb=NULL)
    {
        adjacency_iterator_G I, E;
        for(boost::tie(I, E)=boost::adjacent_vertices(v, g); I!=E; ++I){
            assert(boost::edge(v, *I, g).second);
            if(*I != target){
                bool added=boost::add_edge(target, *I, g).second;
        if(added){
            // rebasing edge from I-v to I-target.
        }else if(cb){
            // did not add, degree will drop by one.
            (*cb)(*I);
        }
            }else{
        }
        }

        //boost::graph_traits<G>::clear_vertex(v, g);
        boost::clear_vertex(v, g);
    assert(!boost::edge(v, target, g).second);
    }

    // vertex v will remain as isolated node, unless erase.
    // if erase, v will be deleted (not collapsed).
    template<typename G>
    void contract_edge(vertex_iterator_G v,
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
} // namespace noboost

#endif
