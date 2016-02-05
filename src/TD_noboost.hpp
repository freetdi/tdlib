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
    template<typename G>
    using vertex_iterator = typename boost::graph_traits<G>::vertex_iterator;
    template<typename G>
    using vertex_descriptor = typename boost::graph_traits<G>::vertex_descriptor;
    template<typename G>
    using adjacency_iterator = typename boost::graph_traits<G>::adjacency_iterator;

    template<typename G>
    void remove_vertex(vertex_iterator<G> u, G &g)
    {
        remove_vertex(*u, g);
    }

    // vertex v will remain as isolated node
    template<typename G>
    void contract_edge(vertex_descriptor<G> v,
                       vertex_descriptor<G> w,
                       G &g, bool erase=true)
    {
        adjacency_iterator<G> I, E;
//        for(boost::tie(I, E)=boost::graph_traits<G>::adjacent_vertices(v, g); I!=E; ++I)
        for(boost::tie(I, E)=boost::adjacent_vertices(v, g); I!=E; ++I){
            if(*I != w){
                //boost::graph_traits<G>::add_edge(w, *I, g);
                boost::add_edge(w, *I, g);
            }
        }

        //boost::graph_traits<G>::clear_vertex(v, g);
        boost::clear_vertex(v, g);
    }

    // vertex v will remain as isolated node, unless erase.
    // if erase, v will be deleted (not collapsed).
    template<typename G>
    void contract_edge(vertex_iterator<G> v,
                       vertex_descriptor<G> into,
                       G &g, bool erase=true)
    {
        contract_edge(*v, into, g);
        if(erase){
            noboost::remove_vertex(v, g);
        }
    }
} // namespace noboost

#endif
// vim:ts=8:sw=4:noet
