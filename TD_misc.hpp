// Lukas Larisch, 2014 - 2015
//
// (c) 2014-2015 Goethe-Universität Frankfurt
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
// Offers miscellaneous algorithms about computing tree decompositions.
//
//
// These functions are most likely to be interesting for outside use:
//
// bool is_valid_decomposition(G_t &G, T_t T)
// void trivial_decomposition(G_t &G, T_t &T)
// int get_width(T_t &T)
// float get_average_bag_size(T_t &T)
//

#ifndef TD_MISC
#define TD_MISC

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include "simple_graph_algos.hpp"

namespace treedec{


/* checks if a tree decomposition is valid with respect to G
 * 
 *  0 = valid
 * -1 = (at least) not a tree
 * -2 = not all vertices or edges covered (but a tree)
 * -3 = there exist vertices, coded in the bags of T, that are not connected
 */
template <typename G_t, typename T_t>
int is_valid_decomposition(G_t &G, T_t T){
    //checks if T is a tree
    std::vector<int> component(boost::num_vertices(T));
    if(boost::connected_components(T, &component[0]) != 1 || boost::num_edges(T) != boost::num_vertices(T)-1)
        return -1;
    //checks if all edges are covered
    std::vector<std::set<unsigned int> > edges;
    typename boost::graph_traits<G_t>::vertex_iterator v, v_end;
    typename boost::graph_traits<G_t>::adjacency_iterator n, n_end;
    for (boost::tie(v, v_end) = boost::vertices(G); v != v_end; v++){
        std::set<unsigned int> edge;
        for(boost::tie(n, n_end) = boost::adjacent_vertices(*v, G); n != n_end; n++){
            edge.insert(G[*v].id);
            edge.insert(G[*n].id);
            edges.push_back(edge);
            edge.clear();
        }
    }
    for(std::vector<std::set<unsigned int> >::iterator it = edges.begin(); it != edges.end(); it++){
        typename boost::graph_traits<T_t>::vertex_iterator vIt, vEnd;
        bool isSubset = false;
        for(boost::tie(vIt, vEnd) = boost::vertices(T); vIt != vEnd; vIt++){  
            if(std::includes(T[*vIt].bag.begin(), T[*vIt].bag.end(), it->begin(), it->end())){
                isSubset = true;
                break;
            }
        }
        if(!isSubset)
            //not all vertices/edges covered, not covered vertices/edges can be found in *it
            return -2;
    }
    std::set<unsigned int> forgotten;

    while(boost::num_vertices(T) != 1){
        typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            if(T[*tIt].bag.size() == 0 || boost::out_degree(*tIt, T) == 0){
                boost::clear_vertex(*tIt, T);
                boost::remove_vertex(*tIt, T);
                break;
            }
            if(boost::out_degree(*tIt, T) == 1){
                if(std::includes(forgotten.begin(), forgotten.end(), T[*tIt].bag.begin(), T[*tIt].bag.end()))
                    //there are coded vertices, that are not connected in T, not connected vertices can be found in T[*tIt].bag
                    return -3;
                
                typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                typename boost::graph_traits<T_t>::vertex_descriptor parent;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*tIt, T); nIt != nEnd; nIt++){
                    parent = *nIt;
                }
                
                std::set_difference(T[*tIt].bag.begin(), T[*tIt].bag.end(), T[parent].bag.begin(), T[parent].bag.end(), std::inserter(forgotten, forgotten.begin()));

                
                boost::clear_vertex(*tIt, T);
                boost::remove_vertex(*tIt, T);
                break;
            }
        }
    }
    
    return 0;
}

template <typename G_t, typename T_t>
void trivial_decomposition(G_t &G, T_t &T){
    std::set<unsigned int> bag;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        bag.insert(G[*vIt].id);
    
    typename boost::graph_traits<T_t>::vertex_descriptor t;
    t = boost::add_vertex(T);
    T[t].bag = bag;
}

template <typename T_t>
int get_width(T_t &T){
    int max = -1;
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++)
        max = ((int)T[*tIt].bag.size() > max)? (int)T[*tIt].bag.size() : max;
    
    return (max-1);
} 

template <typename T_t>
float get_average_bag_size(T_t &T){
    float avg = 0.0;
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++)
        avg += T[*tIt].bag.size();
    
    return avg/boost::num_vertices(T);
} 


template <typename G_t>
void make_index_map(G_t &G, std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &idxMap){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    unsigned int max = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        max = (G[*vIt].id > max)? G[*vIt].id : max;
    
    idxMap.resize(max+1);
    
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        idxMap[G[*vIt].id] = *vIt; 
}


}

#endif
