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
// unsigned int get_adhesion(T_t T)
//

#ifndef TD_MISC
#define TD_MISC

#define DEBUG

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include "simple_graph_algos.hpp"

namespace treedec{


/* checks if a tree decomposition is valid with respect to G
 * 
 *  0 = valid
 * -1 = (at least) not a tree
 * -2 = (at least) not all vertices covered (but a tree)
 * -3 = (at least) not all edges covered (but a tree and all vertices covered)
 * -4 = there exist vertices, coded in the bags of T, that are not connected in T 
 *                           (but T is a tree and all edges/vertices are covered)
 */
template <typename G_t, typename T_t>
int is_valid_treedecomposition(G_t &G, T_t T){
    //checks if T is a tree
    std::vector<int> component(boost::num_vertices(T));
    int num = boost::connected_components(T, &component[0]);
    if(num > 1 || boost::num_edges(T) > boost::num_vertices(T)-1){
#ifdef DEBUG
        std::cout << "decomposition is not a tree!" << std::endl;
#endif
        return -1;
    }

    //checks if exactly the vertices of G are covered
    std::set<unsigned int> coded_vertices;
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++)
        coded_vertices.insert(T[*tIt].bag.begin(), T[*tIt].bag.end());

    std::set<unsigned int> vertices;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        vertices.insert(G[*vIt].id);

    if(coded_vertices != vertices){
        return -2;
    }

    //checks if all edges are covered
    std::vector<std::set<unsigned int> > edges;
    
    for (boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        std::set<unsigned int> edge;
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
            edge.insert(G[*vIt].id);
            edge.insert(G[*nIt].id);
            edges.push_back(edge);
            edge.clear();
        }
    }
    for(std::vector<std::set<unsigned int> >::iterator it = edges.begin(); it != edges.end(); it++){
        bool isSubset = false;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){  
            if(std::includes(T[*tIt].bag.begin(), T[*tIt].bag.end(), it->begin(), it->end())){
                isSubset = true;
                break;
            }
        }
        if(!isSubset){
            //not all edges covered
#ifdef DEBUG
            std::cout << "not all edges covered!" << std::endl;
#endif
            return -3;
        }
    }
    std::set<unsigned int> forgotten;

    while(boost::num_vertices(T) != 1){
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            if(T[*tIt].bag.size() == 0 || boost::out_degree(*tIt, T) == 0){
                boost::clear_vertex(*tIt, T);
                boost::remove_vertex(*tIt, T);
                break;
            }
            if(boost::out_degree(*tIt, T) == 1){
                std::set<unsigned int> intersection;
                std::set_intersection(forgotten.begin(), forgotten.end(), T[*tIt].bag.begin(), T[*tIt].bag.end(), std::inserter(intersection, intersection.begin()));
                if(!intersection.empty()){
                    //there are coded vertices, that are not connected in T
#ifdef DEBUG
                    std::cout << "decomposition not connected!" << std::endl;
#endif
                    return -4;
                }
                
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

template <typename T_t>
unsigned int get_adhesion(T_t T){
    unsigned int max = 0;
    while(boost::num_vertices(T) > 1){
        typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            if(T[*tIt].bag.size() == 0 || boost::out_degree(*tIt, T) == 0){
                boost::clear_vertex(*tIt, T);
                boost::remove_vertex(*tIt, T);
                continue;
            }
            if(boost::out_degree(*tIt, T) == 1){
                typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                typename boost::graph_traits<T_t>::vertex_descriptor parent;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*tIt, T); nIt != nEnd; nIt++)
                    parent = *nIt;

                std::set<unsigned int> intersection;
                std::set_intersection(T[*tIt].bag.begin(), T[*tIt].bag.end(), T[parent].bag.begin(), T[parent].bag.end(), std::inserter(intersection, intersection.begin()));
                
                max = (intersection.size() > max)? intersection.size() : max;
                
                boost::clear_vertex(*tIt, T);
                boost::remove_vertex(*tIt, T);
                break;
            }
        }
    }
    return max;
}

template <typename T_t>
void make_small(T_t &T){
    while(true){
        bool modified = false;
        typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
#ifdef A
            if(T[*tIt].bag.size() == 0 || boost::out_degree(*tIt, T) == 0){
                boost::clear_vertex(*tIt, T);
                boost::remove_vertex(*tIt, T);
                continue;
            }
#endif
            if(boost::out_degree(*tIt, T) == 1){
                typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                boost::tie(nIt, nEnd) = boost::adjacent_vertices(*tIt, T);

                if(std::includes(T[*nIt].bag.begin(), T[*nIt].bag.end(), T[*tIt].bag.begin(), T[*tIt].bag.end())){
                    boost::clear_vertex(*tIt, T);
                    boost::remove_vertex(*tIt, T);
                    modified = true;
                    break;
                }
            }
            else if(boost::out_degree(*tIt, T) == 2){
                 typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                 std::vector<typename boost::graph_traits<T_t>::vertex_descriptor> N;
                 for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*tIt, T); nIt != nEnd; nIt++)
                     N.push_back(*nIt);

                 if(std::includes(T[N[0]].bag.begin(), T[N[0]].bag.end(), T[*tIt].bag.begin(), T[*tIt].bag.end())
                 || std::includes(T[N[1]].bag.begin(), T[N[1]].bag.end(), T[*tIt].bag.begin(), T[*tIt].bag.end())){
                     boost::add_edge(N[0], N[1], T);
                     boost::clear_vertex(*tIt, T);
                     boost::remove_vertex(*tIt, T);
                     modified = true;
                     break;
                 }
            }
        }
        if(!modified)
            return;
    }
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

//glues two "disjoint" decompositions (e.g. decompositions of two components of a graph)
template <typename T_t>
void glue_decompositions(T_t &T1, T_t &T2){
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;

    //copy T2 to T1 and add an edge from root to an arbitrary vertex of T2
    std::vector<typename boost::graph_traits<T_t>::vertex_descriptor> idxMap(boost::num_vertices(T2));
    std::map<typename boost::graph_traits<T_t>::vertex_descriptor, unsigned int> vertex_map;
    unsigned int id = 0;
    for(boost::tie(tIt, tEnd) = boost::vertices(T2); tIt != tEnd; tIt++){
        idxMap[id] = boost::add_vertex(T1);
        vertex_map.insert(std::pair<typename boost::graph_traits<T_t>::vertex_descriptor, unsigned int>(*tIt, id));
        T1[idxMap[id++]].bag = T2[*tIt].bag;
    }

    typename boost::graph_traits<T_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(T2); eIt != eEnd; eIt++){
        typename std::map<typename boost::graph_traits<T_t>::vertex_descriptor, unsigned int>::iterator v, w;
        v = vertex_map.find(boost::source(*eIt, T2));
        w = vertex_map.find(boost::target(*eIt, T2));

        boost::add_edge(idxMap[v->second], idxMap[w->second], T1);
    }

    typename boost::graph_traits<T_t>::vertex_iterator tIt2, tEnd2;

    boost::tie(tIt, tEnd) = boost::vertices(T1);

    boost::add_edge(*tIt, idxMap[0], T1);
}

template <typename G_t>
void reorder_ids_graph(G_t &G, std::vector<unsigned int> &id_map){
    id_map.resize(boost::num_vertices(G));
    unsigned int k = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        id_map[k] = G[*vIt].id;
        G[*vIt].id = k++; 
    }
}

template <typename T_t>
void reorder_ids_decomposition(T_t &T, std::vector<unsigned int> &id_map){
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        std::set<unsigned int> new_bag;
        for(std::set<unsigned int>::iterator sIt = T[*tIt].bag.begin(); sIt != T[*tIt].bag.end(); sIt++)
            new_bag.insert(id_map[*sIt]);

        T[*tIt].bag = new_bag;
    }
}


}

#endif
