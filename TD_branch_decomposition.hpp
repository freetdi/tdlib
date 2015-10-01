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

#ifndef TD_BRANCH_DECOMPOSITION
#define TD_BRANCH_DECOMPOSITION

#include <boost/graph/adjacency_list.hpp>
#include "simple_graph_algos.hpp"

namespace treedec{


/* checks if a branch decomposition is valid with respect to G
 * 
 *  0 = valid
 * -1 = (at least) not a tree
 * -2 = T is not cubic
 * -3 = the bags of the leafs of T are not exactly the edges of G (but T is a tree)
 * -4 = there exists a bag of size 1
 */
template <typename G_t, typename T_t>
int is_valid_branchdecomposition(G_t &G, T_t T){
    if(boost::num_edges(G) == 0 && boost::num_vertices(T) == 0)
        return 0;

    //checks if T is a tree
    std::vector<int> component(boost::num_vertices(T));
    if(boost::connected_components(T, &component[0]) != 1 || boost::num_edges(T) != boost::num_vertices(T)-1)
        return -1;

    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;

    //T cubic?
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        int degree = boost::out_degree(*tIt, T);
        if(degree == 3 || degree == 1 || degree == 0)
            continue;
        else{
            return -2;
        }
    }

    //collect the edges of G
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    std::vector<std::set<unsigned int> > edges;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        std::set<unsigned int> edge;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
            edge.insert(G[*vIt].id);
            edge.insert(G[*nIt].id);
            edges.push_back(edge);
            edge.clear();
        }
    }

    //tests whether the bags of the leafs of T are exactly the edges of G
    for(std::vector<std::set<unsigned int> >::iterator it = edges.begin(); it != edges.end(); it++){
        typename boost::graph_traits<T_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<T_t>::vertex_descriptor t_node;
        unsigned int covered_count = 0;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){  
            if(std::includes(T[*tIt].bag.begin(), T[*tIt].bag.end(), it->begin(), it->end())){
                t_node = *tIt;
                if(*it == T[*tIt].bag)
                    covered_count++;
            }
        }
        if(covered_count != 1){
            return -3;
        }
    }

    //exist bags of size 1?
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){  
        if(T[*tIt].bag.size() == 1)
            return -4;
        
    }

    return 0;
}

template <typename T_t>
void explore_component(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor t, std::set<unsigned int> &V, typename boost::graph_traits<T_t>::vertex_descriptor &parent){
    std::vector<typename boost::graph_traits<T_t>::vertex_descriptor > N;

    typename boost::graph_traits<T_t>::adjacency_iterator  nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(t, T); nIt != nEnd; nIt++){
        if(*nIt != parent)
            N.push_back(*nIt);
    }

    for(std::set<unsigned int>::iterator sIt = T[t].bag.begin(); sIt != T[t].bag.end(); sIt++)
        V.insert(*sIt);

    for(unsigned int i = 0; i < N.size(); i++)
        explore_component(T, N[i], V, t);
}

template <typename T_t>
void compute_cutset(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor t, std::set<unsigned int> &V){
    std::vector<typename boost::graph_traits<T_t>::vertex_descriptor > N;

    typename boost::graph_traits<T_t>::adjacency_iterator  nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(t, T); nIt != nEnd; nIt++)
        N.push_back(*nIt);

    std::set<unsigned int> T1;
    explore_component(T, N[0], T1, t);

    std::set<unsigned int> T2;
    explore_component(T, N[1], T2, t);

    std::set<unsigned int> T3;
    explore_component(T, N[2], T3, t);

    std::set<unsigned int> V1, V2, V3, V_;
    std::set_intersection(T1.begin(), T1.end(), T2.begin(), T2.end(), std::inserter(V1, V1.begin()));
    std::set_intersection(T1.begin(), T1.end(), T3.begin(), T3.end(), std::inserter(V2, V2.begin()));
    std::set_intersection(T2.begin(), T2.end(), T3.begin(), T3.end(), std::inserter(V3, V3.begin()));

    std::set_union(V1.begin(), V1.end(), V2.begin(), V2.end(), std::inserter(V_, V_.begin()));
    std::set_union(V_.begin(), V_.end(), V3.begin(), V3.end(), std::inserter(V, V.begin()));
}


template <typename G_t, typename T_t>
void branch_to_tree_decomposition(G_t &G, T_t &T){
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;

    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        if(boost::out_degree(*tIt, T) == 3){
            std::set<unsigned int> V; 
            compute_cutset(T, *tIt, V);
            T[*tIt].bag = V;
        }
    } 

    //glue isolated vertices with T
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) == 0){
            typename boost::graph_traits<T_t>::vertex_descriptor new_t_node = boost::add_vertex(T);
            std::set<unsigned int> bag;
            bag.insert(G[*vIt].id);
            T[new_t_node].bag = bag;

            if(boost::num_vertices(T) > 1){
                boost::tie(tIt, tEnd) = boost::vertices(T);
                boost::add_edge(*tIt, new_t_node, T);
            }
        }
    }
}


template <typename G_t, typename T_t>
void tree_to_branch_decomposition(G_t &G, T_t &T){
    if(boost::num_edges(G) == 0){
        T.clear();
        return;
    }

    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;

    //collect the edges of G
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    std::vector<std::set<unsigned int> > edges;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        std::set<unsigned int> edge;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
            edge.insert(G[*vIt].id);
            edge.insert(G[*nIt].id);
            edges.push_back(edge);
            edge.clear();
        }
    }

    for(std::vector<std::set<unsigned int> >::iterator it = edges.begin(); it != edges.end(); it++){
        typename boost::graph_traits<T_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<T_t>::vertex_descriptor t_node;
        unsigned int covered_count = 0;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){  
            if(std::includes(T[*tIt].bag.begin(), T[*tIt].bag.end(), it->begin(), it->end())){
                t_node = *tIt;
                if(boost::out_degree(*tIt, T) <= 1 && *it == T[*tIt].bag)
                    covered_count++;
            }
        }
        if(covered_count == 0){
            typename boost::graph_traits<T_t>::vertex_descriptor new_t_node = boost::add_vertex(T);
            T[new_t_node].bag = *it;
            boost::add_edge(t_node, new_t_node, T); 
        }
        else{
            for(unsigned int j = 0; j < covered_count-1; j++){
                for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){  
                    if(*it == T[*tIt].bag){
                        boost::clear_vertex(*tIt, T);
                        boost::remove_vertex(*tIt, T);
                    }
                }
            }
        }
    }

    //remove bags of size 1
    bool exists_singleton = true;
    while(exists_singleton){
        exists_singleton = false;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){  
            if(T[*tIt].bag.size() == 1){
                boost::clear_vertex(*tIt, T);
                boost::remove_vertex(*tIt, T);
                exists_singleton = true;
            }
        }
    }

    //make T cubic

    //trivial cubic tree
    if(boost::num_vertices(T) == 1)
        return;

    bool is_cubic = false;
    while(!is_cubic){
        is_cubic = true;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            int degree = boost::out_degree(*tIt, T);
            if(degree == 3 || degree == 1)
                continue;
            
            std::vector<typename boost::graph_traits<T_t>::vertex_descriptor > N;

            typename boost::graph_traits<T_t>::adjacency_iterator  nIt, nEnd;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*tIt, T); nIt != nEnd; nIt++)
                N.push_back(*nIt);

            if(degree == 2){
                boost::add_edge(N[0], N[1], T);
                boost::clear_vertex(*tIt, T);
                boost::remove_vertex(*tIt, T);
                is_cubic = false;
                break;
            }

            //degree > 3
            typename boost::graph_traits<T_t>::vertex_descriptor new_t_node = boost::add_vertex(T);
            T[new_t_node].bag = T[*tIt].bag;
            boost::add_edge(*tIt, new_t_node, T); 
            boost::remove_edge(*tIt, N[0], T); 
            boost::remove_edge(*tIt, N[1], T); 
            boost::add_edge(new_t_node, N[0], T); 
            boost::add_edge(new_t_node, N[1], T);
            is_cubic = false;
            break;
        }
    }

    //remove "interior" bags
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        if(boost::out_degree(*tIt, T) != 1)
            T[*tIt].bag.clear();
    } 
}

}

#endif
