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
#include "graph.hpp"

namespace treedec{


/* Checks if a branch decomposition is valid with respect to G.
 *
 *  0 = valid
 * -1 = (at least) not a tree
 * -2 = T is not cubic
 * -3 = the bags of the leafs of T are not exactly the edges of G (but T is a tree)
 */
template <typename G_t, typename T_t>
int is_valid_branchdecomposition(G_t &G, T_t &T)
{
    typedef std::vector<typename treedec_traits<T_t>::bag_type> edges_vt;
    typedef typename edges_vt::iterator eIter;

    if(boost::num_edges(G) == 0 && boost::num_vertices(T) == 0){
        return 0;
    }

    //Checks if T is a tree.
    std::vector<int> component(boost::num_vertices(T));
    if(boost::connected_components(T, &component[0]) != 1 || boost::num_edges(T) != boost::num_vertices(T)-1){
        return -1;
    }

    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;

    //Checks if T is cubic.
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        int degree = boost::degree(*tIt, T);
        if(degree == 3 || degree == 1 || degree == 0){
            continue;
        }
        else{
            return -2;
        }
    }

    //Collect the edges of G.
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    std::vector<typename treedec_traits<T_t>::bag_type> edges(boost::num_edges(G));
    unsigned int i = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        typename treedec_traits<T_t>::bag_type edge;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
            if(*vIt >= *nIt){ continue; }

            edge.insert(*vIt);
            edge.insert(*nIt);
            edges[i++] = MOVE(edge);
            edge.clear();
        }
    }

    std::vector<bool> visited(boost::num_vertices(G), false);

    //Checks if the bags of the leafs of T are exactly the edges of G.
    for(eIter it = edges.begin(); it != edges.end(); it++){
        typename boost::graph_traits<T_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<T_t>::vertex_descriptor t_node;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            if(*it == bag(*tIt, T)){
                visited[*tIt] = true;
                break;
            }
        }
    }

    if(std::find(visited.begin(), visited.end(), false) != visited.end()){
        return -3;
    }

    return 0;
}

template <typename T_t>
void explore_component(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor t,
                typename treedec_traits<T_t>::bag_type &V,
                typename boost::graph_traits<T_t>::vertex_descriptor &parent)
{
    V.insert(bag(t, T).begin(), bag(t, T).end());

    typename boost::graph_traits<T_t>::adjacency_iterator  nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(t, T); nIt != nEnd; nIt++){
        if(*nIt != parent){
            explore_component(T, *nIt, V, t);
        }
    }
}

template <typename T_t>
void compute_cutset(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor t,
                typename treedec_traits<T_t>::bag_type &V)
{
    std::vector<typename boost::graph_traits<T_t>::vertex_descriptor > N(3);

    typename boost::graph_traits<T_t>::adjacency_iterator  nIt, nEnd;
    unsigned int i = 0;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(t, T); nIt != nEnd; nIt++){
        N[i++] = *nIt;
    }

    typename treedec_traits<T_t>::bag_type T1;
    explore_component(T, N[0], T1, t);

    typename treedec_traits<T_t>::bag_type T2;
    explore_component(T, N[1], T2, t);

    typename treedec_traits<T_t>::bag_type T3;
    explore_component(T, N[2], T3, t);

    typename treedec_traits<T_t>::bag_type V1, V2, V3, V_;
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
        if(boost::degree(*tIt, T) == 3){
            typename treedec_traits<T_t>::bag_type V;
            compute_cutset(T, *tIt, V);
            bag(*tIt, T) = MOVE(V);
        }
    }

    //Glue isolated vertices with T.
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::degree(*vIt, G) == 0){
            typename boost::graph_traits<T_t>::vertex_descriptor new_t_node = boost::add_vertex(T);
            bag(new_t_node, T).insert(*vIt);

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

    //Collect the edges of G.
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    typedef std::vector<typename treedec_traits<T_t>::bag_type> edges_vt;
    typedef typename edges_vt::iterator eIter;
    edges_vt edges(boost::num_edges(G));
    unsigned int c = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        typename treedec_traits<T_t>::bag_type edge;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
            if(*vIt >= *nIt){ continue; }

            edge.insert(*vIt);
            edge.insert(*nIt);

            edges[c++] = MOVE(edge);
            edge.clear();
        }
    }

    for(eIter it = edges.begin(); it != edges.end(); it++){
        typename boost::graph_traits<T_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<T_t>::vertex_descriptor t_node;
        unsigned int covered_count = 0;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            if(std::includes(bag(*tIt, T).begin(),
                             bag(*tIt, T).end(),
                             it->begin(), it->end())){
                t_node = *tIt;
                if(boost::degree(*tIt, T) <= 1 && *it == bag(*tIt, T)){
                    covered_count++;
                }
            }
        }
        if(covered_count == 0){
            typename boost::graph_traits<T_t>::vertex_descriptor new_t_node = boost::add_vertex(T);
            bag(new_t_node, T) = *it;
            boost::add_edge(t_node, new_t_node, T);
        }
        else{
            for(unsigned int j = 0; j < covered_count-1; j++){
                for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
                    if(*it == bag(*tIt, T)){
                        boost::clear_vertex(*tIt, T);
                        boost::remove_vertex(*tIt, T);
                    }
                }
            }
        }
    }

    //Remove bags of size 1.
    bool exists_singleton = true;
    while(exists_singleton){
        exists_singleton = false;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            if(bag(*tIt, T).size() == 1){
                boost::clear_vertex(*tIt, T);
                boost::remove_vertex(*tIt, T);
                exists_singleton = true;
            }
        }
    }


/* Make T cubic. */

    //Trivial cubic tree.
    if(boost::num_vertices(T) == 1){
        return;
    }

    bool is_cubic = false;
    while(!is_cubic){
        is_cubic = true;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            int degree = boost::degree(*tIt, T);
            if(degree == 3 || degree == 1){
                continue;
            }

            std::vector<typename boost::graph_traits<T_t>::vertex_descriptor > N(degree);

            typename boost::graph_traits<T_t>::adjacency_iterator  nIt, nEnd;
            unsigned int c = 0;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*tIt, T); nIt != nEnd; nIt++){
                N[c++] = *nIt;
            }

            if(degree == 2){
                boost::add_edge(N[0], N[1], T);
                boost::clear_vertex(*tIt, T);
                boost::remove_vertex(*tIt, T);
                is_cubic = false;
                break;
            }

            //degree > 3.
            typename boost::graph_traits<T_t>::vertex_descriptor new_t_node = boost::add_vertex(T);
            bag(new_t_node, T) = bag(*tIt, T);
            boost::add_edge(*tIt, new_t_node, T);
            boost::remove_edge(*tIt, N[0], T);
            boost::remove_edge(*tIt, N[1], T);
            boost::add_edge(new_t_node, N[0], T);
            boost::add_edge(new_t_node, N[1], T);
            is_cubic = false;
            break;
        }
    }

    //Remove "interior" bags.
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        if(boost::degree(*tIt, T) != 1){
            bag(*tIt, T).clear();
        }
    }
}

} //namespace treedec

#endif //TD_BRANCH_DECOMPOSITION

// vim:ts=8:sw=4:et
