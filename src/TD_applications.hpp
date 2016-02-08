// Lukas Larisch, 2014 - 2016
//
// (c) 2014-2016 Goethe-Universität Frankfurt
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
// Offers functionality to solve hard problems with help of treedecompositions
//

#ifndef TD_APPLICATIONS
#define TD_APPLICATIONS

#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>

#include "TD_nice_decomposition.hpp"
#include "TD_simple_graph_algos.hpp"
#include "TD_hard_problems.hpp"
#include "TD_misc.hpp"

namespace treedec{

namespace app{

/* MAX CLIQUE */

template <typename G_t, typename T_t>
void max_clique_with_treedecomposition(G_t &G, T_t &T, std::vector<unsigned int> &global_result, treedec::np::max_clique_base<G_t> &mclb){
    for(unsigned int i = 0; i < boost::num_vertices(T); i++){
        G_t H;
        induced_subgraph(H, G, T[i].bag);
        std::vector<unsigned int> result;
        mclb.max_clique(H, result);
        if(result.size() > global_result.size())
            global_result = result;
    }
}

/* MAX INDEPENDENT SET */

template <typename G_t, typename T_t>
int bottom_up_computation_independent_set(G_t &G, T_t &T, std::vector<std::map<std::set<unsigned int>, int> > &results){
    std::vector<bool> visited(boost::num_vertices(T), false);

    typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
    typename boost::graph_traits<T_t>::vertex_descriptor next = treedec::nice::next_node_postorder(root, T, visited);
    typename boost::graph_traits<T_t>::vertex_descriptor cur = next;

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    make_index_map(G, idxMap);

    while(cur != root){
        cur = next;

        treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

        next = treedec::nice::next_node_postorder(root, T, visited);

        if(node_type == treedec::nice::LEAF){
            results[cur][std::set<unsigned int>()] = 0;
            results[cur][T[cur].bag] = 1;
        }
        else{
            typename boost::graph_traits<T_t>::vertex_descriptor child = *(boost::adjacent_vertices(cur, T).first);

            if(node_type == treedec::nice::JOIN){
                typename boost::graph_traits<T_t>::vertex_descriptor child2 = *(++boost::adjacent_vertices(cur, T).first);

                std::vector<std::set<unsigned int> > subs;
                powerset(T[cur].bag, subs);

                for(unsigned int i = 0; i < subs.size(); i++){
                    if(results[child][subs[i]] < 0 || results[child2][subs[i]] < 0){
                        results[cur][subs[i]] = -1;
                    }
                    else{
                        results[cur][subs[i]] = results[child][subs[i]] + results[child2][subs[i]] - subs[i].size();
                    }
                }
            }

            if(node_type == treedec::nice::INTRODUCE){
                typename boost::graph_traits<G_t>::vertex_descriptor new_vertex = treedec::nice::get_introduced_vertex<G_t, T_t>(cur, T, idxMap);

                results[cur] = results[child];

                for(std::map<std::set<unsigned int>, int>::iterator it = results[child].begin(); it != results[child].end(); it++){
                    bool extensible = true;

                    for(std::set<unsigned int>::iterator sIt = it->first.begin(); sIt != it->first.end(); sIt++){
                        if(boost::edge(idxMap[*sIt], new_vertex, G).second){
                            extensible = false;
                            break;
                        }
                    }

                    std::set<unsigned int> tmp = it->first;
                    tmp.insert(G[new_vertex].id);

                    if(extensible){
                        results[cur][tmp] = results[child][it->first] + 1;
                    }
                    else{
                        results[cur][tmp] = -1;
                    }
                }
            }

            else if(node_type == treedec::nice::FORGET){
                typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex = treedec::nice::get_forgotten_vertex<G_t, T_t>(cur, T, idxMap);

                std::vector<std::set<unsigned int> > subs;
                powerset(T[cur].bag, subs);

                for(unsigned int i = 0; i < subs.size(); i++){
                    std::set<unsigned int> tmp = subs[i];
                    tmp.insert(G[forgotten_vertex].id);

                    int val_with = results[child][subs[i]];
                    int val_without = results[child][tmp];

                    if(val_with >= val_without){
                        results[cur][subs[i]] = val_with;
                    }
                    else if(val_without > val_with){
                        results[cur][subs[i]] = val_without;
                    }
                }
            }
        }
    }

    return results[root].begin()->second;
}

template <typename T_t>
void top_down_computation_independent_set(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor cur, 
                                          std::vector<std::map<std::set<unsigned int>, int> > &results, unsigned int val, 
                                          std::set<unsigned int> &IS, std::set<unsigned int> &IS_comp, 
                                          std::set<unsigned int> have_to_take, unsigned int take_flag){

    treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

    if(node_type == treedec::nice::LEAF){
        for(std::map<std::set<unsigned int>, int>::iterator it = results[cur].begin(); it != results[cur].end(); it++){
            if(it->second == (int)val){
                IS.insert(it->first.begin(), it->first.end());
                return;
            }
        }
    }
    else if(node_type == treedec::nice::INTRODUCE || node_type == treedec::nice::FORGET){
        if(take_flag == 1){
            val = results[cur][have_to_take];
        }

        for(std::map<std::set<unsigned int>, int>::iterator it = results[cur].begin(); it != results[cur].end(); it++){
            if((take_flag == 1 && it->first == have_to_take)
            || (take_flag == 2 && it->second == (int)val && std::includes(it->first.begin(), it->first.end(), have_to_take.begin(), have_to_take.end()))
            || (take_flag == 0 && it->second == (int)val)){
                std::set<unsigned int> intersection;
                std::set_intersection(it->first.begin(), it->first.end(), IS_comp.begin(), IS_comp.end(), 
                                                                           std::inserter(intersection, intersection.begin()));

                if(intersection.size() == 0){
                    std::set_difference(T[cur].bag.begin(), T[cur].bag.end(), it->first.begin(), it->first.end(), 
                                                                           std::inserter(IS_comp, IS_comp.begin()));
                    IS.insert(it->first.begin(), it->first.end());
                    have_to_take = it->first;
                    break;
                }
            }
        }

        if(node_type == treedec::nice::INTRODUCE){
            if(IS.find(treedec::nice::get_introduced_vertex_id(cur, T)) != IS.end()){
                val = val - 1;
            }
        }

        typename boost::graph_traits<T_t>::vertex_descriptor child = *(boost::adjacent_vertices(cur, T).first);

        if(node_type == treedec::nice::FORGET){
            top_down_computation_independent_set(T, child, results, val, IS, IS_comp, have_to_take, 2);
        }
        else{
            have_to_take.erase(treedec::nice::get_introduced_vertex_id(cur, T));
            top_down_computation_independent_set(T, child, results, val, IS, IS_comp, have_to_take, 2);
        }
    }
    else if(node_type == treedec::nice::JOIN){
        for(std::map<std::set<unsigned int>, int>::iterator it = results[cur].begin(); it != results[cur].end(); it++){
            if((take_flag == 1 && it->first == have_to_take)
            || (take_flag == 2 && it->second == (int)val && std::includes(it->first.begin(), it->first.end(), have_to_take.begin(), have_to_take.end()))
            || (take_flag == 0 && it->second == (int)val)){
                std::set<unsigned int> intersection;
                std::set_intersection(it->first.begin(), it->first.end(), IS_comp.begin(), IS_comp.end(), 
                                                                           std::inserter(intersection, intersection.begin()));

                if(intersection.size() == 0){
                    std::set<unsigned int> must_take = it->first;
                    IS.insert(it->first.begin(), it->first.end());

                    std::set_difference(T[cur].bag.begin(), T[cur].bag.end(), it->first.begin(), it->first.end(), 
                                                                               std::inserter(IS_comp, IS_comp.begin()));

                    typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(cur, T); nIt != nEnd; nIt++){
                        top_down_computation_independent_set(T, *nIt, results, results[*nIt][must_take], IS, IS_comp, must_take, 1);
                    }
                    return;
                }
            }
        }

    }
}

template <typename G_t, typename T_t>
void max_independent_set_with_treedecomposition(G_t &G, T_t &T, std::set<unsigned int> &global_result){
    std::vector<std::map<std::set<unsigned int>, int> > results(boost::num_vertices(T)); 

    int max = bottom_up_computation_independent_set(G, T, results);

    if(max > 0){
        std::set<unsigned int> a, b;
        typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
        top_down_computation_independent_set(T, root, results, max, global_result, a, b, 0);
    }
}

/* MIN VERTEX COVER */

template <typename G_t, typename T_t>
void min_vertex_cover_with_treedecomposition(G_t &G, T_t &T, std::set<unsigned int> &global_result){
}

/* MIN DOMINATING SET */

template <typename G_t, typename T_t>
void min_dominating_set_with_treedecomposition(G_t &G, T_t &T, std::set<unsigned int> &global_result){
}

/* MIN COLORING */

template <typename G_t, typename T_t>
void min_coloring_with_treedecomposition(G_t &G, T_t &T, std::vector<std::set<unsigned int> > &global_result){
}

} //namespace app

} //namespace treedec

#endif


