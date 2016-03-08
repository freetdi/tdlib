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

/* Offers functionality to solve hard problems with help of treedecompositions.
 *
 * Provides following functions (namespace treedec::app):
 *
 * - void max_clique_with_treedecomposition(G_t &G, T_t &T, std::vector<unsigned int> &global_result,
 *                                                          treedec::np::max_clique_base<G_t> &mclb)
 * - void max_independent_set_with_treedecomposition(G_t &G, T_t &T, std::set<unsigned int> &global_result)
 * - void min_vertex_cover_with_treedecomposition(G_t &G, T_t &T, std::set<unsigned int> &global_result)
 * - void min_dominating_set_with_treedecomposition(G_t &G, T_t &T, std::set<unsigned int> &global_result)
 *
 */

#ifndef TD_APPLICATIONS
#define TD_APPLICATIONS

#include <map>
#include <set>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>

#include "TD_nice_decomposition.hpp"
#include "TD_simple_graph_algos.hpp"
#include "TD_hard_problems.hpp"
#include "TD_misc.hpp"
#include "TD_noboost.hpp"

namespace treedec{

namespace app{

/* The top-down computation on tree decompositions is equal for different problems.
 *
 *   - 'cur' denotes the current node in T
 *   - 'results' stores the computed table for all nodes of T
 *   - 'val' denotes the value of a set that has to be choosen from tables of
 *       nodes in the subtree of T with root cur
 *   - 'S' denotes the so far choosen set
 *   - 'S_comp' denotes the complement of S with respect to V(G)
 *   - according to (nice treedecomposition-) node types, the choice in the
 *       table of the child of 'cur' must be restricted:
 *         take_flag = 0 -> no restriction
 *         take_flag = 1 -> 'subset' must be choosen
 *         take_flag = 2 -> 'subset' must be a subset
 */
template <typename T_t>
void top_down_computation(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor cur,
                                          std::vector<std::map<std::set<unsigned int>, int> > &results, unsigned int val,
                                          std::set<unsigned int> &S, std::set<unsigned int> &S_comp,
                                          std::set<unsigned int> subset, unsigned int take_flag){

    treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

    if(node_type == treedec::nice::LEAF){
        for(std::map<std::set<unsigned int>, int>::iterator it = results[cur].begin(); it != results[cur].end(); it++){
            if(it->second == (int)val){
                S.insert(it->first.begin(), it->first.end());
                return;
            }
        }
    }
    else if(node_type == treedec::nice::INTRODUCE || node_type == treedec::nice::FORGET){
        if(take_flag == 1){
            val = results[cur][subset];
        }

        for(std::map<std::set<unsigned int>, int>::iterator it = results[cur].begin(); it != results[cur].end(); it++){
            if((take_flag == 1 && it->first == subset)
            || (take_flag == 2 && it->second == (int)val && std::includes(it->first.begin(), it->first.end(), subset.begin(), subset.end()))
            || (take_flag == 0 && it->second == (int)val)){
                std::set<unsigned int> intersection;
                std::set_intersection(it->first.begin(), it->first.end(),
                                      S_comp.begin(), S_comp.end(),
                                      std::inserter(intersection, intersection.begin()));

                if(intersection.size() == 0){
                    std::set_difference(noboost::bag(T, cur).begin(),
                                        noboost::bag(T, cur).end(),
                                        it->first.begin(), it->first.end(),
                                        std::inserter(S_comp, S_comp.begin()));
                    S.insert(it->first.begin(), it->first.end());
                    subset = it->first;
                    break;
                }
            }
        }

        if(node_type == treedec::nice::INTRODUCE){
            if(S.find(treedec::nice::get_introduced_vertex_id(cur, T)) != S.end()){
                val = val - 1;
            }
        }

        typename boost::graph_traits<T_t>::vertex_descriptor child = *(boost::adjacent_vertices(cur, T).first);

        if(node_type == treedec::nice::FORGET){
            top_down_computation(T, child, results, val, S, S_comp, subset, 2);
        }
        else{
            subset.erase(treedec::nice::get_introduced_vertex_id(cur, T));
            top_down_computation(T, child, results, val, S, S_comp, subset, 2);
        }
    }
    else if(node_type == treedec::nice::JOIN){
        for(std::map<std::set<unsigned int>, int>::iterator it = results[cur].begin(); it != results[cur].end(); it++){
            if((take_flag == 1 && it->first == subset)
            || (take_flag == 2 && it->second == (int)val && std::includes(it->first.begin(), it->first.end(), subset.begin(), subset.end()))
            || (take_flag == 0 && it->second == (int)val)){
                std::set<unsigned int> intersection;
                std::set_intersection(it->first.begin(), it->first.end(), S_comp.begin(), S_comp.end(),
                                                                           std::inserter(intersection, intersection.begin()));

                if(intersection.size() == 0){
                    std::set<unsigned int> must_take = it->first;
                    S.insert(it->first.begin(), it->first.end());

                    std::set_difference(noboost::bag(T, cur).begin(),
                                        noboost::bag(T, cur).end(),
                                        it->first.begin(), it->first.end(),
                                        std::inserter(S_comp, S_comp.begin()));

                    typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(cur, T); nIt != nEnd; nIt++){
                        top_down_computation(T, *nIt, results, results[*nIt][must_take], S, S_comp, must_take, 1);
                    }
                    return;
                }
            }
        }
    }
}


/* MAX CLIQUE */

template <typename G_t, typename T_t>
void max_clique_with_treedecomposition(G_t &G, T_t &T, std::vector<unsigned int> &global_result, treedec::np::max_clique_base<G_t> &mclb){
    for(unsigned int i = 0; i < boost::num_vertices(T); i++){
        G_t H;
        induced_subgraph(H, G, noboost::bag(T, i));
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
            results[cur][noboost::bag(T, cur)] = 1;
        }
        else{
            typename boost::graph_traits<T_t>::vertex_descriptor child = *(boost::adjacent_vertices(cur, T).first);

            if(node_type == treedec::nice::JOIN){
                typename boost::graph_traits<T_t>::vertex_descriptor child2 = *(++boost::adjacent_vertices(cur, T).first);

                std::vector<std::set<unsigned int> > subs;
                powerset(noboost::bag(T, cur), subs);

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
                powerset(noboost::bag(T, cur), subs);

                for(unsigned int i = 0; i < subs.size(); i++){
                    std::set<unsigned int> tmp = subs[i];
                    tmp.insert(G[forgotten_vertex].id);

                    int val_with = results[child][subs[i]];
                    int val_without = results[child][tmp];

                    if(val_with >= val_without){
                        results[cur][subs[i]] = val_with;
                    }
                    else{
                        results[cur][subs[i]] = val_without;
                    }
                }
            }
        }
    }

    return results[root].begin()->second;
}



template <typename G_t, typename T_t>
void max_independent_set_with_treedecomposition(G_t &G, T_t &T, std::set<unsigned int> &global_result){
    std::vector<std::map<std::set<unsigned int>, int> > results(boost::num_vertices(T));

    int max = bottom_up_computation_independent_set(G, T, results);

    if(max > 0){
        std::set<unsigned int> a, b;
        typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
        top_down_computation(T, root, results, max, global_result, a, b, 0);
    }
}


/* MIN VERTEX COVER */

template <typename G_t>
bool is_vertex_cover(G_t &G, std::set<unsigned int> &bag, std::set<unsigned int> &set, std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &idxMap){
    for(std::set<unsigned int>::iterator sIt1 = bag.begin(); sIt1 != bag.end(); sIt1++){
        std::set<unsigned int>::iterator sIt2 = sIt1;
        sIt2++;
        for(; sIt2 != bag.end(); sIt2++){
            if(boost::edge(idxMap[*sIt1], idxMap[*sIt2], G).second){
                if(set.find(*sIt1) == set.end() && set.find(*sIt2) == set.end()){
                    return false;
                }
            }
        }
    }
    return true;
}

template <typename G_t>
bool is_vertex_cover(G_t &G, std::set<unsigned int> &bag, const std::set<unsigned int> &set, std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &idxMap){
    for(std::set<unsigned int>::iterator sIt1 = bag.begin(); sIt1 != bag.end(); sIt1++){
        std::set<unsigned int>::iterator sIt2 = sIt1;
        sIt2++;
        for(; sIt2 != bag.end(); sIt2++){
            if(boost::edge(idxMap[*sIt1], idxMap[*sIt2], G).second){
                if(set.find(*sIt1) == set.end() && set.find(*sIt2) == set.end()){
                    return false;
                }
            }
        }
    }
    return true;
}


template <typename G_t, typename T_t>
int bottom_up_computation_vertex_cover(G_t &G, T_t &T, std::vector<std::map<std::set<unsigned int>, int> > &results){
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
            results[cur][noboost::bag(T, cur)] = 1;
        }
        else{
            typename boost::graph_traits<T_t>::vertex_descriptor child = *(boost::adjacent_vertices(cur, T).first);

            if(node_type == treedec::nice::JOIN){
                typename boost::graph_traits<T_t>::vertex_descriptor child2 = *(++boost::adjacent_vertices(cur, T).first);

                std::vector<std::set<unsigned int> > subs;
                powerset(noboost::bag(T, cur), subs);

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

                for(std::map<std::set<unsigned int>, int>::iterator it = results[child].begin(); it != results[child].end(); it++){
                    if(is_vertex_cover(G, noboost::bag(T, cur), it->first, idxMap)){
                        results[cur][it->first] = results[child][it->first];
                    }
                    else{
                        results[cur][it->first] = -1;
                    }
                }

                for(std::map<std::set<unsigned int>, int>::iterator it = results[child].begin(); it != results[child].end(); it++){
                    std::set<unsigned int> tmp = it->first;
                    tmp.insert(G[new_vertex].id);

                    if(it->second != -1){
                        results[cur][tmp] = results[child][it->first] + 1;
                    }
                    else{
                        if(is_vertex_cover(G, noboost::bag(T, cur), tmp, idxMap)){
                            results[cur][tmp] = tmp.size();
                        }
                        else{
                            results[cur][tmp] = -1;
                        }
                    }
                }
            }

            else if(node_type == treedec::nice::FORGET){
                typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex = treedec::nice::get_forgotten_vertex<G_t, T_t>(cur, T, idxMap);

                std::vector<std::set<unsigned int> > subs;
                powerset(noboost::bag(T, cur), subs);

                for(unsigned int i = 0; i < subs.size(); i++){
                    std::set<unsigned int> tmp = subs[i];
                    tmp.insert(G[forgotten_vertex].id);

                    int val_with = results[child][subs[i]];
                    int val_without = results[child][tmp];

                    if(val_with == -1){
                        results[cur][subs[i]] = val_without;
                    }
                    else if(val_without == -1){
                        results[cur][subs[i]] = val_with;
                    }
                    else{
                        if(val_with <= val_without){
                            results[cur][subs[i]] = val_with;
                        }
                        else{
                            results[cur][subs[i]] = val_without;
                        }
                    }
                }
            }
        }
    }

    return results[root].begin()->second;
}


template <typename G_t, typename T_t>
void min_vertex_cover_with_treedecomposition(G_t &G, T_t &T, std::set<unsigned int> &global_result){
    std::vector<std::map<std::set<unsigned int>, int> > results(boost::num_vertices(T)); 

    int max = bottom_up_computation_vertex_cover(G, T, results);

    if(max > 0){
        std::set<unsigned int> a, b;
        typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
        top_down_computation(T, root, results, max, global_result, a, b, 0);
    }
}


/* MIN DOMINATING SET */

template <typename G_t>
bool is_dominating_set(G_t &G, std::set<unsigned int> &bag, const std::set<unsigned int> &set, std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &idxMap){
    for(std::set<unsigned int>::iterator sIt = bag.begin(); sIt != bag.end(); sIt++){
        if(set.find(*sIt) == set.end()){
            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            bool hit = false;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(idxMap[*sIt], G); nIt != nEnd; nIt++){
                if(set.find(G[*nIt].id) != set.end()){
                    hit = true;
                    break;
                }
            }
            if(!hit){
                return false;
            }
        }
    }
    return true;
}

template <typename G_t, typename T_t>
int bottom_up_computation_dominating_set(G_t &G, T_t &T, std::vector<std::map<std::set<unsigned int>, int> > &results){
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
            results[cur][std::set<unsigned int>()] = -1;
            results[cur][noboost::bag(T, cur)] = 1;
        }
        else{
            typename boost::graph_traits<T_t>::vertex_descriptor child = *(boost::adjacent_vertices(cur, T).first);

            if(node_type == treedec::nice::JOIN){
                typename boost::graph_traits<T_t>::vertex_descriptor child2 = *(++boost::adjacent_vertices(cur, T).first);

                std::vector<std::set<unsigned int> > subs;
                powerset(noboost::bag(T, cur), subs);

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

                for(std::map<std::set<unsigned int>, int>::iterator it = results[child].begin(); it != results[child].end(); it++){
                    if(is_dominating_set(G, noboost::bag(T, cur), it->first, idxMap)){
                        results[cur][it->first] = results[child][it->first];
                    }
                    else{
                        results[cur][it->first] = -1;
                    }
                }

                for(std::map<std::set<unsigned int>, int>::iterator it = results[child].begin(); it != results[child].end(); it++){
                    std::set<unsigned int> tmp = it->first;
                    tmp.insert(G[new_vertex].id);

                    if(it->second != -1){
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
                powerset(noboost::bag(T, cur), subs);

                for(unsigned int i = 0; i < subs.size(); i++){
                    std::set<unsigned int> tmp = subs[i];
                    tmp.insert(G[forgotten_vertex].id);

                    int val_with = results[child][subs[i]];
                    int val_without = results[child][tmp];

                    if(val_with == -1){
                        results[cur][subs[i]] = val_without;
                    }
                    else if(val_without == -1){
                        results[cur][subs[i]] = val_with;
                    }
                    else{
                        if(val_with <= val_without){
                            results[cur][subs[i]] = val_with;
                        }
                        else{
                            results[cur][subs[i]] = val_without;
                        }
                    }
                }
            }
        }
    }

    return results[root].begin()->second;
}

template <typename G_t, typename T_t>
void min_dominating_set_with_treedecomposition(G_t &G, T_t &T, std::set<unsigned int> &global_result){
    std::vector<std::map<std::set<unsigned int>, int> > results(boost::num_vertices(T)); 

    int max = bottom_up_computation_dominating_set(G, T, results);

    if(max > 0){
        std::set<unsigned int> a, b;
        typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
        top_down_computation(T, root, results, max, global_result, a, b, 0);
    }
}


/* MIN COLORING */

template <typename G_t, typename T_t>
void min_coloring_with_treedecomposition(G_t &G, T_t &T, std::vector<std::set<unsigned int> > &global_result){
}

} //namespace app

} //namespace treedec

#endif

// vim:ts=8:sw=4:et
