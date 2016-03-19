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
 * - void min_coloring_with_treedecomposition(G_t &G, T_t &T, std::vector<unsigned int> &global_result)
 *
 */

#ifndef TD_APPLICATIONS
#define TD_APPLICATIONS

#include <map>
#include <set>
#include <vector>
#include <algorithm>

#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>

#include "TD_nice_decomposition.hpp"
#include "TD_simple_graph_algos.hpp"
#include "TD_misc.hpp"
#include "TD_noboost.hpp"

namespace treedec{

namespace app{

/* The top-down computation on tree decompositions is equal for some problems.
 *
 *   - 'cur' denotes the current node in T.
 *   - 'results' stores the computed table for all nodes of T.
 *   - 'val' denotes the value of a set that has to be choosen from tables of
 *       nodes in the subtree of 'T' with root 'cur'.
 *   - 'S' denotes the so far choosen set.
 *   - 'S_comp' denotes the complement of 'S' with respect to V(G).
 *   - according to (nice treedecomposition-) node types, the choice in the
 *       table of the child of 'cur' must be restricted:
 *         take_flag = 0 -> no restriction
 *         take_flag = 1 -> 'subset' must be choosen
 *         take_flag = 2 -> 'subset' must be a subset
 */
template <typename T_t>
void top_down_computation(T_t &T,
                    typename boost::graph_traits<T_t>::vertex_descriptor cur,
                    std::vector<std::map<std::set<unsigned int>, int> > &results,
                    unsigned int val, std::set<unsigned int> &S,
                    std::set<unsigned int> &S_comp, std::set<unsigned int> subset,
                    unsigned int take_flag)
{
    treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

    if(node_type == treedec::nice::LEAF){
        for(std::map<std::set<unsigned int>, int>::iterator it =
                results[cur].begin(); it != results[cur].end(); it++)
        {
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

        for(std::map<std::set<unsigned int>, int>::iterator it =
                   results[cur].begin(); it != results[cur].end(); it++)
        {
            if((take_flag == 1 && it->first == subset)
            || (take_flag == 2 && it->second == (int)val 
               && std::includes(it->first.begin(), it->first.end(),
                                subset.begin(), subset.end()))
            || (take_flag == 0 && it->second == (int)val))
            {
                std::set<unsigned int> intersection;
                std::set_intersection(it->first.begin(), it->first.end(),
                                      S_comp.begin(), S_comp.end(),
                                      std::inserter(intersection, intersection.begin()));

                if(intersection.size() == 0){
                    std::set_difference(T[cur].bag.begin(), T[cur].bag.end(),
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

        typename boost::graph_traits<T_t>::vertex_descriptor child =
                               *(boost::adjacent_vertices(cur, T).first);

        if(node_type == treedec::nice::FORGET){
            top_down_computation(T, child, results, val, S, S_comp, subset, 2);
        }
        else{
            subset.erase(treedec::nice::get_introduced_vertex_id(cur, T));
            top_down_computation(T, child, results, val, S, S_comp, subset, 2);
        }
    }
    else if(node_type == treedec::nice::JOIN){
        for(std::map<std::set<unsigned int>, int>::iterator it =
                          results[cur].begin(); it != results[cur].end(); it++)
        {
            if((take_flag == 1 && it->first == subset)
            || (take_flag == 2 && it->second == (int)val
              && std::includes(it->first.begin(), it->first.end(),
                               subset.begin(), subset.end()))
            || (take_flag == 0 && it->second == (int)val))
            {
                std::set<unsigned int> intersection;
                std::set_intersection(it->first.begin(), it->first.end(),
                                      S_comp.begin(), S_comp.end(),
                                      std::inserter(intersection, intersection.begin()));

                if(intersection.size() == 0){
                    std::set<unsigned int> must_take = it->first;
                    S.insert(it->first.begin(), it->first.end());

                    std::set_difference(T[cur].bag.begin(), T[cur].bag.end(),
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
/*
template <typename G_t, typename T_t>
void max_clique_with_treedecomposition(G_t &G, T_t &T,
                               std::vector<unsigned int> &global_result,
                               treedec::np::max_clique_base<G_t> &mclb)
{
    for(unsigned int i = 0; i < boost::num_vertices(T); i++){
        G_t H;
        treedec::induced_subgraph(H, G, T[i].bag);
        std::vector<unsigned int> result;
        mclb.max_clique(H, result);
        if(result.size() > global_result.size())
            global_result = result;
    }
}
*/


/* MAX INDEPENDENT SET */

template <typename G_t, typename T_t>
int bottom_up_computation_independent_set(G_t &G, T_t &T,
                  std::vector<std::map<std::set<unsigned int>, int> > &results)
{
    std::vector<bool> visited(boost::num_vertices(T), false);

    typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
    typename boost::graph_traits<T_t>::vertex_descriptor next = treedec::nice::next_node_postorder(root, T, visited);
    typename boost::graph_traits<T_t>::vertex_descriptor cur = next;

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    treedec::make_index_map(G, idxMap);

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
                treedec::powerset(T[cur].bag, subs);

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
                typename boost::graph_traits<G_t>::vertex_descriptor new_vertex =
                               treedec::nice::get_introduced_vertex<G_t, T_t>(cur, T, idxMap);

                results[cur] = results[child];

                for(std::map<std::set<unsigned int>, int>::iterator it =
                             results[child].begin(); it != results[child].end(); it++)
                {
                    bool extensible = true;

                    for(std::set<unsigned int>::iterator sIt =
                             it->first.begin(); sIt != it->first.end(); sIt++)
                    {
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
                typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                                    treedec::nice::get_forgotten_vertex<G_t, T_t>(cur, T, idxMap);

                std::vector<std::set<unsigned int> > subs;
                treedec::powerset(T[cur].bag, subs);

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
void max_independent_set_with_treedecomposition(G_t &G, T_t &T,
                          std::set<unsigned int> &global_result)
{
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
bool is_vertex_cover(G_t &G, std::set<unsigned int> &bag, std::set<unsigned int> &set,
             std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &idxMap)
{
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
int bottom_up_computation_vertex_cover(G_t &G, T_t &T,
         std::vector<std::map<std::set<unsigned int>, int> > &results)
{
    std::vector<bool> visited(boost::num_vertices(T), false);

    typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
    typename boost::graph_traits<T_t>::vertex_descriptor next = treedec::nice::next_node_postorder(root, T, visited);
    typename boost::graph_traits<T_t>::vertex_descriptor cur = next;

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    treedec::make_index_map(G, idxMap);

    while(cur != root){
        cur = next;

        treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

        next = treedec::nice::next_node_postorder(root, T, visited);

        if(node_type == treedec::nice::LEAF){
            results[cur][std::set<unsigned int>()] = 0;
            results[cur][T[cur].bag] = 1;
        }
        else{
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                     *(boost::adjacent_vertices(cur, T).first);

            if(node_type == treedec::nice::JOIN){
                typename boost::graph_traits<T_t>::vertex_descriptor child2 =
                                     *(++boost::adjacent_vertices(cur, T).first);

                std::vector<std::set<unsigned int> > subs;
                treedec::powerset(T[cur].bag, subs);

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
                typename boost::graph_traits<G_t>::vertex_descriptor new_vertex =
                             treedec::nice::get_introduced_vertex<G_t, T_t>(cur, T, idxMap);

                for(std::map<std::set<unsigned int>, int>::iterator it =
                             results[child].begin(); it != results[child].end(); it++)
                {
                    if(is_vertex_cover(G, T[cur].bag, it->first, idxMap)){
                        results[cur][it->first] = results[child][it->first];
                    }
                    else{
                        results[cur][it->first] = -1;
                    }
                }

                for(std::map<std::set<unsigned int>, int>::iterator it =
                        results[child].begin(); it != results[child].end(); it++)
                {
                    std::set<unsigned int> tmp = it->first;
                    tmp.insert(G[new_vertex].id);

                    if(it->second != -1){
                        results[cur][tmp] = results[child][it->first] + 1;
                    }
                    else{
                        if(is_vertex_cover(G, T[cur].bag, tmp, idxMap)){
                            results[cur][tmp] = tmp.size();
                        }
                        else{
                            results[cur][tmp] = -1;
                        }
                    }
                }
            }

            else if(node_type == treedec::nice::FORGET){
                typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                              treedec::nice::get_forgotten_vertex<G_t, T_t>(cur, T, idxMap);

                std::vector<std::set<unsigned int> > subs;
                treedec::powerset(T[cur].bag, subs);

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
bool is_dominating_set(G_t &G, std::set<unsigned int> &bag, const std::set<unsigned int> &set,
                       std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &idxMap)
{
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


template <typename G_t>
bool is_dominated(G_t &G, typename boost::graph_traits<G_t>::vertex_descriptor new_vertex,
                  const std::set<unsigned int> &set,
                  std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &idxMap)
{
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(new_vertex, G); nIt != nEnd; nIt++){
        if(set.find(G[*nIt].id) != set.end()){
            return true;
        }
    }
    return false;
}

template <typename G_t, typename T_t>
int bottom_up_computation_dominating_set(G_t &G, T_t &T,
                   std::vector<std::map<std::set<unsigned int>, int> > &results)
{
    std::vector<bool> visited(boost::num_vertices(T), false);

    typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
    typename boost::graph_traits<T_t>::vertex_descriptor next = treedec::nice::next_node_postorder(root, T, visited);
    typename boost::graph_traits<T_t>::vertex_descriptor cur = next;

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    treedec::make_index_map(G, idxMap);

    while(cur != root){
        cur = next;

        treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

        next = treedec::nice::next_node_postorder(root, T, visited);

        if(node_type == treedec::nice::LEAF){
            results[cur][std::set<unsigned int>()] = -1;
            results[cur][T[cur].bag] = 1;
        }
        else{
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                        *(boost::adjacent_vertices(cur, T).first);

            if(node_type == treedec::nice::JOIN){
                typename boost::graph_traits<T_t>::vertex_descriptor child2 =
                                        *(++boost::adjacent_vertices(cur, T).first);

                std::vector<std::set<unsigned int> > subs;
                treedec::powerset(T[cur].bag, subs);

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
                std::cout << "INTRODUCE NODE" << std::endl;
                std::cout << "bag: " << std::endl;
                print_results(T[cur].bag);
                typename boost::graph_traits<G_t>::vertex_descriptor new_vertex =
                                    treedec::nice::get_introduced_vertex<G_t, T_t>(cur, T, idxMap);

                for(std::map<std::set<unsigned int>, int>::iterator it =
                                results[child].begin(); it != results[child].end(); it++)
                {
                    std::cout << "old set: ";
                    print_results(it->first);
                    if(is_dominated(G, new_vertex, it->first, idxMap)){
                        std::cout << "is DS: true" << std::endl;
                        results[cur][it->first] = results[child][it->first];
                    }
                    else{
                        std::cout << "is DS: false" << std::endl;
                        results[cur][it->first] = -1;
                    }
                }
                std::cout << std::endl;

                for(std::map<std::set<unsigned int>, int>::iterator it =
                            results[child].begin(); it != results[child].end(); it++)
                {
                    std::set<unsigned int> tmp = it->first;
                    tmp.insert(G[new_vertex].id);

                    if(it->second != -1){
                        results[cur][tmp] = results[child][it->first] + 1;
                    }
                    else{
                        if(is_dominating_set(G, T[cur].bag, tmp, idxMap)){
                            results[cur][tmp] = tmp.size();
                        }
                        else{
                            results[cur][tmp] = -1;
                        }
                    }
                }
            }

            else if(node_type == treedec::nice::FORGET){
                typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                             treedec::nice::get_forgotten_vertex<G_t, T_t>(cur, T, idxMap);

                std::vector<std::set<unsigned int> > subs;
                treedec::powerset(T[cur].bag, subs);

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
        std::cout << "bag: " << std::endl;
        print_results(T[cur].bag);
        std::cout << "table; " << std::endl;
        print_results(results[cur]);
        std::cout << std::endl << std::endl;
    }

    return results[root].begin()->second;
}

unsigned int pow(unsigned int b, unsigned int n){
    unsigned int res = b;
    for(unsigned int i = 0; i < n-1; i++){
        res = res * b;
    }
    return res;
}

/*
template <typename G_t, typename T_t>
int bottom_up_computation_dominating_set(G_t &G, T_t &T, std::vector<std::vector<std::vector<unsigned int> > > &results, std::vector<std::vector<int> > &values){
    std::vector<bool> visited(boost::num_vertices(T), false);

    typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
    typename boost::graph_traits<T_t>::vertex_descriptor next = treedec::nice::next_node_postorder(root, T, visited);
    typename boost::graph_traits<T_t>::vertex_descriptor cur = next;

    std::vector<std::vector<std::vector<std::set<unsigned int> > > > partitions(boost::num_vertices(T)); //values = size() der ersten Menge f"ur SURE

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    treedec::make_index_map(G, idxMap);

    while(cur != root){
        cur = next;

        treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

        next = treedec::nice::next_node_postorder(root, T, visited);

        std::cout << "--------------------------------------" << std::endl;
        std::cout << "bag: " << std::endl;
        print_results(T[cur].bag);
        std::cout << std::endl << std::endl;

        if(node_type == treedec::nice::LEAF){
            //There are two valid colorings for one vertex.
            std::vector<std::set<unsigned int> > P(3);
            std::set<unsigned int> n, y;
            y.insert(*(T[cur].bag.begin()));
            P[NOT] = n;
            P[SURE] = y;
            P[MAYBE] = n;
            partitions[cur].push_back(P);
            P[SURE] = n;
            P[MAYBE] = y;
            partitions[cur].push_back(P);

            print_results_DS(partitions[cur]);
        }
        else{
            typename boost::graph_traits<T_t>::vertex_descriptor child = *(boost::adjacent_vertices(cur, T).first);

            if(node_type == treedec::nice::JOIN){
                //S.619 unten
            }

            if(node_type == treedec::nice::INTRODUCE){
                typename boost::graph_traits<G_t>::vertex_descriptor new_vertex =
                                    treedec::nice::get_introduced_vertex<G_t, T_t>(cur, T, idxMap);

                std::set<unsigned int> N;
                typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(idxMap[new_vertex], G); nIt != nEnd; nIt++){
                    N.insert(G[*nIt].id);
                }

                for(unsigned int i = 0; i < partitions[child].size(); i++){
                    //Modify old colorings.
                    std::vector<std::set<unsigned int> > old_partition = partitions[child][i];
                    std::vector<std::set<unsigned int> > modified_partition(3);
                    modified_partition[SURE] = old_partition[SURE];
                    modified_partition[MAYBE] = old_partition[MAYBE];
                    std::set_difference(old_partition[NOT].begin(), old_partition[NOT].end(),
                                        N.begin(), N.end(),
                                        std::inserter(modified_partition[NOT], modified_partition[NOT].begin()));
                    std::set_intersection(old_partition[NOT].begin(), old_partition[NOT].end(),
                                          N.begin(), N.end(),
                                          std::inserter(modified_partition[MAYBE], modified_partition[MAYBE].begin()));
                    //S.619 oben

                    //todo: coloring -> index in results[child] Bijektiv (gewichtete Summe "uber Eintr"age)

                    std::vector<std::set<unsigned int> > new_partition1 = partitions[child][i];
                    new_partition1[NOT].insert(G[new_vertex].id);
                    results[cur].push_back(new_partition1);
                    //values[cur].push_back(?);

                    std::vector<std::set<unsigned int> > new_partition2 = partitions[child][i];
                    new_partition2[SURE].insert(G[new_vertex].id);
                    results[cur].push_back(new_partition2);
                    //values[cur].push_back(?);

                    std::vector<std::set<unsigned int> > new_partition3 = partitions[child][i];
                    new_partition3[MAYBE].insert(G[new_vertex].id);
                    results[cur].push_back(new_partition3);
                    values[cur].push_back(values[child][i]);
                }
            }

            else if(node_type == treedec::nice::FORGET){
                typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                                      treedec::nice::get_forgotten_vertex<G_t, T_t>(cur, T, idxMap);
                //S.618 oben
                //f"ur alle F"arbungen c des FORGET-Beutel: Wert(c) entspricht min{Wert(c v forgotten=0), Wert(c v forgotten=1)}
            }
        }
    }

    //return results[root].begin()->second;
}
*/

template <typename G_t, typename T_t>
void min_dominating_set_with_treedecomposition(G_t &G, T_t &T, std::set<unsigned int> &global_result){
    std::vector<std::vector<std::vector<unsigned int> > > results(boost::num_vertices(T));
    std::vector<std::vector<int> > values(boost::num_vertices(T));

    //int max = 
    bottom_up_computation_dominating_set(G, T, results, values);
    //std::cout << "max: " << max << std::endl;
/*
    if(max > 0){
        std::set<unsigned int> a, b;
        typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
        top_down_computation(T, root, results, max, global_result, a, b, 0);
    }
*/
}


/* MIN COLORING */

template <typename G_t>
bool is_valid_coloring(G_t &G, std::set<unsigned int> &M, std::vector<int> &coloring){
    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        unsigned int sid = G[boost::source(*eIt, G)].id;
        unsigned int did = G[boost::target(*eIt, G)].id;
        if(M.find(sid) != M.end() && M.find(did) != M.end()){
            assert(coloring.size() > sid && coloring.size() > did);
            if(coloring[sid] == coloring[did]){
                return false;
            }
        }
    }
    return true;
}

template <typename G_t>
bool is_valid_coloring(G_t &G, typename boost::graph_traits<G_t>::vertex_descriptor v, std::vector<int> &coloring){
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    unsigned int vid = G[v].id;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; nIt++){
        unsigned int nid = G[*nIt].id;
        if(coloring[vid] == coloring[nid]){
            return false;
        }
    }
    return true;
}

template <typename G_t>
void all_k_colorings(unsigned int n, unsigned int k, unsigned int max,
        std::set<unsigned int> &M, std::vector<std::vector<int> > &colorings, G_t &G)
{
    std::vector<int> coloring(n, -1);
    std::set<unsigned int>::iterator iM = M.begin();
    while(iM != M.end()){
        coloring[*(iM++)]++;
    }

    iM = M.begin();
    unsigned int c = 0;

    if(treedec::app::is_valid_coloring(G, M, coloring)){
        colorings[c++] = coloring;
    }

    while(iM != M.end() && c < colorings.size()){
        if(coloring[*iM] < k-1){
            if(iM == M.end()){ break; }
            coloring[*iM]++;

            if(treedec::app::is_valid_coloring(G, M, coloring)){
                colorings[c++] = coloring;
            }
        }
        else{
            while(coloring[*iM] == k-1 && iM != M.end()){
                coloring[*iM] = 0;
                iM++;
            }
            if(iM == M.end()){ break; }

            coloring[*iM]++;

            if(treedec::app::is_valid_coloring(G, M, coloring)){
                colorings[c++] = coloring;
            }

            iM = M.begin();
        }
    }

    colorings.resize(c);
}

void _colorings_intersection(std::vector<std::vector<int> > &results_left,
                             std::vector<std::vector<int> > &results_right,
                             std::set<unsigned int> &bag,
                             std::vector<std::vector<int> > &intersection)
{
    for(unsigned int i = 0; i < results_left.size(); i++){
        for(unsigned int j = 0; j < results_right.size(); j++){
            bool is_compatible = true;
            for(std::set<unsigned int>::iterator sIt = bag.begin(); sIt != bag.end(); sIt++){
                if(results_left[i][*sIt] != results_right[j][*sIt]){
                    is_compatible = false;
                    break;
                }
            }
            if(is_compatible){
                std::vector<int> intersection_ = results_left[i];
                for(unsigned int l = 0; l < intersection_.size(); l++){
                    if(intersection_[l] == -1){ intersection_[l] = results_right[j][l]; }
                }
                intersection.push_back(intersection_);
            }
        }
    }
}

template <typename G_t, typename T_t>
bool bottom_up_computation_min_coloring(G_t &G, T_t &T, unsigned int k,
                        std::vector<std::vector<std::vector<int> > > &results)
{
    std::vector<bool> visited(boost::num_vertices(T), false);

    typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
    typename boost::graph_traits<T_t>::vertex_descriptor next = treedec::nice::next_node_postorder(root, T, visited);
    typename boost::graph_traits<T_t>::vertex_descriptor cur = next;

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    treedec::make_index_map(G, idxMap);
#ifdef DEBUG
    std::cout << "k=" << k << std::endl;
#endif
    while(cur != root){
        cur = next;

        treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

        next = treedec::nice::next_node_postorder(root, T, visited);

#ifdef DEBUG
        std::cout << "--------------------------------------" << std::endl;
        std::cout << "bag: " << std::endl;
        print_results(T[cur].bag);
        std::cout << std::endl << std::endl;
#endif
        if(node_type == treedec::nice::LEAF){
            unsigned int leaf = *(T[cur].bag.begin());

#ifdef DEBUG
            std::cout << "LEAF" << std::endl;
            std::cout << "leaf: " << leaf << std::endl;
#endif

            unsigned int entries_count = pow(k, T[cur].bag.size());
            std::vector<std::vector<int> > leaf_results(entries_count);
            treedec::app::all_k_colorings(boost::num_vertices(G), k, entries_count, T[cur].bag, leaf_results, G);

#ifdef DEBUG
            std::cout << "leaf results" << std::endl; print_results(leaf_results); std::cout << std::endl << std::endl;
#endif
            results[cur] = leaf_results;
        }
        else{
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                             *(boost::adjacent_vertices(cur, T).first);

            if(node_type == treedec::nice::JOIN){
                typename boost::graph_traits<T_t>::vertex_descriptor child2 =
                                             *(++boost::adjacent_vertices(cur, T).first);

                std::vector<std::vector<int> > join_results;
                treedec::app::_colorings_intersection(results[child], results[child2], T[cur].bag, join_results);
#ifdef DEBUG
                std::cout << "JOIN" << std::endl;
                std::cout << "join results" << std::endl; print_results(join_results); std::cout << std::endl << std::endl;
#endif
                if(join_results.size() == 0){
                    return false;
                }

                results[cur] = join_results;
            }

            if(node_type == treedec::nice::INTRODUCE){
                typename boost::graph_traits<G_t>::vertex_descriptor new_vertex =
                                    treedec::nice::get_introduced_vertex<G_t, T_t>(cur, T, idxMap);
                unsigned int vid = G[new_vertex].id;

#ifdef DEBUG
                std::cout << "INTRODUCE" << std::endl;
                std::cout << "introduced vertex: " << vid << std::endl;
#endif

                //For all colorings C of the child bag and all possible colors c, test if
                //C extended by c is a valid coloring C'. If this is the case, add C' as
                //a solution to introduce_results.
                std::vector<std::vector<int> > introduce_results;
                for(unsigned int i = 0; i < results[child].size(); i++){
                    std::vector<int> ext_coloring = results[child][i];
                    for(unsigned int j = 0; j < k; j++){
                        ext_coloring[vid] = j;
                        if(treedec::app::is_valid_coloring(G, new_vertex, ext_coloring)){
                            introduce_results.push_back(ext_coloring);
                        }
                    }
                }
#ifdef DEBUG
                std::cout << "introduce results" << std::endl; print_results(introduce_results); std::cout << std::endl << std::endl;
#endif
                if(introduce_results.size() == 0){
                    return false;
                }

                results[cur] = introduce_results;
            }
            else if(node_type == treedec::nice::FORGET){
                typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                                           treedec::nice::get_forgotten_vertex<G_t, T_t>(cur, T, idxMap);
                unsigned int id = G[forgotten_vertex].id;

#ifdef DEBUG
                std::cout << "FORGET" << std::endl;
                std::cout << "forgotten vertex: " << id << std::endl;
#endif

                //Set the 'id'-th entry of an solution of the child bag to -1.
                //Two solutions which are modified in this way have to be
                //consecutive in results[child].
                if(results[child].size() == 1){
                    results[cur].push_back(results[child][0]);
                }

                for(unsigned int i = 1; i < results[child].size(); i++){
                    std::vector<int> last_result = results[child][i-1];
                    std::vector<int> current_result = results[child][i];
                    last_result[id] = -1;
                    current_result[id] = -1;
                    if(current_result != last_result){
                        results[cur].push_back(last_result);
                    }
                    if(i == results[child].size()-1){
                        results[cur].push_back(current_result);
                    }
                }
#ifdef DEBUG
                std::cout << "forget results" << std::endl; print_results(results[cur]); std::cout << std::endl << std::endl;
#endif
            }
        }
    }

    std::cout << "finished! k=" << k << std::endl;

    return true;
}

template <typename G_t, typename T_t>
void top_down_computation_min_coloring(G_t &G, T_t &T,
                                  typename boost::graph_traits<T_t>::vertex_descriptor cur,
                                  std::vector<std::vector<std::vector<int> > > &results,
                                  std::vector<int> &global_result)
{
    treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    treedec::make_index_map(G, idxMap);

#ifdef DEBUG
        std::cout << "--------------------------------------" << std::endl;
        std::cout << "bag: " << std::endl;
        print_results(T[cur].bag);
        std::cout << std::endl << std::endl;

    std::cout << "coloring:" << std::endl;
    print_results(global_result);
#endif

    if(node_type == treedec::nice::LEAF){
#ifdef DEBUG
std::cout << "LEAF" << std::endl;
#endif
    }
    else if(node_type == treedec::nice::INTRODUCE){
#ifdef DEBUG
std::cout << "INTRODUCE" << std::endl;
#endif
        typename boost::graph_traits<T_t>::vertex_descriptor child =
                                 *(boost::adjacent_vertices(cur, T).first);

        top_down_computation_min_coloring(G, T, child, results, global_result);
    }
    else if(node_type == treedec::nice::FORGET){
        typename boost::graph_traits<T_t>::vertex_descriptor child =
                                 *(boost::adjacent_vertices(cur, T).first);

        typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                                 treedec::nice::get_forgotten_vertex<G_t, T_t>(cur, T, idxMap);

#ifdef DEBUG
std::cout << "FORGET" << std::endl;
        std::cout << "forgotten_vertex: " << G[forgotten_vertex].id << std::endl;

        std::cout << "results_forget (bottom-up)" << std::endl;
        print_results(results[child]);
#endif

        //
        for(unsigned int i = 0; i < results[child].size(); i++){
            bool valid = true;
            for(unsigned int j = 0; j < results[child][i].size(); j++){
                if(results[child][i][j] >= 0){
                    if(global_result[j] >= 0 && results[child][i][j] != global_result[j]){
                        valid = false;
                        break;
                    }
                }
            }
            if(valid){
                unsigned int id = G[forgotten_vertex].id;
                global_result[id] = results[child][i][id];
                break;
            }
        }

        if(global_result[G[forgotten_vertex].id] == -1){  std::cerr << "error: forget" << std::endl; exit(-1); }

        top_down_computation_min_coloring(G, T, child, results, global_result);
    }
    else if(node_type == treedec::nice::JOIN){
#ifdef DEBUG
std::cout << "JOIN" << std::endl;
#endif
        typename boost::graph_traits<T_t>::vertex_descriptor lchild =
                                   *(boost::adjacent_vertices(cur, T).first);
        typename boost::graph_traits<T_t>::vertex_descriptor rchild =
                                   *(++boost::adjacent_vertices(cur, T).first);

        top_down_computation_min_coloring(G, T, lchild, results, global_result);
        top_down_computation_min_coloring(G, T, rchild, results, global_result);
    }
}


template <typename G_t, typename T_t>
void min_coloring_with_treedecomposition(G_t &G, T_t &T, std::vector<int> &global_result)
{
    std::vector<std::vector<std::vector<int> > > results(boost::num_vertices(T));

    if(boost::num_vertices(T) == 1){
        return;
    }
    if(boost::num_edges(G) == 0){
        return;
    }

    unsigned int k = 2;

    while(!treedec::app::bottom_up_computation_min_coloring(G, T, k, results)){
        k++;
        results.clear();
        results.resize(boost::num_vertices(T));
    }

    global_result.assign(boost::num_vertices(G), -1);
    typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
    treedec::app::top_down_computation_min_coloring(G, T, root, results, global_result);

    std::cout << "k: " << k << std::endl;
}

} //namespace app

} //namespace treedec

#endif //TD_APPLICATIONS

// vim:ts=8:sw=4:et
