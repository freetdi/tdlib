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
 * - void max_clique_with_treedecomposition(G_t&, T_t&, typename noboost::treedec_traits<T_t>::bag_type &result)
 * - void max_independent_set_with_treedecomposition(G_t&, T_t&, typename noboost::treedec_traits<T_t>::bag_type &result)
 * - void min_vertex_cover_with_treedecomposition(G_t&, T_t&, typename noboost::treedec_traits<T_t>::bag_type &result)
 * - void min_dominating_set_with_treedecomposition(G_t&, T_t&, typename noboost::treedec_traits<T_t>::bag_type &result)
 * - void min_coloring_with_treedecomposition(G_t&, T_t&,
 *                   std::vector<typename noboost::treedec_traits<T_t>::bag_type> &result)
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

namespace detail{

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
                    typename std::vector<typename std::map<typename noboost::treedec_traits<T_t>::bag_type, int> > &results,
                    unsigned int val, typename noboost::treedec_traits<T_t>::bag_type &S,
                    typename noboost::treedec_traits<T_t>::bag_type &S_comp,
                    typename noboost::treedec_traits<T_t>::bag_type subset,
                    unsigned int take_flag)
{
    treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

    if(node_type == treedec::nice::LEAF){
        for(typename std::map<typename noboost::treedec_traits<T_t>::bag_type, int>::iterator it =
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

        for(typename std::map<typename noboost::treedec_traits<T_t>::bag_type, int>::iterator it =
                   results[cur].begin(); it != results[cur].end(); it++)
        {
            if((take_flag == 1 && it->first == subset)
            || (take_flag == 2 && it->second == (int)val 
               && std::includes(it->first.begin(), it->first.end(),
                                subset.begin(), subset.end()))
            || (take_flag == 0 && it->second == (int)val))
            {
                typename noboost::treedec_traits<T_t>::bag_type intersection;
                std::set_intersection(it->first.begin(), it->first.end(),
                                      S_comp.begin(), S_comp.end(),
                                      std::inserter(intersection, intersection.begin()));

                if(intersection.size() == 0){
                    std::set_difference(noboost::bag(cur, T).begin(), noboost::bag(cur, T).end(),
                                        it->first.begin(), it->first.end(),
                                        std::inserter(S_comp, S_comp.begin()));
                    S.insert(it->first.begin(), it->first.end());
                    subset = it->first;
                    break;
                }
            }
        }

        if(node_type == treedec::nice::INTRODUCE){
            if(S.find(treedec::nice::get_introduced_vertex(cur, T)) != S.end()){
                val = val - 1;
            }
        }

        typename boost::graph_traits<T_t>::vertex_descriptor child =
                               *(boost::adjacent_vertices(cur, T).first);

        if(node_type == treedec::nice::FORGET){
            top_down_computation(T, child, results, val, S, S_comp, subset, 2);
        }
        else{
            subset.erase(treedec::nice::get_introduced_vertex(cur, T));
            top_down_computation(T, child, results, val, S, S_comp, subset, 2);
        }
    }
    else if(node_type == treedec::nice::JOIN){
        for(typename std::map<typename noboost::treedec_traits<T_t>::bag_type, int>::iterator it =
                          results[cur].begin(); it != results[cur].end(); it++)
        {
            if((take_flag == 1 && it->first == subset)
            || (take_flag == 2 && it->second == (int)val
              && std::includes(it->first.begin(), it->first.end(),
                               subset.begin(), subset.end()))
            || (take_flag == 0 && it->second == (int)val))
            {
                typename noboost::treedec_traits<T_t>::bag_type intersection;
                std::set_intersection(it->first.begin(), it->first.end(),
                                      S_comp.begin(), S_comp.end(),
                                      std::inserter(intersection, intersection.begin()));

                if(intersection.size() == 0){
                    typename noboost::treedec_traits<T_t>::bag_type must_take = it->first;
                    S.insert(it->first.begin(), it->first.end());

                    std::set_difference(noboost::bag(cur, T).begin(), noboost::bag(cur, T).end(),
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

} //namespace detail


/* MAX CLIQUE */

namespace detail{

template <typename G_t, typename T_t>
bool is_clique(G_t &G, typename noboost::treedec_traits<T_t>::bag_type &C){
    typename noboost::treedec_traits<T_t>::bag_type::iterator sIt1, sIt2;
    for(sIt1 = C.begin(); sIt1 != C.end(); sIt1++){
        sIt2 = sIt1;
        sIt2++;
        for(; sIt2 != C.end(); sIt2++){
            if(!boost::edge(*sIt1, *sIt2, G).second){
                return false;
            }
        }
    }
    return true;
}

} //namespace detail (for max_clique)


//TODO: bottom-up computation would fasten the computation.
template <typename G_t, typename T_t>
unsigned int max_clique_with_treedecomposition(G_t &G, T_t &T,
                               typename noboost::treedec_traits<T_t>::bag_type &global_result)
{
    unsigned int max = 0;

    typename boost::graph_traits<T_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(T); vIt != vEnd; vIt++){
        if(noboost::bag(*vIt, T).size() < max){ continue; }

        std::vector<typename noboost::treedec_traits<T_t>::bag_type> subs;
        treedec::powerset(noboost::bag(*vIt, T), subs);

        for(int i = subs.size()-1; i > 0; i--){
            if(subs[i].size() <= max){ continue; }
            if(treedec::app::detail::is_clique<G_t, T_t>(G, subs[i])){
                max = subs[i].size();
                global_result = subs[i];
                break;
            }
        }
    }

    return max;
}


/* MAX INDEPENDENT SET */

namespace detail{

template <typename G_t, typename T_t>
unsigned int bottom_up_computation_independent_set(G_t &G, T_t &T,
       std::vector<std::map<typename noboost::treedec_traits<T_t>::bag_type, int> > &results)
{
    std::stack<typename boost::graph_traits<T_t>::vertex_descriptor> S;
    treedec::nice::postorder_traversal(T, S);
    typename boost::graph_traits<T_t>::vertex_descriptor cur;

    while(!S.empty()){
        cur = S.top();
        S.pop();

        treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

        if(node_type == treedec::nice::LEAF){
            //Store both possibilities (the empty set and the set containing one vertex).
            results[cur][typename noboost::treedec_traits<T_t>::bag_type()] = 0;
            results[cur][noboost::bag(cur, T)] = 1;
        }
        else if(node_type == treedec::nice::INTRODUCE){
            //For all results S of the child: Store S extended by the introduced vertex with value
            //old_value + 1, if S extended by the introduced vertex is an independent set, else -1.
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                      *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor new_vertex =
                                     treedec::nice::get_introduced_vertex(cur, T);

            results[cur] = results[child];

            for(typename std::map<typename noboost::treedec_traits<T_t>::bag_type, int>::iterator it =
                         results[child].begin(); it != results[child].end(); it++)
            {
                bool extensible = true;

                for(typename noboost::treedec_traits<T_t>::bag_type::iterator sIt =
                         it->first.begin(); sIt != it->first.end(); sIt++)
                {
                    if(boost::edge(*sIt, new_vertex, G).second){
                        extensible = false;
                        break;
                    }
                }

                typename noboost::treedec_traits<T_t>::bag_type tmp = it->first;
                tmp.insert(new_vertex);

                if(extensible){
                    results[cur][tmp] = results[child][it->first] + 1;
                }
                else{
                    results[cur][tmp] = -1;
                }
            }
        }
        else if(node_type == treedec::nice::FORGET){
            //Store the maximum of S and S extended by the forgotten vertex for each subset S
            //of the current bag.
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                      *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                                treedec::nice::get_forgotten_vertex(cur, T);

            std::vector<typename noboost::treedec_traits<T_t>::bag_type> subs;
            treedec::powerset(noboost::bag(cur, T), subs);

            for(unsigned int i = 0; i < subs.size(); i++){
                typename noboost::treedec_traits<T_t>::bag_type tmp = subs[i];
                tmp.insert(forgotten_vertex);

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

        if(node_type == treedec::nice::JOIN){
            //The maximum size of an independent set is the sum of the maximum size of an independent set
            //of the left child and the right child minus the size of the intersection with the current bag,
            //or -1 if the at least one of the child's results is not an independent set.
            typename boost::graph_traits<T_t>::vertex_descriptor child1 =
                                      *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<T_t>::vertex_descriptor child2 =
                                     *(++boost::adjacent_vertices(cur, T).first);

            std::vector<typename noboost::treedec_traits<T_t>::bag_type> subs;
            treedec::powerset(noboost::bag(cur, T), subs);

            for(unsigned int i = 0; i < subs.size(); i++){
                if(results[child1][subs[i]] < 0 || results[child2][subs[i]] < 0){
                    results[cur][subs[i]] = -1;
                }
                else{
                    results[cur][subs[i]] = results[child1][subs[i]] + results[child2][subs[i]] - subs[i].size();
                }
            }
        }
    }

    typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
    return (unsigned int) results[root].begin()->second;
}

} //namespace detail (for independent set)


template <typename G_t, typename T_t>
unsigned int max_independent_set_with_treedecomposition(G_t &G, T_t &T,
                      typename noboost::treedec_traits<T_t>::bag_type &global_result)
{
    std::vector<std::map<typename noboost::treedec_traits<T_t>::bag_type, int> > results(boost::num_vertices(T));

    unsigned int max = treedec::app::detail::bottom_up_computation_independent_set(G, T, results);

    if(max > 0){
        typename noboost::treedec_traits<T_t>::bag_type a, b;
        typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
        treedec::app::detail::top_down_computation(T, root, results, max, global_result, a, b, 0);
    }

    return max;
}


/* MIN VERTEX COVER */

namespace detail{

template <typename G_t>
bool is_vertex_cover(G_t &G, 
    typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &bag,
    const typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &set)
{
    for(typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type::iterator sIt1
             = bag.begin(); sIt1 != bag.end(); sIt1++)
    {
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type::iterator sIt2 = sIt1;
        sIt2++;
        for(; sIt2 != bag.end(); sIt2++){
            if(boost::edge(*sIt1, *sIt2, G).second){
                if(set.find(*sIt1) == set.end() && set.find(*sIt2) == set.end()){
                    return false;
                }
            }
        }
    }
    return true;
}


template <typename G_t, typename T_t>
unsigned int bottom_up_computation_vertex_cover(G_t &G, T_t &T,
         std::vector<std::map<typename noboost::treedec_traits<T_t>::bag_type, int> > &results)
{
    std::stack<typename boost::graph_traits<T_t>::vertex_descriptor> S;
    treedec::nice::postorder_traversal(T, S);
    typename boost::graph_traits<T_t>::vertex_descriptor cur;

    while(!S.empty()){
        cur = S.top();
        S.pop();

        treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

        if(node_type == treedec::nice::LEAF){
            results[cur][typename noboost::treedec_traits<T_t>::bag_type()] = 0;
            results[cur][noboost::bag(cur, T)] = 1;
        }
        else if(node_type == treedec::nice::INTRODUCE){
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                     *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor new_vertex =
                             treedec::nice::get_introduced_vertex(cur, T);

            for(typename std::map<typename noboost::treedec_traits<T_t>::bag_type, int>::iterator it =
                         results[child].begin(); it != results[child].end(); it++)
            {
                if(is_vertex_cover(G, noboost::bag(cur, T), it->first)){
                    results[cur][it->first] = results[child][it->first];
                }
                else{
                    results[cur][it->first] = -1;
                }
            }

            for(typename std::map<typename noboost::treedec_traits<T_t>::bag_type, int>::iterator it =
                    results[child].begin(); it != results[child].end(); it++)
            {
                typename noboost::treedec_traits<T_t>::bag_type tmp = it->first;
                tmp.insert(new_vertex);

                if(it->second != -1){
                    results[cur][tmp] = results[child][it->first] + 1;
                }
                else{
                    if(is_vertex_cover(G, noboost::bag(cur, T), tmp)){
                        results[cur][tmp] = tmp.size();
                    }
                    else{
                        results[cur][tmp] = -1;
                    }
                }
            }
        }
        else if(node_type == treedec::nice::FORGET){
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                     *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                          treedec::nice::get_forgotten_vertex(cur, T);

            std::vector<typename noboost::treedec_traits<T_t>::bag_type> subs;
            treedec::powerset(noboost::bag(cur, T), subs);

            for(unsigned int i = 0; i < subs.size(); i++){
                typename noboost::treedec_traits<T_t>::bag_type tmp = subs[i];
                tmp.insert(forgotten_vertex);
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
        else if(node_type == treedec::nice::JOIN){
            typename boost::graph_traits<T_t>::vertex_descriptor child1 =
                                     *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<T_t>::vertex_descriptor child2 =
                                     *(++boost::adjacent_vertices(cur, T).first);

            std::vector<typename noboost::treedec_traits<T_t>::bag_type> subs;
            treedec::powerset(noboost::bag(cur, T), subs);

            for(unsigned int i = 0; i < subs.size(); i++){
                if(results[child1][subs[i]] < 0 || results[child2][subs[i]] < 0){
                    results[cur][subs[i]] = -1;
                }
                else{
                    results[cur][subs[i]] = results[child1][subs[i]] + results[child2][subs[i]] - subs[i].size();
                }
            }
        }
    }

    typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
    return (unsigned int) results[root].begin()->second;
}

} //namespace detail (min_vertex_cover)


template <typename G_t, typename T_t>
unsigned int min_vertex_cover_with_treedecomposition(G_t &G, T_t &T,
              typename noboost::treedec_traits<T_t>::bag_type &global_result)
{
    std::vector<std::map<typename noboost::treedec_traits<T_t>::bag_type, int> > results(boost::num_vertices(T)); 

    unsigned int max = treedec::app::detail::bottom_up_computation_vertex_cover(G, T, results);

    if(max > 0){
        typename noboost::treedec_traits<T_t>::bag_type a, b;
        typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
        treedec::app::detail::top_down_computation(T, root, results, max, global_result, a, b, 0);
    }

    return max;
}


/* MIN DOMINATING SET */

/*
template <typename G_t>
bool is_dominating_set(G_t &G, std::set<unsigned int> &bag, const std::set<unsigned int> &set){
    for(std::set<unsigned int>::iterator sIt = bag.begin(); sIt != bag.end(); sIt++){
        if(set.find(*sIt) == set.end()){
            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            bool hit = false;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*sIt, G); nIt != nEnd; nIt++){
                if(set.find(*nIt) != set.end()){
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
                  const std::set<unsigned int> &set)
{
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(new_vertex, G); nIt != nEnd; nIt++){
        if(set.find(*nIt) != set.end()){
            return true;
        }
    }
    return false;
}
*/

/*
unsigned int pow(unsigned int b, unsigned int n){
    unsigned int res = b;
    for(unsigned int i = 0; i < n-1; i++){
        res = res * b;
    }
    return res;
}
*/

#if 0

template <typename G_t, typename T_t>
unsigned int bottom_up_computation_dominating_set(G_t &G, T_t &T, std::vector<std::vector<std::vector<unsigned int> > > &results, std::vector<std::vector<int> > &values){
    std::stack<typename boost::graph_traits<T_t>::vertex_descriptor> S;
    treedec::nice::postorder_traversal(T, S);
    typename boost::graph_traits<T_t>::vertex_descriptor cur;

    while(!S.empty()){
        cur = S.top();
        S.pop();

        treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

        if(node_type == treedec::nice::LEAF){
        }
        else if(node_type == treedec::nice::INTRODUCE){
        }
        else if(node_type == treedec::nice::FORGET){
        }
        else if(node_type == treedec::nice::JOIN){
        }
    }
}

template <typename G_t, typename T_t>
unsigned int min_dominating_set_with_treedecomposition(G_t &G, T_t &T,
                  typename noboost::treedec_traits<T_t>::bag_type &global_result)
{
    std::vector<std::vector<std::vector<int> > > results(boost::num_vertices(T));

    unsigned int max = -1;
    //unsigned int max = treedec::app::detail::bottom_up_computation_dominating_set(G, T, results, values);

/*
    if(max > 0){
        std::set<unsigned int> a, b;
        typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
        treedec::app::detail::top_down_computation_dominating_set(T, root, results, max, global_result, a, b, 0);
    }
*/

    return (unsigned int) max;
}

#endif //if 0

/* MIN COLORING */

namespace detail{

template <typename G_t>
bool is_valid_extended_coloring(G_t &G, typename boost::graph_traits<G_t>::vertex_descriptor v, std::vector<int> &coloring){
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    unsigned int vpos = noboost::get_pos(v, G);
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; nIt++){
        unsigned int npos = noboost::get_pos(*nIt, G);
        if(coloring[vpos] == coloring[npos]){
            return false;
        }
    }
    return true;
}

//Computes the "intersection" of two vectors of colorings with respect to 'bag'.
//Two colorings "intersect" with respect to 'bag', if they color the vertices
//in 'bag' with the same colors.
template <typename G_t, typename T_t>
void colorings_intersection(G_t &G, std::vector<std::vector<int> > &results_left,
                            std::vector<std::vector<int> > &results_right,
                            typename noboost::treedec_traits<T_t>::bag_type &bag,
                            std::vector<std::vector<int> > &intersection)
{
    for(unsigned int i = 0; i < results_left.size(); i++){
        for(unsigned int j = 0; j < results_right.size(); j++){
            bool is_compatible = true;
            for(typename noboost::treedec_traits<T_t>::bag_type::iterator sIt 
                    = bag.begin(); sIt != bag.end(); sIt++)
            {
                unsigned int pos = noboost::get_pos(*sIt, G);
                if(results_left[i][pos] != results_right[j][pos]){
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
    std::stack<typename boost::graph_traits<T_t>::vertex_descriptor> S;
    treedec::nice::postorder_traversal(T, S);
    typename boost::graph_traits<T_t>::vertex_descriptor cur;

    while(!S.empty()){
        cur = S.top();
        S.pop();

        treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

        if(node_type == treedec::nice::LEAF){
            //Store all k colorings.
            unsigned int leaf = *(noboost::bag(cur, T).begin());
            unsigned int pos = noboost::get_pos(leaf, G);

            std::vector<int> coloring(boost::num_vertices(G), -1);
            for(unsigned int i = 0; i < k; i++){
                coloring[pos] = i;
                results[cur].push_back(coloring);
           }
        }
        else if(node_type == treedec::nice::INTRODUCE){
            //Store all valid extensions of colorings of the child bag.
            //For all colorings C of the child bag and all possible colors c, test if
            //C extended by c is a valid coloring C'. If this is the case, store C' as
            //a solution.
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                         *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor new_vertex =
                                treedec::nice::get_introduced_vertex(cur, T);
            unsigned int pos = noboost::get_pos(new_vertex, G);

            for(unsigned int i = 0; i < results[child].size(); i++){
                std::vector<int> ext_coloring = results[child][i];
                for(unsigned int j = 0; j < k; j++){
                    ext_coloring[pos] = j;
                    if(treedec::app::detail::is_valid_extended_coloring(G, new_vertex, ext_coloring)){
                        results[cur].push_back(ext_coloring);
                    }
                }
            }
            //No extension possible with k colors.
            if(results[cur].size() == 0){
                return false;
            }
        }
        else if(node_type == treedec::nice::FORGET){
            //If two colorings C1 and C2 of the child bag are equal in positions != 'pos', store
            //the coloring C1 with C1[pos] = -1 (forget the coloring of the forgotten vertex).
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                             *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                                       treedec::nice::get_forgotten_vertex(cur, T);
            unsigned int pos = noboost::get_pos(forgotten_vertex, G);

            std::vector<bool> visited(results[child].size(), false);
            unsigned int colorings_size = boost::num_vertices(G);
            for(unsigned int i = 0; i < results[child].size(); i++){
                if(visited[i]){ continue; }
                visited[i] = true;
                std::vector<int> coloring = results[child][i];
                for(unsigned int j = i+1; j < results[child].size(); j++){
                    bool same = true;
                    for(unsigned int l = 0; l < colorings_size; l++){
                        if(l == pos){ continue; }
                        if(results[child][i][l] != results[child][j][l]){
                            same = false; break;
                        }
                    }
                    if(same){
                        visited[j] = true;
                    }
                }
                coloring[pos] = -1;
                results[cur].push_back(coloring);
            }
        }
        else if(node_type == treedec::nice::JOIN){
            //Store the "intersection" of childs. Two colorings "intersect", if they are equal
            //in the positions of the current bag.
            typename boost::graph_traits<T_t>::vertex_descriptor child1 =
                                         *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<T_t>::vertex_descriptor child2 =
                                         *(++boost::adjacent_vertices(cur, T).first);

            treedec::app::detail::colorings_intersection<G_t, T_t>(G, results[child1], results[child2], noboost::bag(cur, T), results[cur]);

            //No coloring possible with k colors.
            if(results[cur].size() == 0){
                return false;
            }
        }
    }

    return true;
}


template <typename G_t, typename T_t>
void top_down_computation_min_coloring(G_t &G, T_t &T,
                                  typename boost::graph_traits<T_t>::vertex_descriptor cur,
                                  std::vector<std::vector<std::vector<int> > > &results,
                                  std::vector<int> &global_result)
{
    treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

    if(node_type == treedec::nice::LEAF){
        //Nothing to do.
    }
    else if(node_type == treedec::nice::INTRODUCE){
        //Nothing to do.
        typename boost::graph_traits<T_t>::vertex_descriptor child =
                                 *(boost::adjacent_vertices(cur, T).first);

        top_down_computation_min_coloring(G, T, child, results, global_result);
    }
    else if(node_type == treedec::nice::FORGET){
        //There must be a coloring stored in the child's results that uses the the same
        //colors as the current global coloring and that has an entry != -1 in position
        //pos(forgotten_vertex). Search this coloring and update the global coloring
        //in this entry.
        typename boost::graph_traits<T_t>::vertex_descriptor child =
                                 *(boost::adjacent_vertices(cur, T).first);

        typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                                 treedec::nice::get_forgotten_vertex(cur, T);

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
                unsigned int pos = noboost::get_pos(forgotten_vertex, G);
                global_result[pos] = results[child][i][pos];
                break;
            }
        }

        top_down_computation_min_coloring(G, T, child, results, global_result);
    }
    else if(node_type == treedec::nice::JOIN){
        //Nothing to do.
        typename boost::graph_traits<T_t>::vertex_descriptor lchild =
                                   *(boost::adjacent_vertices(cur, T).first);
        typename boost::graph_traits<T_t>::vertex_descriptor rchild =
                                   *(++boost::adjacent_vertices(cur, T).first);

        top_down_computation_min_coloring(G, T, lchild, results, global_result);
        top_down_computation_min_coloring(G, T, rchild, results, global_result);
    }
}

} //namespace detail (min_coloring)


template <typename G_t, typename T_t>
unsigned int min_coloring_with_treedecomposition(G_t &G, T_t &T,
                          std::vector<typename noboost::treedec_traits<T_t>::bag_type> &global_result)
{
    std::vector<std::vector<std::vector<int> > > results(boost::num_vertices(T));


    if(boost::num_vertices(G) == 0){
        return 0;
    }
    if(boost::num_vertices(T) == 1){
        global_result.resize(1);
        global_result[0].insert(*(boost::vertices(G).first));
        return 1;
    }

    unsigned int k = 2;

    while(!treedec::app::detail::bottom_up_computation_min_coloring(G, T, k, results)){
        k++;
        results.clear();
        results.resize(boost::num_vertices(T));
    }

    std::vector<int> global_results_map(boost::num_vertices(G), -1);
    typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
    treedec::app::detail::top_down_computation_min_coloring(G, T, root, results, global_results_map);

    typename std::map<unsigned int, typename boost::graph_traits<G_t>::vertex_descriptor> inv_map;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int pos = noboost::get_pos(*vIt, G);
        inv_map[pos] = *vIt;
    }

    global_result.resize(k);
    for(unsigned int i = 0; i < global_results_map.size(); i++){
        unsigned int pos = global_results_map[i];
        global_result[pos].insert(inv_map[pos]);
    }

    return k;
}

} //namespace app

} //namespace treedec

#endif //TD_APPLICATIONS

// vim:ts=8:sw=4:et
