// Lukas Larisch, 2014 - 2017
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
 * - void min_vertex_cover_with_treedecomposition(G_t&, T_t&, typename treedec_traits<T_t>::bag_type &result)
 *
 * IMPORT NOTE: ensure that the input treedecomposition is directed by
 *              using treedec::make_rooted(undir_t, dir_t)
 *
 */

#ifndef TREEDEC_VERTEX_COVER2_HPP
#define TREEDEC_VERTEX_COVER2_HPP

#include "applications.hpp"
#include <cmath>
#include <climits>

namespace treedec{

namespace app{

namespace detail{

template <typename T_t>
class Intermediate_Results {

}

} //namespace detail

namespace detail{

//bit-encoding of the subset relative to power_set
//TODO: what if |power_set| > available wordlength?
//TODO: other possibility: pos of subset in enumeration of a power set. 
//TODO: whatever works, must be constant amout of work
template <typename T_t>
unsigned encode_set(typename treedec_traits<T_t>::bag_type &power_set, typename treedec_traits<T_t>::bag_type &subset){
    typename treedec_traits<T_t>::bag_type::iterator subIt, powIt;

/*
    std::cout << "encoding" << std::endl;

    for(auto s = subset.begin(); s != subset.end(); s++){
        std::cout << *s << " ";
    } 
    std::cout << std::endl;
*/

    if(subset.size() == 0){
        return 0;
    }

    powIt = power_set.begin();
    subIt = subset.begin();

    unsigned exp = 0;
    unsigned number = 0;

    for(; subIt != subset.end();){
        if(*subIt == *powIt){
            number += pow(2, exp);
            subIt++;
        }
        powIt++;
        exp += 1;
    }
/*
    std::cout << "number: " << number << std::endl;
*/

    return number;
}

template <typename T_t>
void decode_set(typename treedec_traits<T_t>::bag_type &power_set, typename treedec_traits<T_t>::bag_type &subset, unsigned number){    
    for(unsigned i = 0; i < CHAR_BIT*sizeof(number); i++){
        if(number & 1){
            typename treedec_traits<T_t>::bag_type::iterator it = power_set.begin();
            std::advance(it, i);
            subset.insert(*it);
        }
        number >>= 1;
    }

/*
    std::cout << "decoding" << std::endl;

    for(auto s = subset.begin(); s != subset.end(); s++){
        std::cout << *s << " ";
    } 
    std::cout << std::endl;
*/
}

template <typename T_t>
void test_encoding()
{
    std::cout << "test encoding" << std::endl;

    typename treedec_traits<T_t>::bag_type S;
    S.insert(17);
    S.insert(1);
    S.insert(2);
    S.insert(3);
    S.insert(4);

    std::vector<typename treedec_traits<T_t>::bag_type> subs;
    treedec::powerset(S, subs);

    for(unsigned i = 0; i < subs.size(); i++){
        typename treedec_traits<T_t>::bag_type result;

        unsigned number = treedec::app::detail::encode_set<T_t>(S, subs[i]);
        treedec::app::detail::decode_set<T_t>(S, result, number);

        if(subs[i] != result){
            std::cout << "error in encoding!" << std::endl;
        }
    }

    std::cout << "end test encoding" << std::endl;
}

} //namespace detail

namespace detail{

template <typename T_t>
void top_down_computation2(T_t &T,
                    typename boost::graph_traits<T_t>::vertex_descriptor cur,
                    typename std::vector<std::map<unsigned, int> > &results,
                    unsigned int val, typename treedec_traits<T_t>::bag_type &S,
                    typename treedec_traits<T_t>::bag_type &S_comp,
                    typename treedec_traits<T_t>::bag_type subset,
                    unsigned int take_flag)
{
    treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

    if(node_type == treedec::nice::LEAF){
        for(typename std::map<unsigned, int>::iterator it =
                results[cur].begin(); it != results[cur].end(); it++)
        {
            if(it->second == (int)val){
                unsigned encoded = it->first;
                typename treedec_traits<T_t>::bag_type decoded_set;
                treedec::app::detail::decode_set<T_t>(bag(cur, T), decoded_set, encoded);

                S.insert(decoded_set.begin(), decoded_set.end());
                return;
            }
        }
    }
    else if(node_type == treedec::nice::INTRODUCE || node_type == treedec::nice::FORGET){

        //TODO: necessary? subset could be the encoded one here
        unsigned subset_encoded = treedec::app::detail::encode_set<T_t>(bag(cur, T), subset);

        if(take_flag == 1){
            val = results[cur][subset_encoded];
        }

        for(typename std::map<unsigned, int>::iterator it =
                   results[cur].begin(); it != results[cur].end(); it++)
        {
            unsigned encoded = it->first;
            typename treedec_traits<T_t>::bag_type decoded_set;
            treedec::app::detail::decode_set<T_t>(bag(cur, T), decoded_set, encoded);

            if((take_flag == 1 && decoded_set == subset)
              || (take_flag == 2 && it->second == (int)val 
              && std::includes(decoded_set.begin(), decoded_set.end(),
                                subset.begin(), subset.end()))
              || (take_flag == 0 && it->second == (int)val))
            {
                typename treedec_traits<T_t>::bag_type intersection;
                std::set_intersection(decoded_set.begin(), decoded_set.end(),
                                      S_comp.begin(), S_comp.end(),
                                      std::inserter(intersection, intersection.begin()));

                if(intersection.size() == 0){
                    std::set_difference(bag(cur, T).begin(), bag(cur, T).end(),
                                        decoded_set.begin(), decoded_set.end(),
                                        std::inserter(S_comp, S_comp.begin()));
                    S.insert(decoded_set.begin(), decoded_set.end());
                    subset = decoded_set;
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
            top_down_computation2(T, child, results, val, S, S_comp, subset, 2);
        }
        else{
            subset.erase(treedec::nice::get_introduced_vertex(cur, T));
            top_down_computation2(T, child, results, val, S, S_comp, subset, 2);
        }
    }
    else if(node_type == treedec::nice::JOIN){
        for(typename std::map<unsigned, int>::iterator it =
                          results[cur].begin(); it != results[cur].end(); it++)
        {

            unsigned encoded = it->first;
            typename treedec_traits<T_t>::bag_type decoded_set;
            treedec::app::detail::decode_set<T_t>(bag(cur, T), decoded_set, encoded);

            if((take_flag == 1 && decoded_set == subset)
            || (take_flag == 2 && it->second == (int)val
              && std::includes(decoded_set.begin(), decoded_set.end(),
                               subset.begin(), subset.end()))
            || (take_flag == 0 && it->second == (int)val))
            {
                typename treedec_traits<T_t>::bag_type intersection;
                std::set_intersection(decoded_set.begin(), decoded_set.end(),
                                      S_comp.begin(), S_comp.end(),
                                      std::inserter(intersection, intersection.begin()));

                if(intersection.size() == 0){
                    typename treedec_traits<T_t>::bag_type must_take = decoded_set;
                    S.insert(decoded_set.begin(), decoded_set.end());

                    std::set_difference(bag(cur, T).begin(), bag(cur, T).end(),
                                        decoded_set.begin(), decoded_set.end(),
                                        std::inserter(S_comp, S_comp.begin()));

                    typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(cur, T); nIt != nEnd; nIt++){
                        unsigned encoded = treedec::app::detail::encode_set<T_t>(bag(*nIt, T), must_take);
                        top_down_computation2(T, *nIt, results, results[*nIt][encoded], S, S_comp, must_take, 1);
                    }
                    return;
                }
            }
        }
    }
}

} //namespace detail


/* MIN VERTEX COVER */

namespace detail{

template <typename G_t>
bool is_vertex_cover2(G_t &G, 
    typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type &bag,
    const typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type &set)
{
    for(typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type::iterator sIt1
             = bag.begin(); sIt1 != bag.end(); sIt1++)
    {
        typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type::iterator sIt2 = sIt1;
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
unsigned int bottom_up_computation_vertex_cover2(G_t &G, T_t &T,
         std::vector<std::map<unsigned, int> > &results)
{
    std::stack<typename boost::graph_traits<T_t>::vertex_descriptor> S;
    treedec::nice::postorder_traversal(T, S);
    typename boost::graph_traits<T_t>::vertex_descriptor cur;

    while(!S.empty()){
        cur = S.top();
        S.pop();

        treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

        if(node_type == treedec::nice::LEAF){
            results[cur][0] = 0;
            results[cur][1] = 1;
        }
        else if(node_type == treedec::nice::INTRODUCE){
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                     *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor new_vertex =
                             treedec::nice::get_introduced_vertex(cur, T);

            //check for old set in new graph
            for(typename std::map<unsigned, int>::iterator it =
                         results[child].begin(); it != results[child].end(); it++)
            {
                unsigned old_encoded = it->first;
                typename treedec_traits<T_t>::bag_type decoded_set;
                treedec::app::detail::decode_set<T_t>(bag(child, T), decoded_set, old_encoded);

                unsigned new_encoded = treedec::app::detail::encode_set<T_t>(bag(cur, T), decoded_set);

                //TODO: just necessary to check for the new edges introduced by new_vertex
                if(is_vertex_cover2(G, bag(cur, T), decoded_set)){
                    results[cur][new_encoded] = results[child][old_encoded];
                }
                else{
                    results[cur][new_encoded] = -1;
                }
            }

            //check for new set in new graph
            for(typename std::map<unsigned, int>::iterator it =
                    results[child].begin(); it != results[child].end(); it++)
            {
                unsigned old_encoded = it->first;
                typename treedec_traits<T_t>::bag_type decoded_set;
                treedec::app::detail::decode_set<T_t>(bag(child, T), decoded_set, old_encoded);

                typename treedec_traits<T_t>::bag_type new_set = decoded_set;
                new_set.insert(new_vertex);

                unsigned new_encoded = treedec::app::detail::encode_set<T_t>(bag(cur, T), new_set);

                if(it->second != -1){
                    results[cur][new_encoded] = results[child][old_encoded] + 1;
                }
                //TODO: just necessary to check for the new edges introduced by new_vertex
                else{
                    if(is_vertex_cover2(G, bag(cur, T), new_set)){
                        results[cur][new_encoded] = new_set.size();
                    }
                    else{
                        results[cur][new_encoded] = -1;
                    }
                }
            }
        }
        else if(node_type == treedec::nice::FORGET){
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                     *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                          treedec::nice::get_forgotten_vertex(cur, T);

            //TODO: just enumerate it, not store it
            std::vector<typename treedec_traits<T_t>::bag_type> subs;
            treedec::powerset(bag(cur, T), subs);

            for(unsigned int i = 0; i < subs.size(); i++){
                typename treedec_traits<T_t>::bag_type new_set = subs[i];
                new_set.insert(forgotten_vertex);

                unsigned old_encoded1 = treedec::app::detail::encode_set<T_t>(bag(child, T), subs[i]);
                unsigned old_encoded2 = treedec::app::detail::encode_set<T_t>(bag(child, T), new_set);

                int val_with = results[child][old_encoded1];
                int val_without = results[child][old_encoded2];

                unsigned new_encoded = treedec::app::detail::encode_set<T_t>(bag(cur, T), subs[i]);

                if(val_with == -1){
                    results[cur][new_encoded] = val_without;
                }
                else if(val_without == -1){
                    results[cur][new_encoded] = val_with;
                }
                else{
                    if(val_with <= val_without){
                        results[cur][new_encoded] = val_with;
                    }
                    else{
                        results[cur][new_encoded] = val_without;
                    }
                }
            }
        }
        else if(node_type == treedec::nice::JOIN){
            typename boost::graph_traits<T_t>::vertex_descriptor child1 =
                                     *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<T_t>::vertex_descriptor child2 =
                                     *(++boost::adjacent_vertices(cur, T).first);

            //TODO: just enumerate it, not store it
            std::vector<typename treedec_traits<T_t>::bag_type> subs;
            treedec::powerset(bag(cur, T), subs);

            for(unsigned int i = 0; i < subs.size(); i++){
                //TODO: should be all the same, right?! Check this!
                unsigned encoded_left = treedec::app::detail::encode_set<T_t>(bag(child1, T), subs[i]);
                unsigned encoded_right = treedec::app::detail::encode_set<T_t>(bag(child2, T), subs[i]);
                unsigned encoded_cur = treedec::app::detail::encode_set<T_t>(bag(cur, T), subs[i]);

                if(results[child1][encoded_left] < 0 || results[child2][encoded_right] < 0){
                    results[cur][encoded_cur] = -1;
                }
                else{
                    results[cur][encoded_cur] = results[child1][encoded_left] + results[child2][encoded_right] - subs[i].size();
                }
            }
        }
    }

    typename boost::graph_traits<T_t>::vertex_descriptor root = find_root(T);
    return (unsigned int) results[root].begin()->second;
}

} //namespace detail (min_vertex_cover)


template <typename G_t, typename T_t>
unsigned int min_vertex_cover_with_treedecomposition2(G_t &G, T_t &T,
              typename treedec_traits<T_t>::bag_type &global_result)
{

    //treedec::app::detail::test_encoding<T_t>();

    std::vector<std::map<unsigned, int> > results(boost::num_vertices(T));

    unsigned int max = treedec::app::detail::bottom_up_computation_vertex_cover2(G, T, results);

    if(max > 0){
        typename treedec_traits<T_t>::bag_type a, b;
        typename boost::graph_traits<T_t>::vertex_descriptor root = find_root(T);
        treedec::app::detail::top_down_computation2(T, root, results, max, global_result, a, b, 0);
    }

    // assert(treedec::validation::is_valid_vertex_cover(G, result));

    return max;
}




} //namespace app

} //namespace treedec

#endif //TREEDEC_VERTEX_COVER_HPP

// vim:ts=8:sw=4:et
