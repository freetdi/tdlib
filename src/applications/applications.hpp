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
 */

#ifndef TREEDEC_APPLICATIONS_HPP
#define TREEDEC_APPLICATIONS_HPP

#include <map>
#include <set>
#include <vector>
#include <algorithm>

#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>

#include "nice_decomposition.hpp"
#include "simple_graph_algos.hpp"
#include "misc.hpp"
#include "graph.hpp"
#include "validation.hpp"

/*
#include "clique.hpp"
#include "independent_set.hpp"
#include "vertex_cover.hpp"
#include "dominating_set.hpp"
#include "coloring.hpp"
*/

namespace treedec{

static void powerset(std::set<unsigned int> &X, std::vector<std::set<unsigned int> > &subs){
    std::vector<unsigned int> sub;
    for(unsigned int i = 0; i <=X.size(); i++){
        subsets(X, X.size(), i, 0, sub, subs);
    }
}

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


unsigned min(unsigned a, unsigned b){
    return (a < b)? a : b;
}

template <typename I_t, typename B_t>
void make_set(I_t I, B_t &B){
    auto x = (*I).first;
    auto y = (*I).second;

    for(; x != y; ++x){
        B.insert(*x);
    }
}

/*
template <typename T_t>
bool test_encoding(){
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
            std::cerr << "error in encoding!" << std::endl;
            return false;
        }
    }

    return true;
}

*/

template <typename T_t>
class Intermediate_Results {
public: // types
    typedef unsigned Encoded_t; //TODO: should be choosable
    typedef std::set<unsigned> Decoded_t;
    typedef int Value_t; //TODO: for now

    typedef typename boost::graph_traits<T_t>::vertex_descriptor vd_t;

public: // construct
    Intermediate_Results(T_t &t) : _t(t) {
        _results.resize(boost::num_vertices(t));
    }

public: //interface
    void add(vd_t node, Encoded_t key, Value_t value){
        _results[node][key] = value;
    }

    Value_t get(vd_t node, Encoded_t key){
        return _results[node][key];
    }

    bool exists(vd_t node, Encoded_t key){
        return _results[node].find(key) != _results[node].end();
    }

    //bit-encoding of the subset relative to power set
    //TODO: what if |power_set| > available wordlength?
    //TODO: other possibility: pos of subset in enumeration of a power set. 
    //TODO: whatever works, must be constant amout of work
    unsigned encode(vd_t node, Decoded_t key){
        typename Decoded_t::iterator subIt, powIt;

        if(key.size() == 0){
            return 0;
        }

        powIt = bag(node, _t).begin();
        subIt = key.begin();

        unsigned exp = 0;
        unsigned number = 0;

        for(; subIt != key.end();){
            if(*subIt == *powIt){
                number += 1 << exp;
                subIt++;
            }
            powIt++;
            exp += 1;
        }

        return number;
    }

    void decode(vd_t node, Encoded_t key, Decoded_t &result){
        for(unsigned i = 0; i < min(CHAR_BIT*sizeof(key), bag(node, _t).size()); i++){
            if(key & 1){
                typename Decoded_t::iterator it = bag(node, _t).begin();
                std::advance(it, i);
                result.insert(*it);
            }
            key >>= 1;
        }
    }

    T_t &_t;
    typename std::vector<std::map<Encoded_t, Value_t> > _results;

}; //Intermediate_Results


template <typename T_t>
void top_down_computation(T_t &T,
                    typename boost::graph_traits<T_t>::vertex_descriptor cur,
                    typename std::vector<typename std::map<typename treedec_traits<T_t>::bag_type, int> > &results,
                    unsigned int val, typename treedec_traits<T_t>::bag_type &S,
                    typename treedec_traits<T_t>::bag_type &S_comp,
                    typename treedec_traits<T_t>::bag_type subset,
                    unsigned int take_flag)
{
    treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

    if(node_type == treedec::nice::LEAF){
        for(typename std::map<typename treedec_traits<T_t>::bag_type, int>::iterator it =
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

        for(typename std::map<typename treedec_traits<T_t>::bag_type, int>::iterator it =
                   results[cur].begin(); it != results[cur].end(); it++)
        {
            if((take_flag == 1 && it->first == subset)
            || (take_flag == 2 && it->second == (int)val 
               && std::includes(it->first.begin(), it->first.end(),
                                subset.begin(), subset.end()))
            || (take_flag == 0 && it->second == (int)val))
            {
                typename treedec_traits<T_t>::bag_type intersection;
                std::set_intersection(it->first.begin(), it->first.end(),
                                      S_comp.begin(), S_comp.end(),
                                      std::inserter(intersection, intersection.begin()));

                if(intersection.size() == 0){
                    std::set_difference(bag(cur, T).begin(), bag(cur, T).end(),
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
        for(typename std::map<typename treedec_traits<T_t>::bag_type, int>::iterator it =
                          results[cur].begin(); it != results[cur].end(); it++)
        {
            if((take_flag == 1 && it->first == subset)
            || (take_flag == 2 && it->second == (int)val
              && std::includes(it->first.begin(), it->first.end(),
                               subset.begin(), subset.end()))
            || (take_flag == 0 && it->second == (int)val))
            {
                typename treedec_traits<T_t>::bag_type intersection;
                std::set_intersection(it->first.begin(), it->first.end(),
                                      S_comp.begin(), S_comp.end(),
                                      std::inserter(intersection, intersection.begin()));

                if(intersection.size() == 0){
                    typename treedec_traits<T_t>::bag_type must_take = it->first;
                    S.insert(it->first.begin(), it->first.end());

                    std::set_difference(bag(cur, T).begin(), bag(cur, T).end(),
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

template <typename T_t>
void top_down_computation2(T_t &T,
                    typename boost::graph_traits<T_t>::vertex_descriptor cur,
                    treedec::app::detail::Intermediate_Results<T_t> &iRes,
                    unsigned int val, typename treedec_traits<T_t>::bag_type &S,
                    typename treedec_traits<T_t>::bag_type &S_comp,
                    typename treedec_traits<T_t>::bag_type subset,
                    unsigned int take_flag)
{
    treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

    if(node_type == treedec::nice::LEAF){
        for(typename std::map<unsigned, int>::iterator it =
                iRes._results[cur].begin(); it != iRes._results[cur].end(); it++)
        {
            if(it->second == (int)val){
                unsigned encoded = it->first;
                typename treedec_traits<T_t>::bag_type decoded_set;
                iRes.decode(cur, encoded, decoded_set);

                S.insert(decoded_set.begin(), decoded_set.end());
                return;
            }
        }
    }
    else if(node_type == treedec::nice::INTRODUCE || node_type == treedec::nice::FORGET){

        //TODO: necessary? subset could be the encoded one here
        unsigned subset_encoded = iRes.encode(cur, subset);

        if(take_flag == 1){
            val = iRes.get(cur, subset_encoded);
        }

        for(typename std::map<unsigned, int>::iterator it =
                   iRes._results[cur].begin(); it != iRes._results[cur].end(); it++)
        {
            unsigned encoded = it->first;
            typename treedec_traits<T_t>::bag_type decoded_set;
            iRes.decode(cur, encoded, decoded_set);

            //TODO: cleanup and comment
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
            top_down_computation2(T, child, iRes, val, S, S_comp, subset, 2);
        }
        else{
            subset.erase(treedec::nice::get_introduced_vertex(cur, T));
            top_down_computation2(T, child, iRes, val, S, S_comp, subset, 2);
        }
    }
    else if(node_type == treedec::nice::JOIN){
        for(typename std::map<unsigned, int>::iterator it =
                          iRes._results[cur].begin(); it != iRes._results[cur].end(); it++)
        {

            unsigned encoded = it->first;
            typename treedec_traits<T_t>::bag_type decoded_set;
            iRes.decode(cur, encoded, decoded_set);

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
                        unsigned encoded = iRes.encode(*nIt, must_take);
                        top_down_computation2(T, *nIt, iRes, iRes.get(*nIt, encoded), S, S_comp, must_take, 1);
                    }
                    return;
                }
            }
        }
    }
}

} //namespace detail


} //namespace app

} //namespace treedec

#endif //TREEDEC_APPLICATIONS_HPP

// vim:ts=8:sw=4:et
