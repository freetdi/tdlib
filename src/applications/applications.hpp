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


unsigned min(unsigned a, unsigned b){
    return (a < b)? a : b;
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

template <typename I_t>
class encoded_iterator{
public:
    encoded_iterator(unsigned number, I_t superIt, I_t superEnd)
     : _num(number), _supIt(superIt), _supEnd(superEnd){
        if(_num == 0){
             _supIt = _supEnd;
        }
        else{
            while(_num != 0){
                if(_num & 1){
                    _last = *_supIt;
                    _num >>= 1;
                    break;
                }
                _num >>= 1;
                _supIt++;
            }
        }      
    }

    bool operator==(const encoded_iterator& o) const{
        return _supIt == o._supIt;
    }

    bool operator!=(const encoded_iterator& o) const{
        return !operator==(o);
    }

    bool operator==(const I_t& o) const{
        return _supIt == o;
    }

    bool operator!=(const I_t& o) const{
        return !operator==(o);
    }

    template <typename R_t>
    bool operator==(const R_t& o) const{
        return _supIt == o;
    }

    template <typename R_t>
    bool operator!=(const R_t& o) const{
        return !operator==(o);
    }

    void operator++(){
        if(_num == 0){
             _supIt = _supEnd;
        }
        else{
            while(_num != 0){
                if(_num & 1){
                    _supIt++;
                    _last = *_supIt;
                    _num >>= 1;
                    break;
                }
                _num >>= 1;
                _supIt++;
            }
        }
    }

    unsigned operator*() const{
        return _last;
    }

private:
    unsigned _num;
    unsigned _last;

    I_t _supIt;
    I_t _supEnd;
};


template <typename R_t>
unsigned get_size(R_t sIt, R_t sEnd){
    unsigned n;
    while(sIt != sEnd){
        sIt++;
        n++;
    }
    return n;
}



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

#if 0
    bool exists(vd_t node, Encoded_t key){
        return _results[node].find(key) != _results[node].end();
    }
#endif

    //bit-encoding of the subset relative to power set
    //TODO: what if |power_set| > available wordlength?
    //TODO: other possibility: pos of subset in enumeration of a power set. 
    //TODO: whatever works, must be constant amout of work
    //encode from set
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

    //encode from iterator range
    template <typename R_t>
    unsigned encode(vd_t node, R_t sIt, R_t sEnd){
        typename Decoded_t::iterator subIt, powIt;

        if(sIt == sEnd){
            return 0;
        }

        powIt = bag(node, _t).begin();
        
        unsigned number = 0;
        unsigned summand = 1;

        for(; sIt != sEnd;){
            if(*sIt == *powIt){
                number += summand;
                sIt++;
            }
            powIt++;
            summand <<= 1;
        }

        return number;
    }

    //different supersets
    unsigned encode(vd_t node, vd_t old_node, encoded_iterator<typename treedec_traits<T_t>::bag_type::iterator> encIt){
        auto powIt = bag(node, _t).begin();

        if(encIt == bag(old_node, _t).end()){
            return 0;
        }
        
        unsigned number = 0;
        unsigned summand = 1;

        for(; encIt != bag(old_node, _t).end();){
            if(*encIt == *powIt){
                number += summand;
                ++encIt;
            }
            powIt++;
            summand <<= 1;
        }

        return number;
    }

    //different supersets and new vertex (same as e = encode(new_node, old_node, encIt); update_encoding(new_node, e, new_vertex))
    unsigned encode(vd_t node, vd_t old_node, encoded_iterator<typename treedec_traits<T_t>::bag_type::iterator> encIt, vd_t new_vertex){
        auto powIt = bag(node, _t).begin();
        
        unsigned number = 0;
        unsigned summand = 1;
        bool found = false;

        for(; encIt != bag(old_node, _t).end();){
            if(*encIt == *powIt){
                number += summand;
                ++encIt;
            }
            else if(!found && new_vertex == *powIt){
                number += summand;
                found = true;
            }
            powIt++;
            summand <<= 1;
        }

        if(!found){
            for(; powIt != bag(node, _t).end();){
                if(new_vertex == *powIt){
                    number += summand;
                    break;
                }
                powIt++;
                summand <<= 1;
            }
        }

        return number;
    }

    //same superset, but new vertex (find bitpos and set it to 1)
    unsigned update_encoding(vd_t node, Encoded_t key, vd_t new_vertex){
        unsigned number = key;
        unsigned summand = 1;

        auto powIt = bag(node, _t).begin();
        for(; powIt != bag(node, _t).end(); ++powIt){
            if(new_vertex == *powIt){
                number += summand;
                break;
            }
            summand <<= 1;
        }

        return number;
    }

    void decode(vd_t node, Encoded_t key, Decoded_t &result){
        typename Decoded_t::iterator it = bag(node, _t).begin();
        for(unsigned i = 0; i < min(CHAR_BIT*sizeof(key), bag(node, _t).size()); i++){
            if(key & 1){
                result.insert(*it);
            }
            ++it;
            key >>= 1;
        }
    }

    unsigned get_size(vd_t node, Encoded_t key){
        unsigned s = 0;
        for(unsigned i = 0; i < min(CHAR_BIT*sizeof(key), bag(node, _t).size()); i++){
            if(key & 1){
                s++;
            }
            key >>= 1;
        }
        return s;
    }

    T_t &_t;
    typename std::vector<std::map<Encoded_t, Value_t> > _results;

}; //Intermediate_Results


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
void top_down_computation2(T_t &T,
                    typename boost::graph_traits<T_t>::vertex_descriptor cur,
                    treedec::app::detail::Intermediate_Results<T_t> &iRes,
                    unsigned int val, typename treedec_traits<T_t>::bag_type &solution,
                    unsigned parents_choice=0, bool minimizing=true)
{
    treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);
    treedec::nice::enum_node_type parent_type = treedec::nice::get_type_parent(cur, T);

    if(parent_type == treedec::nice::INVALID){
        assert(bag(cur, T).size() == 0);

        typename boost::graph_traits<T_t>::vertex_descriptor child =
                               *(boost::adjacent_vertices(cur, T).first);

        top_down_computation2(T, child, iRes, val, solution, parents_choice, minimizing);
        return;
    }

    typename boost::graph_traits<T_t>::vertex_descriptor parent =
                           boost::source(*(boost::in_edges(cur, T).first), T);

    if(node_type == treedec::nice::LEAF){
        if(val == 1){
             solution.insert(*bag(cur, T).begin());
        }
        return;
    }
    else if(parent_type == treedec::nice::FORGET){
        unsigned new_vertex =
                                 treedec::nice::get_forgotten_vertex(parent, T);


        encoded_iterator<typename treedec_traits<T_t>::bag_type::iterator> encIt(parents_choice, bag(parent, T).begin(), bag(parent, T).end());
        
        unsigned encoded_without = iRes.encode(cur, parent, encIt);
        unsigned encoded_with = iRes.encode(cur, parent, encIt, new_vertex);

        int val_without = iRes.get(cur, encoded_without);
        int val_with = iRes.get(cur, encoded_with);

        if(val_with == -1){
            parents_choice = encoded_without;
        }
        else if(val_without == -1){
            parents_choice = encoded_with;
            solution.insert(new_vertex);
        }
        else{
            if(val_with <= val_without){
                parents_choice = minimizing? encoded_with : encoded_without;
            }
            else{
                parents_choice = minimizing? encoded_without : encoded_with;
            }
        }

        if(node_type == treedec::nice::JOIN){
            auto adjIt = boost::adjacent_vertices(cur, T).first;

            auto child1 = *(adjIt++);
            auto child2 = *adjIt;

            top_down_computation2(T, child1, iRes, val, solution, parents_choice, minimizing);
            top_down_computation2(T, child2, iRes, val, solution, parents_choice, minimizing);
            return;

        }
        else{
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                   *(boost::adjacent_vertices(cur, T).first);

            top_down_computation2(T, child, iRes, val, solution, parents_choice, minimizing);
        }
    }




    else if(parent_type == treedec::nice::INTRODUCE){
        unsigned new_vertex = treedec::nice::get_introduced_vertex(parent, T);

        encoded_iterator<typename treedec_traits<T_t>::bag_type::iterator> encIt(parents_choice, bag(parent, T).begin(), bag(parent, T).end());
        
        unsigned encoded_without = iRes.encode(cur, parent, encIt);
        unsigned encoded_with = iRes.encode(cur, parent, encIt, new_vertex);

        int val_without = iRes.get(cur, encoded_without);
        int val_with = iRes.get(cur, encoded_with);


//TODO: choice, add to solution..

        if(node_type == treedec::nice::JOIN){
            auto adjIt = boost::adjacent_vertices(cur, T).first;

            auto child1 = *(adjIt++);
            auto child2 = *adjIt;

            top_down_computation2(T, child1, iRes, val, solution, parents_choice, minimizing);
            top_down_computation2(T, child2, iRes, val, solution, parents_choice, minimizing);
            return;

        }
        else{
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                   *(boost::adjacent_vertices(cur, T).first);

            top_down_computation2(T, child, iRes, val, solution, parents_choice, minimizing);
        }

        incomplete();
    }


    else if(parent_type == treedec::nice::JOIN){
        auto adjIt = boost::adjacent_vertices(cur, T).first;

        auto child1 = *(adjIt++);
        auto child2 = *adjIt;

        top_down_computation2(T, child1, iRes, val, solution, parents_choice, minimizing);
        top_down_computation2(T, child2, iRes, val, solution, parents_choice, minimizing);
        return;
    }
}



template <typename T_t>
void top_down_computation2_old(T_t &T,
                    typename boost::graph_traits<T_t>::vertex_descriptor cur,
                    treedec::app::detail::Intermediate_Results<T_t> &iRes,
                    unsigned int val, typename treedec_traits<T_t>::bag_type &S,
                    typename treedec_traits<T_t>::bag_type &S_comp,
                    typename treedec_traits<T_t>::bag_type subset,
                    unsigned int take_flag, unsigned parents_choice=0)
{
    treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

    if(node_type == treedec::nice::LEAF){
        if(val == 1){
             S.insert(*bag(cur, T).begin());
        }
        return;
    }
    else if(node_type == treedec::nice::INTRODUCE){
        //TODO: necessary? subset could be the encoded one here
        unsigned subset_encoded = iRes.encode(cur, subset);

        if(take_flag == 1){ //child of a join node
            val = iRes.get(cur, subset_encoded);
        }

        unsigned encoded;

        for(typename std::map<unsigned, int>::iterator it =
                   iRes._results[cur].begin(); it != iRes._results[cur].end(); it++)
        {
            encoded = it->first;
            typename treedec_traits<T_t>::bag_type decoded_set;
            iRes.decode(cur, encoded, decoded_set);

            bool condition = false;

            if(take_flag == 0 && it->second == (int)val){
                std::cout << "case 1" << std::endl;
                condition = true;
            }
            else if(take_flag == 1 && decoded_set == subset){
                std::cout << "case 2" << std::endl;
                condition = true;
            }
            else if(take_flag == 2 && it->second == (int)val && std::includes(decoded_set.begin(), decoded_set.end(), subset.begin(), subset.end())){
                std::cout << "case 3" << std::endl;
                condition = true;
            }

            if(condition){
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

        if(S.find(treedec::nice::get_introduced_vertex(cur, T)) != S.end()){
            val = val - 1;
        }

        typename boost::graph_traits<T_t>::vertex_descriptor child =
                               *(boost::adjacent_vertices(cur, T).first);


        subset.erase(treedec::nice::get_introduced_vertex(cur, T));
        top_down_computation2_old(T, child, iRes, val, S, S_comp, subset, 2, encoded);
    }
    else if(node_type == treedec::nice::FORGET){
        //TODO: necessary? subset could be the encoded one here
        unsigned subset_encoded = iRes.encode(cur, subset);

        if(take_flag == 1){ //child of a join node
            val = iRes.get(cur, subset_encoded);
        }

        unsigned encoded;

        for(typename std::map<unsigned, int>::iterator it =
                   iRes._results[cur].begin(); it != iRes._results[cur].end(); it++)
        {
            encoded = it->first;
            typename treedec_traits<T_t>::bag_type decoded_set;
            iRes.decode(cur, encoded, decoded_set);

            bool condition = false;

            if(take_flag == 0 && it->second == (int)val){
                std::cout << "case 1" << std::endl;
                condition = true;
            }
            else if(take_flag == 1 && decoded_set == subset){
                std::cout << "case 2" << std::endl;
                condition = true;
            }
            else if(take_flag == 2 && it->second == (int)val && std::includes(decoded_set.begin(), decoded_set.end(), subset.begin(), subset.end())){
                std::cout << "case 3" << std::endl;
                condition = true;
            }

            if(condition){
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

        typename boost::graph_traits<T_t>::vertex_descriptor child =
                               *(boost::adjacent_vertices(cur, T).first);

        top_down_computation2_old(T, child, iRes, val, S, S_comp, subset, 2, encoded);
    }
    else if(node_type == treedec::nice::JOIN){
        for(typename std::map<unsigned, int>::iterator it =
                          iRes._results[cur].begin(); it != iRes._results[cur].end(); it++)
        {

/* why not?!
            if(it->second != (int)val){
                continue;
            }
*/

            unsigned encoded = it->first;
            typename treedec_traits<T_t>::bag_type decoded_set;
            iRes.decode(cur, encoded, decoded_set);

            bool condition = false;

            if(take_flag == 0 && it->second == (int)val){
                std::cout << "case 1" << std::endl;
                condition = true;
            }
            else if(take_flag == 1 && decoded_set == subset){
                std::cout << "case 2" << std::endl;
                condition = true;
            }
            else if(take_flag == 2 && it->second == (int)val && std::includes(decoded_set.begin(), decoded_set.end(), subset.begin(), subset.end())){
                std::cout << "case 3" << std::endl;
                condition = true;
            }

            if(condition){
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


                    auto adjIt = boost::adjacent_vertices(cur, T).first;

                    auto child1 = *(adjIt++);
                    auto child2 = *adjIt;

                    top_down_computation2_old(T, child1, iRes, it->second, S, S_comp, must_take, 1, encoded);
                    top_down_computation2_old(T, child2, iRes, it->second, S, S_comp, must_take, 1, encoded);
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
