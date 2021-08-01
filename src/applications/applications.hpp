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
// Foundation, 51 Franklin Street - Suite 500, Boston, MA 02110-1335, USA.
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


namespace treedec{

namespace app{

enum optimization_type { MINIMIZING, MAXIMIZING };

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
     : _num(number), _last(0), _supIt(superIt), _supEnd(superEnd) {
        if(_num == 0){
             _supIt = _supEnd;
        }else{
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
        if(_results[node].find(key) == _results[node].end()){
            return -1;
        }
        return _results[node][key];
    }

    bool exists(vd_t node, Encoded_t key){
        return _results[node].find(key) != _results[node].end();
    }

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

    unsigned encode_more(vd_t node, vd_t old_node, encoded_iterator<typename treedec_traits<T_t>::bag_type::iterator> encIt, vd_t new_vertex){
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

    unsigned encode_less(vd_t node, vd_t old_node, encoded_iterator<typename treedec_traits<T_t>::bag_type::iterator> encIt, vd_t skip_vertex){
        auto powIt = bag(node, _t).begin();
        
        unsigned number = 0;
        unsigned summand = 1;

        for(; encIt != bag(old_node, _t).end();){
            if(*encIt == skip_vertex){
                ++encIt;
                continue;
            }

            if(*encIt == *powIt){
                number += summand;
                ++encIt;
            }
            powIt++;
            summand <<= 1;
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


/* The top-down computation on tree decompositions.
 *  Traces back the decisions that have been made in bottom-up run
 *  Running time O(n*(width(T)+1))
 */

template <typename T_t>
void top_down_computation(T_t &T,
                    typename boost::graph_traits<T_t>::vertex_descriptor cur,
                    treedec::app::detail::Intermediate_Results<T_t> &iRes,
                    unsigned int val, typename treedec_traits<T_t>::bag_type &solution,
                    unsigned parents_choice, optimization_type opt_type)
{
    treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);
    treedec::nice::enum_node_type parent_type = treedec::nice::get_type_parent(cur, T);
 
    if(parent_type == treedec::nice::INVALID){
        assert(bag(cur, T).size() == 0);

        typename boost::graph_traits<T_t>::vertex_descriptor child =
                               *(boost::adjacent_vertices(cur, T).first);

        top_down_computation(T, child, iRes, val, solution, parents_choice, opt_type);
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
        unsigned new_vertex = treedec::nice::get_forgotten_vertex(parent, T);

        encoded_iterator<typename treedec_traits<T_t>::bag_type::iterator> encIt(parents_choice, bag(parent, T).begin(), bag(parent, T).end());
        
        unsigned encoded_without = iRes.encode(cur, parent, encIt);
        unsigned encoded_with = iRes.encode_more(cur, parent, encIt, new_vertex);

        int val_without = iRes.get(cur, encoded_without);
        int val_with = iRes.get(cur, encoded_with);

        if(val_with == -1){
            parents_choice = encoded_without;
            val = val_without;
        }
        else if(val_without == -1){
            parents_choice = encoded_with;
            val = val_with;
            solution.insert(new_vertex);
        }
        else{
            if(val_with <= val_without){
                if(opt_type==treedec::app::MINIMIZING){
                    parents_choice = encoded_with;
                    solution.insert(new_vertex);
                    val = val_with;
                }
                else{
                    parents_choice = encoded_without;
                    val = val_without;
                }
            }
            else{
                if(opt_type==treedec::app::MINIMIZING){
                    parents_choice = encoded_without;
                    val = val_without;
                }
                else{
                    parents_choice = encoded_with;
                    solution.insert(new_vertex);
                    val = val_with;
                }
            }
        }
    }
    else if(parent_type == treedec::nice::INTRODUCE){
        unsigned new_vertex = treedec::nice::get_introduced_vertex(parent, T);

        encoded_iterator<typename treedec_traits<T_t>::bag_type::iterator> encIt(parents_choice, bag(parent, T).begin(), bag(parent, T).end());
        unsigned encoded_without = iRes.encode_less(cur, parent, encIt, new_vertex);

        int val_without = iRes.get(cur, encoded_without);

        val = val_without;

        parents_choice = encoded_without;
    }
    else if(parent_type == treedec::nice::JOIN){
    }


    if(node_type == treedec::nice::INTRODUCE){
            unsigned intro_vertex = treedec::nice::get_introduced_vertex(cur, T);
            val = (solution.find(intro_vertex) != solution.end())? val-1u : val;
    }

    if(node_type == treedec::nice::JOIN){
        auto adjIt = boost::adjacent_vertices(cur, T).first;

        auto child1 = *(adjIt++);
        auto child2 = *adjIt;

        top_down_computation(T, child1, iRes, iRes.get(child1, parents_choice), solution, parents_choice, opt_type);
        top_down_computation(T, child2, iRes, iRes.get(child2, parents_choice), solution, parents_choice, opt_type);
    }
    else{
        typename boost::graph_traits<T_t>::vertex_descriptor child =
                               *(boost::adjacent_vertices(cur, T).first);

        top_down_computation(T, child, iRes, val, solution, parents_choice, opt_type);
    }
}

} //namespace detail

} //namespace app

} //namespace treedec

#endif //TREEDEC_APPLICATIONS_HPP

// vim:ts=8:sw=4:et
