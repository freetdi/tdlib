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

#ifndef TREEDEC_VERTEX_COVER_HPP
#define TREEDEC_VERTEX_COVER_HPP

#include <cmath>
#include <climits>

#include "applications.hpp"
#include "../validation.hpp"

namespace treedec{

namespace app{


/* MIN VERTEX COVER */

namespace detail{

#if 0
//TODO: more efficient
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
#endif

//TODO: marker?
template <typename G_t>
bool is_valid_extension(G_t &G,
    typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type &bag,
    const typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type &old_VC,
    typename boost::graph_traits<G_t>::vertex_descriptor new_vertex)
{
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(new_vertex, G); nIt != nEnd; nIt++){
        if(bag.find(*nIt) != bag.end()){ //restriction to the current bag
            if(old_VC.find(*nIt) == old_VC.end()){ //the edge {new_vertex, *nIt} is not covered by old_VC
                return false;
            }
        }
    }

    return true;
}

/* TODO
 *   
 *  -currently, '-1' entries are stored. In that case, one can use a vector instead of a map.
 *   then, the encoding of a set is the pos in the vector. also, using subset iters is not
 *   necessary, just count from 0 to 2^bagsize
 *
 *  -otherwise, get rid of the '-1' entriess
*/

template <typename G_t, typename T_t>
unsigned int bottom_up_computation_vertex_cover(G_t &G, T_t &T,
         treedec::app::detail::Intermediate_Results<T_t> &iRes)
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
            iRes.add(cur, 0, 0);
            iRes.add(cur, 1, 1);
        }
        else if(node_type == treedec::nice::INTRODUCE){
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                     *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor new_vertex =
                             treedec::nice::get_introduced_vertex(cur, T);

            //check if old VC is still valid, if new_vertex is added to the considered graph
            for(typename std::map<unsigned, int>::iterator it =
                         iRes._results[child].begin(); it != iRes._results[child].end(); it++)
            {
                unsigned old_encoded = it->first;
                typename treedec_traits<T_t>::bag_type decoded_set; //TODO: get rid of this
                iRes.decode(child, old_encoded, decoded_set);
                unsigned new_encoded = iRes.encode(cur, decoded_set);

                if(is_valid_extension(G, bag(cur, T), decoded_set, new_vertex)){
                    iRes.add(cur, new_encoded, iRes.get(child, old_encoded));
                }
                else{
                    iRes.add(cur, new_encoded, -1);
                }
            }

            //check for new set in new graph
            for(typename std::map<unsigned, int>::iterator it =
                    iRes._results[child].begin(); it != iRes._results[child].end(); it++)
            {
                unsigned old_encoded = it->first;

                encoded_iterator<typename treedec_traits<T_t>::bag_type::iterator> encIt(old_encoded, bag(child, T).begin(), bag(child, T).end());
                unsigned new_encoded = iRes.encode_more(cur, child, encIt, new_vertex);

                if(it->second != -1){
                    iRes.add(cur, new_encoded, it->second + 1);
                }
                else{ //TODO: is this really correct?!
/*
                    typename treedec_traits<T_t>::bag_type new_set;
                    iRes.decode(child, old_encoded, new_set);
                    new_set.insert(new_vertex);
                    if(is_vertex_cover2(G, bag(cur, T), new_set)){
                        iRes.add(cur, new_encoded, new_set.size());
                    }
                    else{
*/
                        iRes.add(cur, new_encoded, -1);
//                    }
                }
            }
        }
        else if(node_type == treedec::nice::FORGET){
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                     *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                          treedec::nice::get_forgotten_vertex(cur, T);

            BOOST_AUTO(subset_iters_it, make_subsets_range(bag(cur, T).begin(), bag(cur, T).end(), 0, bag(cur, T).size()).first);

            for(; subset_iters_it != bag(cur, T).end(); ++subset_iters_it){
                unsigned old_encoded1 = iRes.encode(child, (*subset_iters_it).first, (*subset_iters_it).second);
                unsigned old_encoded2 = iRes.update_encoding(child, old_encoded1, forgotten_vertex);

                int val_without = iRes.get(child, old_encoded1);
                int val_with = iRes.get(child, old_encoded2);

                unsigned new_encoded = iRes.encode(cur, (*subset_iters_it).first, (*subset_iters_it).second);

                if(val_with == -1){
                    iRes.add(cur, new_encoded, val_without);
                }
                else if(val_without == -1){
                    iRes.add(cur, new_encoded, val_with);
                }
                else{
                    if(val_with <= val_without){
                        iRes.add(cur, new_encoded, val_with);
                    }
                    else{
                        iRes.add(cur, new_encoded, val_without);
                    }
                }
            }
        }
        else if(node_type == treedec::nice::JOIN){
            typename boost::graph_traits<T_t>::vertex_descriptor child1 =
                                     *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<T_t>::vertex_descriptor child2 =
                                     *(++boost::adjacent_vertices(cur, T).first);

            for(typename std::map<unsigned, int>::iterator it =
                    iRes._results[child1].begin(); it != iRes._results[child1].end(); it++)
            {
                unsigned encoded = it->first;

                if(iRes.get(child1, encoded) < 0 || iRes.get(child2, encoded) < 0){
                    iRes.add(cur, encoded, -1);
                }
                else{
                    unsigned subset_size = iRes.get_size(child1, encoded);
                    iRes.add(cur, encoded, iRes.get(child1, encoded) + iRes.get(child2, encoded) - subset_size);
                }
            }
        }
    }

    typename boost::graph_traits<T_t>::vertex_descriptor root = find_root(T);
    return (unsigned int) iRes.get(root, 0);
}

} //namespace detail (min_vertex_cover)



template <typename G_t, typename T_t>
unsigned int min_vertex_cover_with_treedecomposition(G_t &G, T_t &T,
              typename treedec_traits<T_t>::bag_type &global_result, bool certificate=true)
{

    assert(treedec::is_undirected_type(G));
    assert(treedec::is_edge_set_type(G));
    assert(treedec::no_loops(G));

    assert(treedec::is_valid_treedecomposition(G, T));


    treedec::app::detail::Intermediate_Results<T_t> iRes(T);

    unsigned int min = treedec::app::detail::bottom_up_computation_vertex_cover(G, T, iRes);

    if(certificate && min > 0){
        typename boost::graph_traits<T_t>::vertex_descriptor root = find_root(T);
        treedec::app::detail::top_down_computation2(T, root, iRes, min, global_result, 0, treedec::app::MINIMIZING);
    }

    assert(certificate && treedec::validation::is_valid_vertex_cover(G, global_result));

    return min;
}




} //namespace app

} //namespace treedec

#endif //TREEDEC_VERTEX_COVER_HPP

// vim:ts=8:sw=4:et
