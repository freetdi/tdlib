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
 * - void max_independent_set_with_treedecomposition(G_t&, T_t&, typename treedec_traits<T_t>::bag_type &result)
 *
 * IMPORT NOTE: ensure that the input treedecomposition is directed by
 *              using treedec::make_rooted(undir_t, dir_t)
 *
 */

#ifndef TREEDEC_INDEPENDENT_SET2_HPP
#define TREEDEC_INDEPENDENT_SET2_HPP

#include "applications.hpp"

namespace treedec{

namespace app{



/* MAX INDEPENDENT SET */

namespace detail{

template <typename G_t, typename T_t>
unsigned int bottom_up_computation_independent_set2(G_t &G, T_t &T,
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
            //For all results S of the child: Store S extended by the introduced vertex with value
            //old_value + 1, if S extended by the introduced vertex is an independent set, else -1.
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                      *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor new_vertex =
                                     treedec::nice::get_introduced_vertex(cur, T);

            for(typename std::map<unsigned, int>::iterator it =
                         iRes._results[child].begin(); it != iRes._results[child].end(); it++)
            {
                bool extensible = true;

                unsigned old_encoded = it->first;
                encoded_iterator<typename treedec_traits<T_t>::bag_type::iterator> encIt(old_encoded, bag(child, T).begin(), bag(child, T).end());

                for(; encIt != bag(child, T).end(); ++encIt){
                    if(boost::edge(*encIt, new_vertex, G).second){
                        extensible = false;
                        break;
                    }
                }

                encoded_iterator<typename treedec_traits<T_t>::bag_type::iterator> encIt2(old_encoded, bag(child, T).begin(), bag(child, T).end());
                unsigned new_encoded1 = iRes.encode(cur, child, encIt2);

                encoded_iterator<typename treedec_traits<T_t>::bag_type::iterator> encIt3(old_encoded, bag(child, T).begin(), bag(child, T).end());
                unsigned new_encoded2 = iRes.encode(cur, child, encIt3, new_vertex);

                if(extensible){
                    iRes.add(cur, new_encoded2, it->second + 1);
                }
                else{
                    iRes.add(cur, new_encoded2, -1);
                }

                iRes.add(cur, new_encoded1, it->second);
            }
        }
        else if(node_type == treedec::nice::FORGET){
            //Store the maximum of S and S extended by the forgotten vertex for each subset S
            //of the current bag.
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                      *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                                treedec::nice::get_forgotten_vertex(cur, T);

            BOOST_AUTO(subset_iters_it, make_subsets_range(bag(cur, T).begin(), bag(cur, T).end(), 0, bag(cur, T).size()).first);

            for(; subset_iters_it != bag(cur, T).end(); ++subset_iters_it){
//                typename treedec_traits<T_t>::bag_type subset((*subset_iters_it).first, (*subset_iters_it).second);

//                typename treedec_traits<T_t>::bag_type new_set = subset;
//                new_set.insert(forgotten_vertex);


                unsigned old_encoded1 = iRes.encode(child, (*subset_iters_it).first, (*subset_iters_it).second);

                unsigned old_encoded2 = iRes.update_encoding(child, old_encoded1, forgotten_vertex);
//                unsigned old_encoded2 = iRes.encode(child, new_set);

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
                        iRes.add(cur, new_encoded, val_without);
                    }
                    else{
                        iRes.add(cur, new_encoded, val_with);
                    }
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


            for(typename std::map<unsigned, int>::iterator it =
                         iRes._results[child1].begin(); it != iRes._results[child1].end(); it++)
            {
                unsigned encoded = it->first;

                if(it->second < 0 || iRes.get(child2, encoded) < 0){
                    iRes.add(cur, encoded, -1);
                }
                else{
                    unsigned subset_size = iRes.get_size(cur, encoded);
                    iRes.add(cur, encoded, it->second + iRes.get(child2, encoded) - subset_size);
                }
            }
        }
    }

    typename boost::graph_traits<T_t>::vertex_descriptor root = find_root(T);
    return (unsigned int) iRes.get(root, 0);
}

} //namespace detail


template <typename G_t, typename T_t>
unsigned int max_independent_set_with_treedecomposition2(G_t &G, T_t &T,
                      typename treedec_traits<T_t>::bag_type &global_result, bool certificate=true)
{

    assert(treedec::is_undirected_type(G));
    assert(treedec::is_edge_set_type(G));
    assert(treedec::no_loops(G));

    assert(treedec::is_valid_treedecomposition(G, T));

    if(boost::num_edges(G) == 0){
        if(boost::num_vertices(G) > 0){
            global_result.insert(boost::vertices(G).first, boost::vertices(G).second);
            return boost::num_vertices(G);
        }
        else{
            return 0;
        }
    }

//    std::vector<std::map<typename treedec_traits<T_t>::bag_type, int> > results(boost::num_vertices(T));

    treedec::app::detail::Intermediate_Results<T_t> iRes(T);

//    unsigned int max = treedec::app::detail::bottom_up_computation_independent_set2(G, T, results);

    unsigned int max = treedec::app::detail::bottom_up_computation_independent_set2(G, T, iRes);


    if(certificate && max > 0){
        typename treedec_traits<T_t>::bag_type a, b;
        typename boost::graph_traits<T_t>::vertex_descriptor root = find_root(T);
//        treedec::app::detail::top_down_computation(T, root, results, max, global_result, a, b, 0);
        treedec::app::detail::top_down_computation2(T, root, iRes, max, global_result, a, b, 0);  
    }

    assert(treedec::validation::is_valid_independent_set(G, global_result));

    return max;
}

} //namespace app

} //namespace treedec

#endif //TREEDEC_INDEPENDENT_SET2_HPP

// vim:ts=8:sw=4:et
