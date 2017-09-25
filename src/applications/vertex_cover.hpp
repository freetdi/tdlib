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

#ifndef TD_VERTEX_COVER
#define TD_VERTEX_COVER

#include "applications.hpp"

namespace treedec{

namespace app{


/* MIN VERTEX COVER */

namespace detail{

template <typename G_t>
bool is_vertex_cover(G_t &G, 
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
unsigned int bottom_up_computation_vertex_cover(G_t &G, T_t &T,
         std::vector<std::map<typename treedec_traits<T_t>::bag_type, int> > &results)
{
    std::stack<typename boost::graph_traits<T_t>::vertex_descriptor> S;
    treedec::nice::postorder_traversal(T, S);
    typename boost::graph_traits<T_t>::vertex_descriptor cur;

    while(!S.empty()){
        cur = S.top();
        S.pop();

        treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

        if(node_type == treedec::nice::LEAF){
            results[cur][typename treedec_traits<T_t>::bag_type()] = 0;
            results[cur][bag(cur, T)] = 1;
        }
        else if(node_type == treedec::nice::INTRODUCE){
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                     *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor new_vertex =
                             treedec::nice::get_introduced_vertex(cur, T);

            for(typename std::map<typename treedec_traits<T_t>::bag_type, int>::iterator it =
                         results[child].begin(); it != results[child].end(); it++)
            {
                if(is_vertex_cover(G, bag(cur, T), it->first)){
                    results[cur][it->first] = results[child][it->first];
                }
                else{
                    results[cur][it->first] = -1;
                }
            }

            for(typename std::map<typename treedec_traits<T_t>::bag_type, int>::iterator it =
                    results[child].begin(); it != results[child].end(); it++)
            {
                typename treedec_traits<T_t>::bag_type tmp = it->first;
                tmp.insert(new_vertex);

                if(it->second != -1){
                    results[cur][tmp] = results[child][it->first] + 1;
                }
                else{
                    if(is_vertex_cover(G, bag(cur, T), tmp)){
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

            std::vector<typename treedec_traits<T_t>::bag_type> subs;
            treedec::powerset(bag(cur, T), subs);

            for(unsigned int i = 0; i < subs.size(); i++){
                typename treedec_traits<T_t>::bag_type tmp = subs[i];
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

            std::vector<typename treedec_traits<T_t>::bag_type> subs;
            treedec::powerset(bag(cur, T), subs);

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
              typename treedec_traits<T_t>::bag_type &global_result)
{
    std::vector<std::map<typename treedec_traits<T_t>::bag_type, int> > results(boost::num_vertices(T));

    unsigned int max = treedec::app::detail::bottom_up_computation_vertex_cover(G, T, results);

    if(max > 0){
        typename treedec_traits<T_t>::bag_type a, b;
        typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
        treedec::app::detail::top_down_computation(T, root, results, max, global_result, a, b, 0);
    }

    // assert(treedec::validation::is_valid_vertex_cover(G, result));

    return max;
}

} //namespace app

} //namespace treedec

#endif //TD_VERTEX_COVER

// vim:ts=8:sw=4:et
