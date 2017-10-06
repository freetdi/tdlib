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
 * - void min_coloring_with_treedecomposition(G_t&, T_t&, std::vector<typename treedec_traits<T_t>::bag_type> &result)
 *
 * IMPORT NOTE: ensure that the input treedecomposition is directed by
 *              using treedec::make_rooted(undir_t, dir_t)
 *
 */

#ifndef TD_COLORING
#define TD_COLORING

#include "applications.hpp"

namespace treedec{

namespace app{


/* MIN COLORING */

namespace detail{

template <typename G_t>
bool is_valid_extended_coloring(G_t &G, typename boost::graph_traits<G_t>::vertex_descriptor v, std::vector<int> &coloring){
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    auto vpos=boost::get(boost::vertex_index, G, v);
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; nIt++){
        auto npos=boost::get(boost::vertex_index, G, *nIt);
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
                            typename treedec_traits<T_t>::bag_type &bag,
                            std::vector<std::vector<int> > &intersection)
{
    for(unsigned int i = 0; i < results_left.size(); i++){
        for(unsigned int j = 0; j < results_right.size(); j++){
            bool is_compatible = true;
            for(typename treedec_traits<T_t>::bag_type::iterator sIt 
                    = bag.begin(); sIt != bag.end(); sIt++)
            {
                auto pos=boost::get(boost::vertex_index, G, *sIt);
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
            unsigned int leaf = *(bag(cur, T).begin());
            unsigned int pos = get_pos(leaf, G);

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
            unsigned int pos = get_pos(new_vertex, G);

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
            unsigned int pos = get_pos(forgotten_vertex, G);

            std::vector<BOOL> visited(results[child].size(), false);
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

            treedec::app::detail::colorings_intersection<G_t, T_t>(G, results[child1], results[child2], bag(cur, T), results[cur]);

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
                unsigned int pos = get_pos(forgotten_vertex, G);
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
                          std::vector<typename treedec_traits<T_t>::bag_type> &global_result)
{
    std::vector<std::vector<std::vector<int> > > results(boost::num_vertices(T));


    if(boost::num_vertices(G) == 0){
        return 0;
    }
    if(boost::num_edges(G) == 0){
        global_result.resize(1);
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            global_result[0].insert(*vIt);
        }
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
        auto pos=boost::get(boost::vertex_index, G, *vIt);
        inv_map[pos] = *vIt;
    }

    global_result.resize(k);
    for(unsigned int i = 0; i < global_results_map.size(); i++){
        unsigned int col = global_results_map[i];
        global_result[col].insert(inv_map[i]);
    }

    // assert(treedec::validation::is_valid_coloring(G, result));

    return k;
}

} //namespace app

} //namespace treedec

#endif //TD_COLORING

// vim:ts=8:sw=4:et
