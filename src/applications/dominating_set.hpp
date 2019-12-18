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
 * Provides following functions (namespace treedec::app):
 *
 * - void min_dominating_set_with_treedecomposition(G_t&, T_t&, typename treedec_traits<T_t>::bag_type &result)
 *
 * IMPORT NOTE: ensure that the input treedecomposition is directed by
 *              using treedec::make_rooted(undir_t, dir_t)
 *
 */

#ifndef TREEDEC_DOMINATING_SET_HPP
#define TREEDEC_DOMINATING_SET_HPP

#include "applications.hpp"

#define TREEDEC_GET_POS(a,b) ( boost::get(boost::vertex_index, b, a) )

namespace treedec{

namespace app{


/* MIN DOMINATING SET */

namespace detail{

static void all_k_colorings(unsigned int /*n*/, unsigned int k,
        std::set<unsigned int> &M, std::vector<std::vector<int> > &colorings,
        std::vector<int> &pattern)
{
    if(M.size() == 0){
        return;
    }

    std::vector<int> coloring(pattern);
    std::set<unsigned int>::iterator iM = M.begin();
    while(iM != M.end()){
        coloring[*(iM++)]++;
    }

    iM = M.begin();
    unsigned int c = 0;

    colorings[c++] = coloring;

    while(iM != M.end() && c < colorings.size()){
        if(coloring[*iM] < (int)k-1){
            if(iM == M.end()){ break; }
            coloring[*iM]++;

            colorings[c++] = coloring;
        }
        else{
            while(coloring[*iM] == (int)k-1 && iM != M.end()){
                coloring[*iM] = 0;
                iM++;
            }
            if(iM == M.end()){ break; }

            coloring[*iM]++;

            colorings[c++] = coloring;

            iM = M.begin();
        }
    }

    colorings.resize(c);
}

static void all_k_colorings(unsigned int n, unsigned int k,
        std::set<unsigned int> &M, std::vector<std::vector<int> > &colorings)
{
    std::vector<int> pattern(n, -1);
    all_k_colorings(n, k, M, colorings, pattern);
}

static unsigned int pow(unsigned int b, unsigned int n){
    unsigned int res = b;
    for(unsigned int i = 0; i < n-1; i++){
        res = res * b;
    }
    return res;
}

template <typename G_t, typename T_t>
unsigned int bottom_up_computation_dominating_set(G_t &G, T_t &T,
          std::vector<std::map<std::vector<int>, boost::tuple<int, std::vector<int>, std::vector<int> > > > &results,
          typename std::map<unsigned int, typename boost::graph_traits<G_t>::vertex_descriptor> &inv_map)
{
    std::stack<typename boost::graph_traits<T_t>::vertex_descriptor> S;
    treedec::nice::postorder_traversal(T, S);
    typename boost::graph_traits<T_t>::vertex_descriptor cur;

    while(!S.empty()){
        cur = S.top();
        S.pop();

        treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

        if(node_type == treedec::nice::LEAF){
            auto leaf=*(bag(cur, T).begin());
            auto pos=boost::get(boost::vertex_index, G, leaf);

            std::vector<int> result(boost::num_vertices(G), -1);
            for(unsigned int i = 0; i < 3; i++){
                result[pos] = i;
                results[cur][result] = boost::make_tuple((i == 2) ? 1 : (i == 1) ? 0 : -1, std::vector<int>(), std::vector<int>());
            }
        }
        else if(node_type == treedec::nice::INTRODUCE){
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                         *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor new_vertex =
                                treedec::nice::get_introduced_vertex(cur, T);
            auto pos=boost::get(boost::vertex_index, G, new_vertex);

            //If x has a neighbour in the current bag that is dominating .. in the coloring C, store the coloring
            //C' formed by coloring according to C and coloring x as dominated by a vertex.
            for(std::map<std::vector<int>, boost::tuple<int, std::vector<int>, std::vector<int> > >::iterator it
                         = results[child].begin(); it != results[child].end(); it++)
            {
                typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
                bool applied = false;
                boost::tie(nIt, nEnd) = boost::adjacent_vertices(new_vertex, G);
                for(; nIt != nEnd; nIt++){
                    unsigned int posn = TREEDEC_GET_POS(*nIt, G);
                    if(it->first[posn] == 2){
                        std::vector<int> result(it->first);
                        result[pos] = 0;
                        results[cur][result] = boost::tuple<int, std::vector<int>, std::vector<int> >
                                                 (boost::get<0>(it->second), it->first, std::vector<int>());
                        applied = true;
                        break;
                    }
                }
                if(!applied){
                    std::vector<int> result(it->first);
                    result[pos] = 0;
                    results[cur][result] = boost::tuple<int, std::vector<int>, std::vector<int> >
                                                 (-1, it->first, std::vector<int>());
                }
            }

            //Add the coloring C' formed by coloring according to C and coloring x as dominating.
            for(std::map<std::vector<int>, boost::tuple<int, std::vector<int>, std::vector<int> > >::iterator it
                     = results[child].begin(); it != results[child].end(); it++)
            {
                std::vector<int> result(it->first);
                result[pos] = 2;
                std::vector<int> result2(result);
                for(unsigned int i = 0; i < result2.size(); i++){
                    if(result2[i] == 0){
                        if(boost::edge(new_vertex, inv_map[i], G).second){
                            result2[i] = 1;
                        }
                    }
                }

                result2[pos] = -1;

                if(boost::get<0>(results[child][result2]) == -1){
                    results[cur][result] = boost::tuple<int, std::vector<int>, std::vector<int> >
                                                 (-1, it->first, std::vector<int>());
                }
                else{
                    results[cur][result] = boost::tuple<int, std::vector<int>, std::vector<int> >
                                                 (boost::get<0>(results[child][result2])+1,
                                                  result2, std::vector<int>());
                }

            }

            //Add the coloring C' formed by coloring according to C and coloring x as not known yet.
            for(std::map<std::vector<int>, boost::tuple<int, std::vector<int>, std::vector<int> > >::iterator it =
                     results[child].begin(); it != results[child].end(); it++)
            {
                std::vector<int> result(it->first);
                result[pos] = 1;
                results[cur][result] = boost::make_tuple(boost::get<0>(results[child][it->first]),
                        it->first, std::vector<int>());
            }
        }
        else if(node_type == treedec::nice::FORGET){
            typename boost::graph_traits<T_t>::vertex_descriptor child =
                                             *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                                       treedec::nice::get_forgotten_vertex(cur, T);
            unsigned int pos = TREEDEC_GET_POS(forgotten_vertex, G);

            for(std::map<std::vector<int>, boost::tuple<int, std::vector<int>, std::vector<int> > >::iterator it
                        = results[child].begin(); it != results[child].end(); it++)
            {
                if(it->first[pos] == 2){
                    int val1 = it->second.get<0>();
                    std::vector<int> result(it->first);

                    std::vector<int> choice1(result);
                    result[pos] = 0;
                    std::vector<int> choice2(result);

                    int val2 = boost::get<0>(results[child][result]);
                    result[pos] = -1;

                    if(val1 == -1){
                        results[cur][result] = boost::make_tuple(val2, choice2, std::vector<int>());
                    }
                    else if(val2 == -1){
                        results[cur][result] = boost::make_tuple(val1, choice1, std::vector<int>());
                    }
                    else if(val1 < val2){
                        results[cur][result] = boost::make_tuple(val1, choice1, std::vector<int>());
                    }
                    else{
                        results[cur][result] = boost::make_tuple(val2, choice2, std::vector<int>());
                    }
                }
            }
        }
        else if(node_type == treedec::nice::JOIN){
            typename boost::graph_traits<T_t>::vertex_descriptor child1 =
                                         *(boost::adjacent_vertices(cur, T).first);

            typename boost::graph_traits<T_t>::vertex_descriptor child2 =
                                         *(++boost::adjacent_vertices(cur, T).first);

            std::set<unsigned int> M;
            for(typename treedec_traits<T_t>::bag_type::iterator bIt =
                        bag(cur, T).begin(); bIt != bag(cur, T).end(); bIt++)
            {
                unsigned int pos = TREEDEC_GET_POS(*bIt, G);
                M.insert(pos);
            }
            std::vector<std::vector<int> > colorings(pow(3, M.size()));
            all_k_colorings(boost::num_vertices(G), 3, M, colorings);

            std::vector<std::vector<int> > modified_colorings;

            for(unsigned int i = 0; i < colorings.size(); i++){
                std::set<unsigned int> L;
                for(unsigned int j = 0; j < colorings[i].size(); j++){
                    if(colorings[i][j] == 0){
                        L.insert(j);
                    }
                }

                std::vector<std::vector<int> > colorings2;
                unsigned int s = (L.size() == 0)? 0 : pow(2, L.size());
                if(s > 0){
                    colorings2.resize(s);
                }

                all_k_colorings(boost::num_vertices(G), 2, L, colorings2);

                if(L.size() == 0){
                    colorings2.push_back(colorings[i]);
                }

                for(unsigned int u = 0; u < colorings2.size(); u++){
                    for(unsigned int v = 0; v < colorings2[u].size(); v++){
                        if(colorings[i][v] == 1){ colorings2[u][v] = 1; }
                        else if(colorings[i][v] == 2){ colorings2[u][v] = 2; }
                    }
                }

                unsigned int minimum = UINT_MAX;
                unsigned int min_g = 0;
                unsigned int min_h = 0;
                unsigned int ones = std::count(colorings[i].begin(), colorings[i].end(), 2);

                for(unsigned int g = 0; g < colorings2.size(); g++){
                    for(unsigned int h = 0; h < colorings2.size(); h++){
                        //colorings[i][l] == 0 -> colorings2[g][l] == 0 or colorings2[h][l] == 0
                        bool invalid = false;
                        for(unsigned int l = 0; l < colorings2[g].size(); l++){
                            if(colorings[i][l] == 0 && colorings2[g][l] != 0 && colorings2[h][l] != 0){
                                invalid = true;
                                break;
                            }
                        }
                        if(invalid){
                            continue;
                        }

                        int val_g = boost::get<0>(results[child1][colorings2[g]]);
                        int val_h = boost::get<0>(results[child2][colorings2[h]]);

                        if(val_g != -1 && val_h != -1){
                            if(val_g + val_h - ones < minimum){
                                min_g = g;
                                min_h = h;
                                minimum = val_g + val_h - ones;
                            }
                        }
                    }
                }

                results[cur][colorings[i]] = boost::tuple<int, std::vector<int>, std::vector<int> >
                                                 (minimum, colorings2[min_g], colorings2[min_h]);
            }
        }
    }


    typename boost::graph_traits<T_t>::vertex_descriptor root = find_root(T);
    return (unsigned int) boost::get<0>(results[root].begin()->second);
}


template <typename G_t, typename T_t>
void top_down_computation_min_dominating_set(G_t &G, T_t &T,
              typename boost::graph_traits<T_t>::vertex_descriptor cur,
              std::vector<std::map<std::vector<int>,
                                   boost::tuple<int, std::vector<int>, std::vector<int> > > > &results,
              typename treedec_traits<T_t>::bag_type &global_result,
              std::vector<int> &have_to_take)
{
    treedec::nice::enum_node_type node_type = treedec::nice::get_type(cur, T);

    if(node_type == treedec::nice::LEAF){
        //Nothing to do.
    }
    else if(node_type == treedec::nice::INTRODUCE){
        typename boost::graph_traits<T_t>::vertex_descriptor child =
                                 *(boost::adjacent_vertices(cur, T).first);

//        typename boost::graph_traits<G_t>::vertex_descriptor introduced_vertex =
                                 treedec::nice::get_introduced_vertex(cur, T);

        std::vector<int> next_htt(boost::get<1>(results[cur][have_to_take]));

        top_down_computation_min_dominating_set(G, T, child, results, global_result, next_htt);
    }
    else if(node_type == treedec::nice::FORGET){
        typename boost::graph_traits<T_t>::vertex_descriptor child =
                                 *(boost::adjacent_vertices(cur, T).first);

        typename boost::graph_traits<G_t>::vertex_descriptor forgotten_vertex =
                                 treedec::nice::get_forgotten_vertex(cur, T);

        unsigned int pos = TREEDEC_GET_POS(forgotten_vertex, G);

        std::vector<int> htt(boost::get<1>(results[cur][have_to_take]));
        int assignment = htt[pos];

        if(assignment == 2){
            global_result.insert(forgotten_vertex);
        }

        top_down_computation_min_dominating_set(G, T, child, results, global_result, htt);
    }
    else if(node_type == treedec::nice::JOIN){
        typename boost::graph_traits<T_t>::vertex_descriptor lchild =
                                   *(boost::adjacent_vertices(cur, T).first);
        typename boost::graph_traits<T_t>::vertex_descriptor rchild =
                                   *(++boost::adjacent_vertices(cur, T).first);

        std::vector<int> htt1(boost::get<1>(results[cur][have_to_take]));
        std::vector<int> htt2(boost::get<2>(results[cur][have_to_take]));

        top_down_computation_min_dominating_set(G, T, lchild, results, global_result, htt1);
        top_down_computation_min_dominating_set(G, T, rchild, results, global_result, htt2);
    }
}

} //namespace detail (min_dominating_set)


template <typename G_t, typename T_t>
unsigned int min_dominating_set_with_treedecomposition(G_t &G, T_t &T,
                  typename treedec_traits<T_t>::bag_type &global_result, bool certificate=true)
{
    if(boost::num_vertices(G) == 0){
        return 0;
    }

    typename std::map<unsigned int, typename boost::graph_traits<G_t>::vertex_descriptor> inv_map;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        auto pos=boost::get(boost::vertex_index, G, *vIt);
        inv_map[pos] = *vIt;
    }

    std::vector<std::map<std::vector<int>, boost::tuple<int, std::vector<int>, std::vector<int> > > > results(boost::num_vertices(T));

    unsigned int min = treedec::app::detail::bottom_up_computation_dominating_set(G, T, results, inv_map);

    if(certificate && min > 0){
        std::vector<int> domset(boost::num_vertices(G), -1);
        typename boost::graph_traits<T_t>::vertex_descriptor root = find_root(T);
        std::vector<int> have_to_take(boost::num_vertices(G), -1);
        treedec::app::detail::top_down_computation_min_dominating_set(G, T, root, results, global_result, have_to_take);
    }

    assert(certificate && treedec::validation::is_valid_dominating_set(G, global_result));

    return (unsigned int) min;
}

} //namespace app

} //namespace treedec

#undef TREEDEC_GET_POS

#endif //TREEDEC_DOMINATING_SET_HPP

// vim:ts=8:sw=4:et
