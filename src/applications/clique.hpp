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
 * - void max_clique_with_treedecomposition(G_t&, T_t&, typename treedec_traits<T_t>::bag_type &result)
 *
 * IMPORT NOTE: ensure that the input treedecomposition is directed by
 *              using treedec::make_rooted(undir_t, dir_t)
 *
 */

#ifndef TREEDEC_CLIQUE_HPP
#define TREEDEC_CLIQUE_HPP

#include "applications.hpp"

namespace treedec{

namespace app{


/* MAX CLIQUE */

namespace detail{

template <typename G_t, typename A_t, typename B_t>
bool is_clique(G_t &G, A_t A, B_t B){
    BOOST_AUTO(p1, A);
    for(; p1 != B; ++p1){
        BOOST_AUTO(p2, p1);
        p2++;
        for(; p2 != B; ++p2){
            if(!boost::edge(*p1, *p2, G).second){
                return false;
            }
        }
    }
    return true;
}


//for all vertices included in the test set, increment
//a counter for all neighbours
//check if counter is size-1 for all vertices in that set
//at the end.
template <typename G_t, typename A_t, typename B_t>
bool is_clique2(G_t &G, A_t A, B_t B, unsigned size, typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &vdMap){
    std::vector<unsigned> counter(boost::num_vertices(G), 0);

    BOOST_AUTO(vIt, A);
    for(; vIt != B; ++vIt){
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(vdMap[*vIt], G); nIt != nEnd; ++nIt){
            if(boost::degree(*nIt, G) < size-1u){
                return false;
            }

            counter[*nIt]++;
        }
    }

    vIt = A;
    for(; vIt != B; ++vIt){
        if(counter[vdMap[*vIt]] != size-1u){
            return false;
        }
    }


    return true;
}


//for all vertices included in the test set, increment
//a counter for all neighbours
//check if counter is size-1 for all vertices in that set
//at the end.
template <typename G_t, typename A_t, typename B_t>
bool is_clique3(G_t &G, A_t A, B_t B, unsigned size){
    std::vector<unsigned> counter(boost::num_vertices(G), 0);

    BOOST_AUTO(vIt, A);
    for(; vIt != B; ++vIt){
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; ++nIt){
            if(boost::degree(*nIt, G) < size-1u){
                return false;
            }

            counter[*nIt]++;
        }
    }

    vIt = A;
    for(; vIt != B; ++vIt){
        if(counter[*vIt] != size-1u){
            return false;
        }
    }


    return true;
}

} //namespace detail (for max_clique)


template <typename G_t, typename T_t>
unsigned int max_clique_with_treedecomposition(G_t &G, T_t &T,
                               typename treedec_traits<T_t>::bag_type &global_result)
{

    assert(treedec::is_undirected_type(G));
    assert(treedec::is_edge_set_type(G));
    assert(treedec::no_loops(G));

    assert(treedec::is_valid_treedecomposition(G, T));

    if(boost::num_edges(G) == 0){
        if(boost::num_vertices(G) > 0){
            global_result.insert(*(boost::vertices(G).first));
            return 1;
        }
        else{
            return 0;
        }
    }

    unsigned int max = 1;

    typename boost::graph_traits<T_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(T); vIt != vEnd; vIt++){
        //We wouldn't find a larger clique.
        if(bag(*vIt, T).size() <= max){
            continue;
        }


        //there should be a combination of width and n, such that using an induced subgraph is faster..
        bool use_induced_subgraph = false;
        G_t H;
        typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> vdMap;
        typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> vdMap2;

#if 0
        if(bag(*vIt, T).size() > 20 && boost::num_vertices(G) > 100){
            use_induced_subgraph = true;
            treedec::copy_induced_subgraph(H, G, bag(*vIt, T), &vdMap, &vdMap2);
        }
#endif 

        //Search for a clique of size greater than 'size' by inspecting all subsets
        //of size exactly 'size' for size = max+1,max+2,..
        for(unsigned int size = max+1; size <= bag(*vIt, T).size(); size++){
            BOOST_AUTO(P, make_subsets_range(bag(*vIt, T).begin(), bag(*vIt, T).end(), size, size));
            BOOST_AUTO(I, P.first);
            bool changed = false;

            for(; I != bag(*vIt, T).end(); ++I){
                if((use_induced_subgraph && treedec::app::detail::is_clique2(H, (*I).first, (*I).second, size, vdMap2))
                || treedec::app::detail::is_clique(G, (*I).first, (*I).second)){
//                  if(treedec::app::detail::is_clique(G, (*I).first, (*I).second)){
//                  if(treedec::app::detail::is_clique3(G, (*I).first, (*I).second, size)){
                    max = size;

                    global_result.clear();
                    BOOST_AUTO(p, (*I).first);

                    for(; p != (*I).second; p++){
                        global_result.insert(*p);
                    }

                    changed = true;

                    //We wouldn't find a larger clique in this loop.
                    break;
                }
            }
            //This bag doesn't contain a clique larger that 'max'.
            if(!changed){
                break;
            }
        }
    }

    assert(treedec::validation::is_valid_clique(G, global_result));

    return max;
}

} //namespace app

} //namespace treedec

#endif //TREEDEC_CLIQUE_HPP

// vim:ts=8:sw=4:et
