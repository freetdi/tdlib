// Lukas Larisch, 2014 - 2017
//
// (c) 2014-2016 Goethe-Universität Frankfurt
// (c) 2017      King Abdullah University of Science and Technology
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

#ifndef TD_INDEPENDENT_SET
#define TD_INDEPENDENT_SET

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

} //namespace detail (for max_clique)


template <typename G_t, typename T_t>
unsigned int max_clique_with_treedecomposition(G_t &G, T_t &T,
                               typename treedec_traits<T_t>::bag_type &global_result)
{
    unsigned int max = 0;

    typename boost::graph_traits<T_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(T); vIt != vEnd; vIt++){
        //We wouldn't find a larger clique.
        if(bag(*vIt, T).size() <= max){ continue; }

        //Search for a clique of size at least 'size' by inspecting all subsets
        //of size exactly 'size' for size = max+1,max+2,..
        for(unsigned int size = max+1; size <= bag(*vIt, T).size(); size++){
            BOOST_AUTO(P, make_subsets_range(bag(*vIt, T).begin(), bag(*vIt, T).end(), size, size));
            BOOST_AUTO(I, P.first);
            bool changed = false;

            for(; I != bag(*vIt, T).end(); ++I){
                if(treedec::app::detail::is_clique(G, (*I).first, (*I).second)){
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

    // assert(treedec::validation::is_valid_clique(G, result));

    return max;
}

} //namespace app

} //namespace treedec

#endif //TD_INDEPENDENT_SET

// vim:ts=8:sw=4:et
