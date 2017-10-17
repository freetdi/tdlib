// Lukas Larisch, 2014 - 2016
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
// Offers functionality to compute upper bounds on the tree width of a graph.
//

/*
 These functions are most likely to be interesting for outside use:

 - void minDegree(G_t G&)
 - void fillIn(G_t G&)

*/

#ifndef TD_UPPER_BOUNDS
#define TD_UPPER_BOUNDS

#include "elimination_orderings.hpp"

namespace treedec{

namespace ub{

template <typename G_t>
unsigned int minDegree(G_t& G){
    typedef typename treedec::graph_traits<G_t>::treedec_type T_t;

    if(boost::num_vertices(G) == 0){
        return 0;
    }

    typedef typename std::vector<typename treedec_chooser<G_t>::value_type> O_t;
    T_t T; // dummy

    impl::minDegree<G_t, T_t, O_t> MD(G, &T, (O_t*)NULL, -1u, false);
    MD.do_it();

    return MD.get_bagsize();
}


template <typename G_t>
unsigned int fillIn(G_t& G)
{
    typedef typename treedec::graph_traits<G_t>::treedec_type T_t;

    typedef typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> O_t;
    impl::fillIn<G_t, T_t, O_t> FI(G, (T_t*)NULL, (O_t*)NULL, 0);
    FI.do_it();

    return FI.get_bagsize();
}

} //namespace ub

} //namespace treedec

#endif //TD_UPPER_BOUNDS

// vim:ts=8:sw=4:et
