// Lukas Larisch, 2014 - 2016
//
// (c) 2014-2016 Goethe-Universität Frankfurt
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

#include "TD_elimination_orderings.hpp"

namespace treedec{

namespace ub{

// remove later (underscore prefix, to be considered detail/private)
template <typename G_t>
unsigned int _minDegree(G_t &G)
{
    return treedec::impl::minDegree_decomp(G);
}

template <typename G_t>
unsigned int minDegree(G_t& G){
    return _minDegree(G);
}

// remove later (underscore prefix, to be considered detail/private)
template <typename G_t>
unsigned int _fillIn(G_t &G)
{
    return treedec::impl::fillIn_decomp(G);
}

template <typename G_t>
unsigned int fillIn(G_t& G){
    return _fillIn(G);
}

} //namespace ub

} //namespace treedec

#endif //TD_UPPER_BOUNDS

// vim:ts=8:sw=4:et
