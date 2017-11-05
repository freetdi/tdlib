// Lukas Larisch, 2014 - 2017
// Felix Salfelder 2016 - 2017
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


#ifndef TREEDEC_PREPROCESSING_EDGE_INSERTION_HPP
#define TREEDEC_PREPROCESSING_EDGE_INSERTION_HPP

#include <vector>
#include <set>
#include <boost/graph/copy.hpp>

#include "algo.hpp"
#include "marker.hpp"
#include "config_traits.hpp"
#include "graph.hpp"

#include "lower_bounds.hpp"
#include "preprocessing.hpp"

namespace treedec{

template <typename G_t>
void edge_insertion_neighbors(G_t & G, unsigned ub_bagsize){
    assert(ub_bagsize > 0);

    treedec::k_neighbour_improved_graph(G, ub_bagsize-1);
    //call BothSimplicial until no reduction is possible
}

template <typename G_t>
void edge_insertion_paths(G_t & G, unsigned ub_bagsize){
    assert(ub_bagsize > 0);

    treedec::k_path_improved_graph(G, ub_bagsize-1);
    //call BothSimplicial until no reduction is possible
}






} //namespace treedec

#endif //guard
