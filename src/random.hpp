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
// Foundation, 51 Franklin Street - Suite 500, Boston, MA 02110-1335, USA.
//
//

#ifndef TREEDEC_RANDOM_HPP
#define TREEDEC_RANDOM_HPP

#include <vector>
#include <algorithm>    // std::random_shuffle

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace treedec{

namespace random{

template <typename G_t>
int try_some_elimination_orderings(G_t &G,
       std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &best_elim_ordering,
       unsigned int count = 10)
{
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> elim_ordering(boost::num_vertices(G));
    for(unsigned int i=0; i<elim_ordering.size(); ++i){
        elim_ordering[i] = i;
    }

    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > elimination_orderings(count);
    for(unsigned int i=0; i<elimination_orderings.size(); ++i){
        //Use the built-in random generator.
        std::random_shuffle(elim_ordering.begin(), elim_ordering.end());
        elimination_orderings[i] = elim_ordering;
    }

    int min_width = INT_MAX;

#ifdef OPENMP
    #pragma omp parallel for
#endif
    for(unsigned int i = 0; i < count; i++){
        G_t H;
        boost::copy_graph(G, H);
        int width_i = treedec::get_width_of_elimination_ordering(H, elimination_orderings[i]);
        if(width_i < min_width){
            min_width = width_i;
            best_elim_ordering = MOVE(elimination_orderings[i]);
        }
    }

    return min_width;
}

} //namespace random

} //namespace treedec

#endif //TREEDEC_RANDOM_HPP
