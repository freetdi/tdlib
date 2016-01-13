// Lukas Larisch, 2014 - 2015
//
// (c) 2014-2015 Goethe-Universität Frankfurt
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
// Offers functionality to solve hard problems with help of treedecompositions
//

#define HAVE_CLIQUER

#ifndef TD_APPLICATIONS
#define TD_APPLICATIONS

#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include "TD_hard_problems.hpp"


namespace treedec{

namespace app{

template <typename G_t, typename T_t>
void max_clique_with_treedecomposition(G_t &G, T_t &T, std::vector<unsigned int> &global_result, treedec::np::max_clique_base<G_t> &mclb){
    for(unsigned int i = 0; i < boost::num_vertices(T); i++){
        G_t H;
        induced_subgraph(H, G, T[i].bag);
        std::vector<unsigned int> result;
        mclb.max_clique(H, result);
        if(result.size() > global_result.size())
            global_result = result;
    }
}

template <typename G_t, typename T_t>
void max_independent_set_with_treedecomposition(G_t &G, T_t &T, std::vector<unsigned int> &global_result, treedec::np::max_independent_set_base<G_t> &misb){
}

template <typename G_t, typename T_t>
void min_vertex_cover_with_treedecomposition(G_t &G, T_t &T, std::vector<unsigned int> &global_result, treedec::np::min_vertex_cover_base<G_t> &mvcb){
}

template <typename G_t, typename T_t>
void min_coloring_with_treedecomposition(G_t &G, T_t &T, std::vector<std::vector<unsigned int> > &global_result, treedec::np::min_coloring_base<G_t> &mcob){
}

template <typename G_t, typename T_t>
void min_dominating_set_with_treedecomposition(G_t &G, T_t &T, std::vector<unsigned int> &global_result, treedec::np::min_dominating_set_base<G_t> &mvcb){
}

} //namespace app

} //namespace treedec

#endif


