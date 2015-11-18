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

#ifndef TD_APPLICATIONS
#define TD_APPLICATIONS

#include <vector>
#include <boost/graph/adjacency_list.hpp>

namespace treedec{
    
namespace app{

template <typename G_t, typename T_t>
void max_clique(G_t &G, T_t &T){
}

template <typename G_t, typename T_t>
void max_vertex_cover(G_t &G, T_t &T){
}

template <typename G_t, typename T_t>
void max_independent_set(G_t &G, T_t &T){
}

template <typename G_t, typename T_t>
void max_coloring(G_t &G, T_t &T){
}

template <typename G_t, typename T_t>
void hamiltonian_cycle(G_t &G, T_t &T){
}

#endif


