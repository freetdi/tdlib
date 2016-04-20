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

#ifndef TD_APPROXIMATE_SEPERATOR
#define TD_APPROXIMATE_SEPERATOR

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include "simple_graph_algos.hpp"


template <typename G_t>
void approximate_vertex_seperator(G_t &G, std::vector<bool> &disabled, typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &X, typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &S){
}

#endif //ifdef TD_APPROXIMATE_SEPERATOR

// vim:ts=8:sw=4:et
