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
// Offers functionality to compute upper bounds on tree width
//

#ifndef TD_UPPER_BOUNDS
#define TD_UPPER_BOUNDS

#include <vector>
#include <boost/graph/adjacency_list.hpp>

namespace treedec{

namespace ub{

template <typename G_t>
unsigned int _minDegree(G_t &G){
    unsigned int upper_bound = 0;

    while(boost::num_edges(G) > 0){
        //Search a minimum degree vertex.
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        boost::tie(vIt, vEnd) = boost::vertices(G);
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex = *vIt;
        unsigned int min_degree = boost::num_vertices(G);

        for(; vIt != vEnd; vIt++){
            unsigned int degree = boost::out_degree(*vIt, G);
            if(degree != 0 && degree < min_degree){
                min_degree = degree;
                min_vertex = *vIt;
            }
        }

        if(min_degree > upper_bound){
            upper_bound = min_degree;
        }

        //Collect the neighbours of min_vertex.
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> neighbours(boost::out_degree(min_vertex, G));

        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;

        unsigned int i = 0;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, G); nIt != nEnd; nIt++){
            neighbours[i++] = *nIt;
        }

        //Make the neighbours of min_vertex a clique.
        for(unsigned int i = 0; i < neighbours.size(); i++){
            for(unsigned int j = i+1; j < neighbours.size(); j++){
                boost::add_edge(neighbours[i], neighbours[j], G);
            }
        }

        boost::clear_vertex(min_vertex, G);
    }

    return upper_bound;
}

template <typename G_t>
unsigned int minDegree(G_t G){
    return _minDegree(G);
}

} //namespace ub

} //namespace treedec

#endif //ifdef TD_UPPER_BOUNDS

// vim:ts=8:sw=4:et
