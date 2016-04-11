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

        unsigned int min_degree = UINT_MAX;
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

        noboost::make_clique(boost::adjacent_vertices(min_vertex, G), G);

        boost::clear_vertex(min_vertex, G);
    }

    return upper_bound;
}

template <typename G_t>
unsigned int minDegree(G_t& G){
    return _minDegree(G);
}

template <typename G_t>
unsigned int _minFill(G_t &G){
    unsigned int upper_bound = 0;

    while(boost::num_edges(G) > 0){
        //Search a vertex v such that least edges are missing for making the neighbourhood of v a clique.
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        boost::tie(vIt, vEnd) = boost::vertices(G);
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex = *vIt;

        unsigned int min_fill = UINT_MAX;
        for(; vIt != vEnd; vIt++){
            if(boost::out_degree(*vIt, G) == 0){
                continue;
            }

            unsigned int current_fill = 0;

            typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
            for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(*vIt, G); nIt1 != nEnd; nIt1++){
                nIt2 = nIt1;
                nIt2++;
                for(; nIt2 != nEnd; nIt2++){
                    if(!boost::edge(*nIt1, *nIt2, G).second){
                        current_fill++;
                    }
                }
            }

            if(current_fill < min_fill){
                min_fill = current_fill;
                min_vertex = *vIt;
                if(current_fill == 0){
                    break;
                }
            }
        }

        if(boost::degree(min_vertex, G) > upper_bound){
            upper_bound = boost::degree(min_vertex, G);
        }

        noboost::make_clique(boost::adjacent_vertices(min_vertex, G), G);

        boost::clear_vertex(min_vertex, G);
    }

    return upper_bound;
}

template <typename G_t>
unsigned int minFill(G_t& G)
{
    return _minFill(G);
}

} //namespace ub

} //namespace treedec

#endif //TD_UPPER_BOUNDS

// vim:ts=8:sw=4:et
