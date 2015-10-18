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
        //search a minimum degree vertex
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;
        unsigned int min_degree = boost::num_vertices(G);
    
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int degree = boost::out_degree(*vIt, G);
            if(degree != 0 && degree < min_degree){
                min_degree = degree;
                min_vertex = *vIt;
            }
        }

        if(min_degree > upper_bound)
            upper_bound = min_degree;
    
        //collect the neighbours
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> neighbours;
        
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;    

        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, G); nIt != nEnd; nIt++)
            neighbours.push_back(*nIt);

	//make the neighbours a clique
        for(unsigned int i = 0; i < neighbours.size(); i++){
            for(unsigned int j = i+1; j < neighbours.size(); j++)
                boost::add_edge(neighbours[i], neighbours[j], G);	
        }

        boost::clear_vertex(min_vertex, G);
    }

    return upper_bound;
}

//quite fast if the container for edges is vecS
template <typename G_t>
unsigned int _minDegree_fast(G_t &G){
    unsigned int upper_bound = 0;

    while(boost::num_edges(G) > 0){
        //search a minimum degree vertex
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;
        unsigned int min_degree = boost::num_vertices(G);
    
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int degree = boost::out_degree(*vIt, G);
            if(degree != 0 && degree < min_degree){
                min_degree = degree;
                min_vertex = *vIt;
            }
        }

        if(min_degree > upper_bound)
            upper_bound = min_degree;
    
        //collect the neighbours
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> neighbours;
        
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;    

        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, G); nIt != nEnd; nIt++)
            neighbours.push_back(*nIt);

	//make the neighbours a clique
        for(unsigned int i = 0; i < neighbours.size(); i++){
            for(unsigned int j = i+1; j < neighbours.size(); j++){
                if(!boost::edge(neighbours[i], neighbours[j], G).second)
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

}
}

#endif


