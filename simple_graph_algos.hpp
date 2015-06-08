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
// Offers functionality to compute the connected components of a graph,
// vertex/edge deletion and other simple operations on graphs
//

#ifndef TD_SIMPLE_GRAPH_ALGOS
#define TD_SIMPLE_GRAPH_ALGOS

#include <boost/graph/adjacency_list.hpp>
#include <set>

template <typename G_t>
G_t graph_after_deletion(G_t G, std::set<unsigned int> &X){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> to_delete;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
       if(X.find(G[*vIt].id) != X.end())
           to_delete.push_back(*vIt);
    }
    
    for(unsigned int i = 0; i < to_delete.size(); i++){
        boost::clear_vertex(to_delete[i], G);
        boost::remove_vertex(to_delete[i], G);
    }
    return G;
}


template <typename G_t>
void TD_copy_graph(G_t &G, G_t &H){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    unsigned int max = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        max = (G[*vIt].id > max)? G[*vIt].id : max;

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap(max+1); 
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        idxMap[G[*vIt].id] = boost::add_vertex(H); 
        H[idxMap[G[*vIt].id]].id = G[*vIt].id; 
    }
    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++)
        boost::add_edge(idxMap[G[boost::source(*eIt, G)].id], idxMap[G[boost::target(*eIt, G)].id], H); 
}


template <typename G_t>
void delete_vertices(G_t &G, std::set<unsigned int> &X){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> to_delete;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
       if(X.find(G[*vIt].id) != X.end())
           to_delete.push_back(*vIt);
    }
    
    for(unsigned int i = 0; i < to_delete.size(); i++){
        boost::clear_vertex(to_delete[i], G);
        boost::remove_vertex(to_delete[i], G);
    }
}

template <typename G_t>
void delete_edges(G_t &G, std::vector<std::vector<unsigned int> > &edges){
    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(unsigned int i = 0; i < edges.size(); i++){
        for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
            if((G[boost::source(*eIt, G)].id == edges[i][0] && G[boost::target(*eIt, G)].id == edges[i][1])
             ||(G[boost::source(*eIt, G)].id == edges[i][1] && G[boost::target(*eIt, G)].id == edges[i][0])){
                
                boost::remove_edge(boost::source(*eIt, G), boost::target(*eIt, G), G);
                break;
            }
        }
    }
}
                
template <typename G_t>
G_t get_induced_subgraph(G_t &G, std::set<unsigned int> &X){
    std::set<unsigned int> complement;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        complement.insert(G[*vIt].id);
    for(std::set<unsigned int>::iterator sIt = X.begin(); sIt != X.end(); sIt++)
        complement.erase(*sIt);
    
    return graph_after_deletion(G, complement);
}

template <typename G_t>
bool is_edge_between_sets(G_t &G, std::set<unsigned int> &X, std::set<unsigned int> &Y){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::adjacency_iterator  nIt, nEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(X.find(G[*vIt].id) != X.end()){
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
                if(Y.find(G[*nIt].id) != Y.end())
                    return true;
            }
        }
    }
    return false;
}

template <typename G_t>
void get_neighbourhood(G_t &G, std::set<unsigned int> &X, std::set<unsigned int> &S_X){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::adjacency_iterator  nIt, nEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(X.find(G[*vIt].id) != X.end()){
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
                if(X.find(G[*nIt].id) == X.end())
                    S_X.insert(G[*nIt].id);
            }
        }
    }
}

template <typename G_t> 
void t_search_components(G_t &G, typename boost::graph_traits<G_t>::vertex_descriptor vertex, std::vector<unsigned int> &visited, std::vector<std::set<unsigned int> > &components, int comp_idx){
    visited[G[vertex].id] = true;
    std::vector<unsigned int> N;
    typename boost::graph_traits<G_t>::adjacency_iterator  nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(vertex, G); nIt != nEnd; nIt++){
        if(!visited[G[*nIt].id]){
            components[comp_idx].insert(G[*nIt].id);
            t_search_components(G, *nIt, visited, components, comp_idx);
        }
    }
}

template <typename G_t> 
void get_components(G_t &G, std::vector<std::set<unsigned int> > &components){
    unsigned int max = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        max = (G[*vIt].id > max)? G[*vIt].id : max;
    
    std::vector<unsigned int> visited(max+1);
    
    for(unsigned int i = 0; i < max+1; i++)
        visited[i] = false;
    
    int comp_idx = -1;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(!visited[G[*vIt].id]){
            components.resize(components.size()+1);
            comp_idx++;
            
            components[comp_idx].insert(G[*vIt].id);
            t_search_components(G, *vIt, visited, components, comp_idx);
        }
    }
}

#endif
