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

#include <set>
#include <boost/graph/adjacency_list.hpp>
#include "TD_noboost.hpp"

template <typename G_t>
void delete_edges(G_t &G, std::vector<std::vector<unsigned int> > &edges){
    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(unsigned int i = 0; i < edges.size(); i++){
        for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
            unsigned sid=noboost::get_id(G, boost::source(*eIt, G));
            unsigned tid=noboost::get_id(G, boost::target(*eIt, G));
            if((sid == edges[i][0] && tid == edges[i][1])
             ||(sid == edges[i][1] && tid == edges[i][0])){
                boost::remove_edge(boost::source(*eIt, G), boost::target(*eIt, G), G);
                break;
            }
        }
    }
}

template <typename G_t>
void induced_subgraph(G_t &H, G_t &G, std::set<unsigned int> &X){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap(boost::num_vertices(G));
    for(std::set<unsigned int>::iterator sIt = X.begin(); sIt != X.end(); sIt++){
       idxMap[*sIt] = boost::add_vertex(H);
       H[idxMap[*sIt]].id = *sIt;
    }

    std::vector<bool> disabled(boost::num_vertices(G), true);
    for(std::set<unsigned int>::iterator sIt = X.begin(); sIt != X.end(); sIt++){
        disabled[*sIt] = false;
    }

    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        unsigned sid=noboost::get_id(G, boost::source(*eIt, G));
        unsigned tid=noboost::get_id(G, boost::target(*eIt, G));
        if(!disabled[sid] && !disabled[tid]){
            boost::add_edge(idxMap[sid], idxMap[tid], H);
        }
    }
}

template <typename G_t>
bool is_edge_between_sets(G_t &G, typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &X, typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &Y){
    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt1 = X.begin(); sIt1 != X.end(); sIt1++){
        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt2 = Y.begin(); sIt2 != Y.end(); sIt2++){
            if(boost::edge(*sIt1, *sIt2, G).second){
                return true;
            }
        }
    }
    return false;
}

template <typename G_t>
void get_neighbourhood(G_t &G, std::vector<bool> &disabled, std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &X, std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &S_X){
    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = X.begin(); sIt != X.end(); sIt++){
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*sIt, G); nIt != nEnd; nIt++){
           unsigned id = noboost::get_id(G, *nIt);
           if(!disabled[id] && X.find(*nIt) == X.end()){
               S_X.insert(*nIt);
           }
        }
    }
}

template <typename G_t>
void t_search_components(G_t &G, typename boost::graph_traits<G_t>::vertex_descriptor vertex, std::vector<bool> &visited, std::vector<std::set<unsigned int> > &components, int comp_idx){
    unsigned id = noboost::get_id(G, vertex);
    visited[id] = true;
    std::vector<unsigned int> N;
    typename boost::graph_traits<G_t>::adjacency_iterator  nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(vertex, G); nIt != nEnd; nIt++){
        unsigned nid = noboost::get_id(G, *nIt);
        if(!visited[nid]){
            components[comp_idx].insert(nid);
            t_search_components(G, *nIt, visited, components, comp_idx);
        }
    }
}

template <typename G_t>
void t_search_components(G_t &G,
        typename boost::graph_traits<G_t>::vertex_descriptor vertex,
        std::vector<bool> &visited,
        std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > &components,
        int comp_idx)
{
    unsigned id=noboost::get_id(G, vertex);
    visited[id] = true;
    std::vector<unsigned int> N;
    typename boost::graph_traits<G_t>::adjacency_iterator  nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(vertex, G); nIt != nEnd; nIt++){
        unsigned nid=noboost::get_id(G, *nIt);
        if(!visited[nid]){
            components[comp_idx].insert(*nIt);
            t_search_components(G, *nIt, visited, components, comp_idx);
        }
    }
}

template <typename G_t>
void get_components(G_t &G, std::vector<std::set<unsigned int> > &components){
    unsigned int max = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned id = noboost::get_id(G, *vIt);
        max = (id > max)? id : max;
    }

    std::vector<bool> visited(max+1);

    for(unsigned int i = 0; i < max+1; i++){
        visited[i] = false;
    }

    int comp_idx = -1;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned id = noboost::get_id(G, *vIt);
        if(!visited[id]){
            components.resize(components.size()+1);
            comp_idx++;

            components[comp_idx].insert(id);
            t_search_components(G, *vIt, visited, components, comp_idx);
        }
    }
}


template <typename G_t>
void get_components_provided_map(G_t &G, std::vector<std::set<unsigned int> > &components, std::vector<bool> &visited){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    int comp_idx = -1;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned id = noboost::get_id(G, *vIt);
        if(!visited[id]){
            components.resize(components.size()+1);
            comp_idx++;

            components[comp_idx].insert(id);
            t_search_components(G, *vIt, visited, components, comp_idx);
        }
    }
}

template <typename G_t>
void get_components_provided_map(G_t &G, std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > &components, std::vector<bool> &visited){  
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    int comp_idx = -1;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned id=noboost::get_id(G, *vIt);
        if(!visited[id]){
            components.resize(components.size()+1);
            comp_idx++;

            components[comp_idx].insert(*vIt);
            t_search_components(G, *vIt, visited, components, comp_idx);
        }
    }
}

template <typename G_t>
void make_index_map(G_t &G, std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &idxMap){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    unsigned int max = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned id=noboost::get_id(G, *vIt);
        max = (id > max)? id : max;
    }

    idxMap.resize(max+1);

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned id=noboost::get_id(G, *vIt);
        idxMap[id] = *vIt;
    }
}

#endif //ifdef TD_SIMPLE_GRAPH_ALGOS

// vim:ts=8:sw=4:et
