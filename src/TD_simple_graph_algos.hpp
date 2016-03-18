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
// vertex/edge deletion and other simple operations on graphs.
//

#ifndef TD_SIMPLE_GRAPH_ALGOS
#define TD_SIMPLE_GRAPH_ALGOS

namespace treedec{

#include <set>
#include <boost/graph/adjacency_list.hpp>
#include "TD_noboost.hpp"

//Provides a mapping (version for bags).
template <typename G_t>
void induced_subgraph(G_t &H, G_t &G,
                      typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &X,
                      typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &vdMap)
{
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> internal_map(boost::num_vertices(G));
    std::vector<bool> disabled(boost::num_vertices(G), true);
    vdMap.resize(X.size());

    for(typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type::iterator sIt
         = X.begin(); sIt != X.end(); sIt++)
    {
       unsigned int pos1 = noboost::get_pos(*sIt, G);
       internal_map[pos1] = boost::add_vertex(H);
       disabled[pos1] = false;

       unsigned int pos2 = noboost::get_pos(internal_map[pos1], H);
       vdMap[pos2] = *sIt;
    }

    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        unsigned int spos=noboost::get_pos(boost::source(*eIt, G), G);
        unsigned int dpos=noboost::get_pos(boost::target(*eIt, G), G);
        if(!disabled[spos] && !disabled[dpos]){
            boost::add_edge(internal_map[spos], internal_map[dpos], H);
        }
    }
}

//Provides a mapping (version for set<descriptor>)
template <typename G_t>
void induced_subgraph(G_t &H, G_t &G,
                      typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &X,
                      typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &vdMap)
{
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> internal_map(boost::num_vertices(G));
    std::vector<bool> disabled(boost::num_vertices(G), true);
    vdMap.resize(X.size());

    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt
         = X.begin(); sIt != X.end(); sIt++)
    {
       unsigned int pos1 = noboost::get_pos(*sIt, G);
       internal_map[pos1] = boost::add_vertex(H);
       disabled[pos1] = false;

       unsigned int pos2 = noboost::get_pos(internal_map[pos1], H);
       vdMap[pos2] = *sIt;
    }

    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        unsigned int spos=noboost::get_pos(boost::source(*eIt, G), G);
        unsigned int dpos=noboost::get_pos(boost::target(*eIt, G), G);
        if(!disabled[spos] && !disabled[dpos]){
            boost::add_edge(internal_map[spos], internal_map[dpos], H);
        }
    }
}

//Provides a mapping (version for bags).
template <typename G_t>
void induced_subgraph(G_t &H, G_t &G,
                      typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &X)
{
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> internal_map(boost::num_vertices(G));
    std::vector<bool> disabled(boost::num_vertices(G), true);

    for(typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type::iterator sIt
         = X.begin(); sIt != X.end(); sIt++)
    {
       unsigned int pos1 = noboost::get_pos(*sIt, G);
       internal_map[pos1] = boost::add_vertex(H);
       disabled[pos1] = false;
    }

    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        unsigned int spos=noboost::get_pos(boost::source(*eIt, G), G);
        unsigned int dpos=noboost::get_pos(boost::target(*eIt, G), G);
        if(!disabled[spos] && !disabled[dpos]){
            boost::add_edge(internal_map[spos], internal_map[dpos], H);
        }
    }
}

template <typename G_t>
bool is_edge_between_sets(G_t &G,
                 typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &X,
                 typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &Y)
{
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
void get_neighbourhood(G_t &G, std::vector<bool> &disabled,
             std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &X,
             std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &S_X)
{
    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = X.begin(); sIt != X.end(); sIt++){
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*sIt, G); nIt != nEnd; nIt++){
           unsigned int pos = noboost::get_pos(*nIt, G);
           if(!disabled[pos] && X.find(*nIt) == X.end()){
               S_X.insert(*nIt);
           }
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
    unsigned int pos = noboost::get_pos(vertex, G); std::cout << "pos2: " << pos << std::endl;
    visited[pos] = true;
    typename boost::graph_traits<G_t>::adjacency_iterator  nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(vertex, G); nIt != nEnd; nIt++){
        unsigned int npos = noboost::get_pos(*nIt, G); std::cout << "npos: " << npos << std::endl;
        if(!visited[npos]){
            components[comp_idx].insert(*nIt);
            t_search_components(G, *nIt, visited, components, comp_idx);
        }
    }
}


template <typename G_t>
void get_components_provided_map(G_t &G,
             std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > &components,
             std::vector<bool> &visited){

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    int comp_idx = -1;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int pos = noboost::get_pos(*vIt, G); std::cout << "pos: " << pos << std::endl;
        if(!visited[pos]){
            components.resize(components.size()+1);
            comp_idx++;

            components[comp_idx].insert(*vIt);
            t_search_components(G, *vIt, visited, components, comp_idx);
        }
    }
}

template <typename G_t>
void get_components(G_t &G,
             std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > &components)
{
    std::vector<bool> visited(boost::num_vertices(G), false);
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    int comp_idx = -1;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int pos = noboost::get_pos(*vIt, G);
        if(!visited[pos]){
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
        unsigned int pos = noboost::get_pos(*vIt, G);
        max = (pos > max)? pos : max;
    }

    idxMap.resize(max+1);

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int pos = noboost::get_pos(*vIt, G);
        idxMap[pos] = *vIt;
    }
}

} //namespace treedec

#endif //TD_SIMPLE_GRAPH_ALGOS

// vim:ts=8:sw=4:et
