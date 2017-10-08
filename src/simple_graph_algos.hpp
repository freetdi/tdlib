// Lukas Larisch, 2014 - 2016
// Felix Salfelder, 2016
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

/*
 * Offers functionality to compute the connected components of a graph,
 * vertex/edge deletion and other simple operations on graphs.
 */

#ifndef TD_SIMPLE_GRAPH_ALGOS
#define TD_SIMPLE_GRAPH_ALGOS

#include <set>
#include <limits>

#include <boost/graph/adjacency_list.hpp>
#include "graph.hpp"
#include "marker_util.hpp"
#include "induced_subgraph.hpp"

namespace treedec{

//This function is used in the minimalChordal algorithm.
template <typename G_t, typename E_t>
void induced_subgraph_omit_edges(G_t &H, const G_t &G,
                      typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &X,
                      E_t &edges,
                      typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &vdMap)
{
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> internal_map(boost::num_vertices(G));
    std::vector<BOOL> disabled(boost::num_vertices(G), true);
    vdMap.resize(X.size());

    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt
         = X.begin(); sIt != X.end(); sIt++)
    {
       auto pos1=boost::get(boost::vertex_index, G, *sIt);
       internal_map[pos1] = boost::add_vertex(H);
       disabled[pos1] = false;

       auto pos2=boost::get(boost::vertex_index, H, internal_map[pos1]);
       vdMap[pos2] = *sIt;
    }

    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        auto s=boost::source(*eIt, G);
        auto t=boost::target(*eIt, G);
        auto spos=boost::get(boost::vertex_index, G, s);
        auto tpos=boost::get(boost::vertex_index, G, t);
        if(!disabled[spos] && !disabled[tpos]){
            bool omit = false;
            for(unsigned int i = 0; i < edges.size(); i++){
                if((edges[i].first == boost::source(*eIt, G)
                 && edges[i].second == boost::target(*eIt, G))
                || (edges[i].second == boost::source(*eIt, G)
                 && edges[i].first == boost::target(*eIt, G)))
                {
                    omit = true;
                    break;
                }
            }
            if(!omit){
                boost::add_edge(internal_map[spos], internal_map[tpos], H);
            }else{
            }
        }else{
        }
    }
}

// todo: prototype in graph.hpp
template <typename H_t, typename G_t, class S_t, class M_t>
void copy_induced_subgraph(H_t &H, G_t const &G, S_t const& X, M_t* vdMap)
{
    assert(boost::num_vertices(H)==0);
    typedef typename boost::graph_traits<G_t>::vertex_descriptor G_vertex_descriptor;
    typedef typename boost::graph_traits<H_t>::vertex_iterator H_vertex_iterator;
    std::vector<G_vertex_descriptor> internal_map(boost::num_vertices(G));
    std::vector<BOOL> disabled(boost::num_vertices(G), true);
    if(vdMap){
        vdMap->resize(X.size());
    }
    H = MOVE(G_t(X.size()));
    H_vertex_iterator h=boost::vertices(H).first;
    size_t pos2=0;

    for(typename S_t::const_iterator sIt=X.begin(); sIt!=X.end(); ++sIt){
       unsigned int pos1 = get_pos(*sIt, G);
       internal_map[pos1] = *h;
       ++h;
       disabled[pos1] = false;

       if(vdMap){
           (*vdMap)[pos2] = *sIt;
           ++pos2;
       }
    }
    assert(h==boost::vertices(H).second);

    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt!=eEnd; ++eIt){
        unsigned int spos=get_pos(boost::source(*eIt, G), G);
        unsigned int dpos=get_pos(boost::target(*eIt, G), G);
        if(!disabled[spos] && !disabled[dpos]){
            boost::add_edge(internal_map[spos], internal_map[dpos], H);
        }
    }
}

// paste subgraph of G induced by X into H.
// store map V(H) -> X \subset V(G) in vdMap.
//
// TODO: misleading name. paste_induced_subgraph?
template <typename G_t, class S_t, class M_t>
void induced_subgraph(G_t &H, G_t const &G, S_t const& X, M_t* vdMap)
{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    if(boost::num_vertices(H)==0){
        return copy_induced_subgraph(H, G, X, vdMap);
    }
    else{
        throw exception_invalid_precondition();
    }
    std::vector<vertex_descriptor> internal_map(boost::num_vertices(G));
    std::vector<BOOL> disabled(boost::num_vertices(G), true);
    if(vdMap){
        vdMap->resize(X.size());
    }

    for(typename S_t::const_iterator sIt=X.begin(); sIt!=X.end(); ++sIt){
       unsigned int pos1 = get_pos(*sIt, G);
       internal_map[pos1] = boost::add_vertex(H);
       disabled[pos1] = false;

       if(vdMap){
           unsigned int pos2 = get_pos(internal_map[pos1], H);
           (*vdMap)[pos2] = *sIt;
       }
    }

    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt!=eEnd; ++eIt){
        unsigned int spos=get_pos(boost::source(*eIt, G), G);
        unsigned int dpos=get_pos(boost::target(*eIt, G), G);
        if(!disabled[spos] && !disabled[dpos]){
            boost::add_edge(internal_map[spos], internal_map[dpos], H);
        }
    }
}

template <typename G_t, class S_t, class M_t>
void induced_subgraph(G_t &H, G_t const &G, S_t const& X, M_t& vdMap)
{
    return induced_subgraph(H, G, X, &vdMap);
}

// same, but without vdMap
template <typename G_t, class S_t, class M_t>
void induced_subgraph(G_t &H, G_t const &G, S_t const& X)
{
    return induced_subgraph(H, G, X, NULL);
}

// TODO: different containers?
template <typename G_t, typename i1, typename i2>
std::pair<typename boost::graph_traits<G_t>::edge_descriptor, bool>
   edge(i1 Xit, i1 Xend, i2 Yit, i2 Yend, G_t const& G)
{
    for(; Xit!=Xend; ++Xit){
        for(; Yit!=Yend; ++Yit){
            BOOST_AUTO(P, boost::edge(*Xit, *Yit, G));
            if(P.second){
                return P;
            }
        }
    }
    return std::make_pair(typename boost::graph_traits<G_t>::edge_descriptor(),false);
}

// wrapper, use begin()/end()
template <typename G_t, typename vertex_set>
bool is_edge_between_sets(G_t &G, vertex_set const& X, vertex_set const& Y)
{
    return edge(G, X.begin(), X.end(), Y.begin(), Y.end()).second;
}

template <typename G_t, typename It>
inline void get_neighbourhood(G_t const &G, std::vector<BOOL> &disabled,
             It Xit, It Xend,
             std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &S_X)
{
    for(;Xit!=Xend; ++Xit){
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*Xit, G); nIt != nEnd; nIt++){
           unsigned int pos = get_pos(*nIt, G);

           if(!disabled[pos]){
               S_X.insert(*nIt);
           }
        }
    }
}

// wrapper
// vertex into S_X if
//  - it is adjacent to a vertex in X
//  - if it is not disabled (by position)
//  - if it is not an element of X
//
template <typename G_t>
inline void get_neighbourhood(G_t &G, std::vector<BOOL> &disabled,
             std::set<typename boost::graph_traits<G_t>::vertex_descriptor> const &X,
             std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &S_X)
{
    return get_neighbourhood(G, disabled, X.begin(), X.end(), S_X);
}

template <typename G_t, typename B_t>
void t_search_components(G_t const &G,
        typename boost::graph_traits<G_t>::vertex_descriptor vertex,
        std::vector<BOOL> &visited,
        std::vector<B_t> &components,
        int comp_idx)
{
    unsigned int pos = get_pos(vertex, G);
    visited[pos] = true;
    typename boost::graph_traits<G_t>::adjacency_iterator  nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(vertex, G); nIt != nEnd; nIt++){
        unsigned int npos = get_pos(*nIt, G);
        if(!visited[npos]){
            components[comp_idx].insert(*nIt);
            t_search_components(G, *nIt, visited, components, comp_idx);
        }
    }
}

template <typename G_t, typename VB_t>
void get_components_provided_map(G_t const &G,
             VB_t &components,
             std::vector<BOOL> &visited){

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    int comp_idx = -1;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int pos = get_pos(*vIt, G);
        if(!visited[pos]){
            components.resize(components.size()+1);
            comp_idx++;

            components[comp_idx].insert(*vIt);
            t_search_components(G, *vIt, visited, components, comp_idx);
        }
    }
}

#if 0 // duplicates?

template <typename G_t, typename B_t>
void t_search_components(G_t &G,
        typename boost::graph_traits<G_t>::vertex_descriptor vertex,
        std::vector<BOOL> &visited,
        std::vector<B_t> &components,
        int comp_idx)
{
    unsigned int pos = get_pos(vertex, G);
    visited[pos] = true;
    typename boost::graph_traits<G_t>::adjacency_iterator  nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(vertex, G); nIt != nEnd; nIt++){
        unsigned int npos = get_pos(*nIt, G);
        if(!visited[npos]){
            components[comp_idx].insert(*nIt);
            t_search_components(G, *nIt, visited, components, comp_idx);
        }
    }
}

template <typename G_t>
void get_components_provided_map(G_t &G,
             std::vector<typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type> &components,
             std::vector<BOOL> &visited){

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    int comp_idx = -1;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int pos = get_pos(*vIt, G);
        if(!visited[pos]){
            components.resize(components.size()+1);
            comp_idx++;

            components[comp_idx].insert(*vIt);
            t_search_components(G, *vIt, visited, components, comp_idx);
        }
    }
}
#endif

template <typename G_t>
void get_components(G_t &G,
             std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > &components)
{
    std::vector<BOOL> visited(boost::num_vertices(G), false);
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    int comp_idx = -1;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        auto pos=boost::get(boost::vertex_index, G, *vIt);
        if(!visited[pos]){
            components.resize(components.size()+1);
            comp_idx++;

            components[comp_idx].insert(*vIt);
            t_search_components(G, *vIt, visited, components, comp_idx);
        }
    }
}

// find a neighbour x of v with the least common neighbours
template <typename G_t, class M>
inline typename boost::graph_traits<G_t>::vertex_descriptor
   get_least_common_vertex(const typename boost::graph_traits<G_t>::vertex_descriptor v,
           M& marker, const G_t &G)
{
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;

    auto min_common=std::numeric_limits<vertices_size_type>::max();

    marker.clear();

    auto p=boost::adjacent_vertices(v, G);
    typename boost::graph_traits<G_t>::vertex_descriptor w = *p.first;
    mark_range(p.first, p.second, marker);

    auto q=boost::adjacent_vertices(v, G);
    for(; q.first != q.second; q.first++){
        vertices_size_type cnt_common=0;
        auto p=boost::adjacent_vertices(*q.first, G);
        for(; p.first!=p.second; ++p.first){
            if(marker.is_marked(*p.first)){
                cnt_common++;
            }else{
            }
        }
        if(cnt_common < min_common){
            w = *q.first;
            min_common = cnt_common;
        }
    }

    return w;
}

} //namespace treedec

#endif //TD_SIMPLE_GRAPH_ALGOS

// vim:ts=8:sw=4:et
