// Lukas Larisch, 2014 - 2016
//
// (c) 2014-2016 Goethe-Universität Frankfurt
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
// Nice tree decomposition-related stuff.
//
//
// These functions are most likely to be interesting for outside use:
//
// void nicify(T_t &T)
//

#ifndef TD_NICE_DECOMPOSITION
#define TD_NICE_DECOMPOSITION

#include <set>
#include <vector>
#include <boost/graph/adjacency_list.hpp>

namespace treedec{

namespace nice{

template <typename G_t, typename T_t>
typename boost::graph_traits<G_t>::vertex_descriptor get_introduced_vertex(typename boost::graph_traits<T_t>::vertex_descriptor v, T_t &T,
          std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &idxMap){
    if(T[v].bag.size() == 1){
        return idxMap[*(T[v].bag.begin())];
    }
    else{
        typename boost::graph_traits<T_t>::vertex_descriptor child = *(boost::adjacent_vertices(v, T).first);
        std::set<unsigned int>::iterator sIt1, sIt2, sEnd1, sEnd2;
        sIt1 = T[v].bag.begin();
        sIt2 = T[child].bag.begin();
        sEnd1 = T[v].bag.end();
        sEnd2 = T[child].bag.end();

        //If *sIt1 != *sIt2, then *sIt1 must be the introduced vertex.
        for(; sIt1 != sEnd1 && sIt2 != sEnd2; ){
            if(*sIt1 == *sIt2){
                sIt1++;
                sIt2++;
            }
            else{
                return idxMap[*sIt1];
            }
        }
        return idxMap[*(T[v].bag.rbegin())];
    }
}

template <typename T_t>
unsigned int get_introduced_vertex_id(typename boost::graph_traits<T_t>::vertex_descriptor v, T_t &T){
    if(T[v].bag.size() == 1){
        return *(T[v].bag.begin());
    }
    else{
        typename boost::graph_traits<T_t>::vertex_descriptor child = *(boost::adjacent_vertices(v, T).first);
        std::set<unsigned int>::iterator sIt1, sIt2, sEnd1, sEnd2;
        sIt1 = T[v].bag.begin();
        sIt2 = T[child].bag.begin();
        sEnd1 = T[v].bag.end();
        sEnd2 = T[child].bag.end();

        //If *sIt1 != *sIt2, then *sIt1 must be the introduced vertex.
        for(; sIt1 != sEnd1 && sIt2 != sEnd2; ){
            if(*sIt1 == *sIt2){
                sIt1++;
                sIt2++;
            }
            else{
                return *sIt1;
            }
        }
        return *(T[v].bag.rbegin());
    }
}

template <typename G_t, typename T_t>
typename boost::graph_traits<G_t>::vertex_descriptor get_forgotten_vertex(typename boost::graph_traits<G_t>::vertex_descriptor v, T_t &T,
          std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &idxMap){
    typename boost::graph_traits<T_t>::vertex_descriptor child = *(boost::adjacent_vertices(v, T).first);

    if(T[child].bag.size() == 1){
        return idxMap[*(T[child].bag.begin())];
    }

    std::set<unsigned int>::iterator sIt1, sIt2, sEnd1, sEnd2;
    sIt1 = T[v].bag.begin();
    sIt2 = T[child].bag.begin();
    sEnd1 = T[v].bag.end();
    sEnd2 = T[child].bag.end();

    //if *sIt1 != *sIt2, then *sIt2 must be the forgotten vertex.
    for(; sIt1 != sEnd1 && sIt2 != sEnd2; ){
        if(*sIt1 == *sIt2){
            sIt1++;
            sIt2++;
        }
        else{
            return idxMap[*sIt2];
        }
    }
    return idxMap[*(T[child].bag.rbegin())];
}


enum enum_node_type { LEAF, INTRODUCE, FORGET, JOIN };

//Returns the type of a node in a nice tree decomposition.
template <typename T_t>
enum_node_type get_type(typename boost::graph_traits<T_t>::vertex_descriptor v, T_t &T){
    if(boost::out_degree(v, T) == 2){
        return JOIN;
    }
    else if(boost::out_degree(v, T) == 1){
        typename boost::graph_traits<T_t>::vertex_descriptor child = *(boost::adjacent_vertices(v, T).first);

        if(T[v].bag.size() > T[child].bag.size()){
            return INTRODUCE;
        }
        else if(T[v].bag.size() < T[child].bag.size()){
            return FORGET;
        }
        else{
        }
    }
    else if(boost::out_degree(v, T) == 0){
        return LEAF;
    }
    else{
    }
}

template <typename T_t>
typename boost::graph_traits<T_t>::vertex_descriptor next_node_preorder(typename boost::graph_traits<T_t>::vertex_descriptor v, T_t &T, std::vector<bool> &visited){
    typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
    if(!visited[v]){
        visited[v] = true;
        return v;
    }
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, T); nIt != nEnd; nIt++){
        if(!visited[*nIt]){
            visited[*nIt] = true;
            return *nIt; 
        }
        typename boost::graph_traits<T_t>::vertex_descriptor rtn = next_node_preorder(*nIt, T, visited);
        if(rtn != -1){
            visited[rtn] = true;
            return rtn;
        }
        else{
            continue;
        }
    }
    return -1;
}

template <typename T_t>
typename boost::graph_traits<T_t>::vertex_descriptor next_node_postorder(typename boost::graph_traits<T_t>::vertex_descriptor v, T_t &T, std::vector<bool> &visited){
    typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, T); nIt != nEnd; nIt++){
        if(!visited[*nIt]){
            return next_node_postorder(*nIt, T, visited);
        }
    }
    visited[v] = true;
    return v;
}


// author: Philipp Klaus Krause, philipp@informatik.uni-frankfurt.de, pkk@spth.de, 2010 - 2011
// functions: nicify_joins, nicify_diffs, nicify_diffs_more, find_root, nicify

// Ensure that all joins are at proper join nodes: Each node that has two children has the same bag as its children.
// Complexity: Linear in the number of vertices of T.
template <class T_t>
void nicify_joins(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor t){
    typename boost::graph_traits<T_t>::adjacency_iterator c, c_end;
    typename boost::graph_traits<T_t>::vertex_descriptor c0, c1;

    boost::tie(c, c_end) = boost::adjacent_vertices(t, T);

    switch (boost::out_degree(t, T)){
        case 0:
            return;
        case 1:
            nicify_joins(T, *c);
            return;
        case 2:
            break;
        default:
            c0 = *c++;
            c1 = *c;
            typename boost::graph_traits<T_t>::vertex_descriptor d = boost::add_vertex(T);
            boost::add_edge(d, c0, T);
            boost::add_edge(d, c1, T);

            boost::remove_edge(t, c0, T);
            boost::remove_edge(t, c1, T);

            T[d].bag = T[t].bag;
            boost::add_edge(t, d, T);

            nicify_joins(T, t);
            return;
    }

    c0 = *c++;
    c1 = *c;

    nicify_joins(T, c0);

    if(T[t].bag != T[c0].bag){
        typename boost::graph_traits<T_t>::vertex_descriptor d = boost::add_vertex(T);
        boost::add_edge(d, c0, T);
        boost::add_edge(t, d, T);
        boost::remove_edge(t, c0, T);
        T[d].bag = T[t].bag;
    }

    nicify_joins(T, c1);

    if(T[t].bag != T[c1].bag){
        typename boost::graph_traits<T_t>::vertex_descriptor d = boost::add_vertex(T);
        boost::add_edge(d, c1, T);
        boost::add_edge(t, d, T);
        boost::remove_edge(t, c1, T);
        T[d].bag = T[t].bag;
    }
}

// Ensure that all nodes' bags are either a subset or a superset of their successors'.
// Complexity: Linear in the number of vertices of T.
template <class T_t>
void nicify_diffs(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor t){
    typename boost::graph_traits<T_t>::adjacency_iterator c, c_end;
    typename boost::graph_traits<T_t>::vertex_descriptor c0, c1;

    boost::tie(c, c_end) = boost::adjacent_vertices(t, T);

    switch(boost::out_degree(t, T)){
        case 0:
/*
            //ensures that leafs have empty bags
            if(T[t].bag.size())
                boost::add_edge(t, boost::add_vertex(T), T);
*/
                return;
        case 1:
            break;
        case 2:
            c0 = *c++;
            c1 = *c;

            nicify_diffs(T, c0);
            nicify_diffs(T, c1);
            return;
        default:
            //an error occured
            return;
    }

    c0 = *c;
    nicify_diffs(T, c0);

    if(std::includes(T[t].bag.begin(), T[t].bag.end(), T[c0].bag.begin(), T[c0].bag.end())
    || std::includes(T[c0].bag.begin(), T[c0].bag.end(), T[t].bag.begin(), T[t].bag.end()))
        return;

    typename boost::graph_traits<T_t>::vertex_descriptor d = boost::add_vertex(T);

    boost::add_edge(d, c0, T);
    boost::add_edge(t, d, T);
    boost::remove_edge(t, c0, T);

    std::set_intersection(T[t].bag.begin(), T[t].bag.end(), T[c0].bag.begin(), T[c0].bag.end(), std::inserter(T[d].bag, T[d].bag.begin()));
}

// // Ensure that all bag sizes of adjacent bags differ by at most one.
template <class T_t>
void nicify_diffs_more(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor t){
    typename boost::graph_traits<T_t>::adjacency_iterator c, c_end;
    typename boost::graph_traits<T_t>::vertex_descriptor c0, c1;

    boost::tie(c, c_end) = boost::adjacent_vertices(t, T);

    switch(boost::out_degree(t, T)){
        case 0:
            if(T[t].bag.size() > 1){
                typename boost::graph_traits<T_t>::vertex_descriptor d = boost::add_vertex(T);
                T[d].bag = T[t].bag;
                T[d].bag.erase(T[d].bag.begin());
                boost::add_edge(t, d, T);

                nicify_diffs_more(T, t);
            }
            return;
        case 1:
            break;
        case 2:
            c0 = *c++;
            c1 = *c;

            nicify_diffs_more(T, c0);
            nicify_diffs_more(T, c1);
            return;
        default:
            //an error occured
            return;
    }

    c0 = *c;

    size_t c0_size, t_size;
    t_size = T[t].bag.size();
    c0_size = T[c0].bag.size();

    if(t_size <= c0_size + 1 && t_size + 1 >= c0_size){
        nicify_diffs_more(T, c0);
        return;
    }

    typename boost::graph_traits<T_t>::vertex_descriptor d = boost::add_vertex(T);
    boost::add_edge(d, c0, T);
    boost::add_edge(t, d, T);
    boost::remove_edge(t, c0, T);

    T[d].bag = T[t_size > c0_size ? t : c0].bag;
    std::set<unsigned int>::iterator i;

    for(i = T[d].bag.begin(); T[t_size < c0_size ? t : c0].bag.find(*i) != T[t_size < c0_size ? t : c0].bag.end(); ++i);

    T[d].bag.erase(i);

    nicify_diffs_more(T, t);
}

// Find a root of an acyclic graph T
// Complexity: Linear in the number of vertices of T.
template <class T_t>
typename boost::graph_traits<T_t>::vertex_descriptor find_root(T_t &T){
    typename boost::graph_traits<T_t>::vertex_descriptor t = *(boost::vertices(T).first);
    typename boost::graph_traits<T_t>::in_edge_iterator e, e_end;

    for(boost::tie(e, e_end) = boost::in_edges(t, T); e != e_end; boost::tie(e, e_end) = boost::in_edges(t, T))
        t = boost::source(*e, T);

    return t;
}

// Transform a tree decomposition into a nice tree decomposition.
template <class T_t>
void nicify(T_t &T){
    typename boost::graph_traits<T_t>::vertex_descriptor t = find_root(T);

    // Ensure we have an empty bag at the root.
    if(T[t].bag.size() > 0){
        typename boost::graph_traits<T_t>::vertex_descriptor d = t;
        t = boost::add_vertex(T);
        boost::add_edge(t, d, T);
    }

    nicify_joins(T, t);
    nicify_diffs(T, t);
    nicify_diffs_more(T, t);
}

} //namespace nice

} //namespace treedec

#endif
