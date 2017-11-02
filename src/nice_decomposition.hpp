// Lukas Larisch, 2014 - 2016
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
 * Nice tree decomposition-related stuff.
 *
 * These functions are most likely to be interesting for outside use:
 *
 * - void nicify(T_t &T)
 *
 */

#ifndef TD_NICE_DECOMPOSITION
#define TD_NICE_DECOMPOSITION

#include <set>
#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "graph.hpp"
#include "misc.hpp"
#include "algo.hpp"

#define BAG_(x, y) boost::get(treedec::bag_t(), y, x)

namespace treedec{

namespace nice{

template <typename T_t>
typename treedec_traits<T_t>::bag_type::value_type
          get_introduced_vertex(
              typename boost::graph_traits<T_t>::vertex_descriptor v, T_t &T)
{
    if(BAG_(v, T).size() == 1){
        return *(BAG_(v, T).begin());
    }
    else{
        typename boost::graph_traits<T_t>::vertex_descriptor child =
                                  *(boost::adjacent_vertices(v, T).first);
        typename treedec_traits<T_t>::bag_type::iterator sIt1, sIt2, sEnd1, sEnd2;
        sIt1 = BAG_(v, T).begin();
        sIt2 = BAG_(child, T).begin();
        sEnd1 = BAG_(v, T).end();
        sEnd2 = BAG_(child, T).end();

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
        return *(BAG_(v, T).rbegin());
    }
}

template <typename T_t>
typename treedec_traits<T_t>::bag_type::value_type
          get_forgotten_vertex(
              typename boost::graph_traits<T_t>::vertex_descriptor v, T_t &T)
{
    typename boost::graph_traits<T_t>::vertex_descriptor child =
                            *(boost::adjacent_vertices(v, T).first);

    if(BAG_(child, T).size() == 1){
        return *(BAG_(child, T).begin());
    }

    typename treedec_traits<T_t>::bag_type::iterator sIt1, sIt2, sEnd1, sEnd2;
    sIt1 = BAG_(v, T).begin();
    sIt2 = BAG_(child, T).begin();
    sEnd1 = BAG_(v, T).end();
    sEnd2 = BAG_(child, T).end();

    //If *sIt1 != *sIt2, then *sIt2 must be the forgotten vertex.
    for(; sIt1 != sEnd1 && sIt2 != sEnd2; ){
        if(*sIt1 == *sIt2){
            sIt1++;
            sIt2++;
        }
        else{
            return *sIt2;
        }
    }
    return *(BAG_(child, T).rbegin());
}


enum enum_node_type { LEAF, INTRODUCE, FORGET, JOIN, INVALID };

//Returns the type of a node in a nice tree decomposition.
template <typename T_t>
enum_node_type get_type(typename boost::graph_traits<T_t>::vertex_descriptor v, T_t &T){
    if(boost::out_degree(v, T) == 2){
        return JOIN;
    }
    else if(boost::out_degree(v, T) == 1){
        typename boost::graph_traits<T_t>::vertex_descriptor child = *(boost::adjacent_vertices(v, T).first);

        if(BAG_(v, T).size() > BAG_(child, T).size()){
            return INTRODUCE;
        }
        else if(BAG_(v, T).size() < BAG_(child, T).size()){
            return FORGET;
        }
        else{
            return INVALID;
        }
    }
    else if(boost::out_degree(v, T) == 0){
        return LEAF;
    }
    else{
        return INVALID;
    }
}

//Find a root of an acyclic graph T
//Complexity: Linear in the number of vertices of T.
template <class T_t>
typename boost::graph_traits<T_t>::vertex_descriptor find_root(T_t &T){
    typename boost::graph_traits<T_t>::vertex_descriptor t = *(boost::vertices(T).first);
    typename boost::graph_traits<T_t>::in_edge_iterator e, e_end;

    for(boost::tie(e, e_end) = boost::in_edges(t, T); e != e_end; boost::tie(e, e_end) = boost::in_edges(t, T)){
        t = boost::source(*e, T);
    }

    return t;
}

template <typename T_t>
void postorder_traversal(T_t &T, std::stack<typename boost::graph_traits<T_t>::vertex_descriptor> &S){
    std::stack<typename boost::graph_traits<T_t>::vertex_descriptor> S_tmp;

    std::vector<BOOL> visited(boost::num_vertices(T), false);

    //The root can be chosen freely.
    typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::nice::find_root(T);
    S_tmp.push(root);
    visited[root] = true;

    while(!S_tmp.empty()){
        typename boost::graph_traits<T_t>::vertex_descriptor v = S_tmp.top();
        S_tmp.pop();
        S.push(v);

        typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, T); nIt != nEnd; nIt++){
            if(!visited[*nIt]){
                S_tmp.push(*nIt);
                visited[*nIt] = true;
            }
        }
    }
}

//author: Philipp Klaus Krause, philipp@informatik.uni-frankfurt.de, pkk@spth.de, 2010 - 2011
//functions: nicify_joins, nicify_diffs, nicify_diffs_more, find_root, nicify

//Ensure that all joins are at proper join nodes: Each node that has two children has the same bag as its children.
//Complexity: Linear in the number of vertices of T.
template <class T_t>
void nicify_joins(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor t){
    typename boost::graph_traits<T_t>::adjacency_iterator c, c_end;
    typename boost::graph_traits<T_t>::vertex_descriptor c0, c1;

    boost::tie(c, c_end) = boost::adjacent_vertices(t, T);

    switch (boost::out_degree(t, T)){
        case 0:
            return;
        case 1:
            treedec::nice::nicify_joins(T, *c);
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

            boost::get(treedec::bag_t(), T, d) = boost::get(treedec::bag_t(), T, t);
            boost::add_edge(t, d, T);

            treedec::nice::nicify_joins(T, t);
            return;
    }

    c0 = *c++;
    c1 = *c;

    treedec::nice::nicify_joins(T, c0);

    if(BAG_(t, T) != BAG_(c0, T)){
        typename boost::graph_traits<T_t>::vertex_descriptor d = boost::add_vertex(T);
        boost::add_edge(d, c0, T);
        boost::add_edge(t, d, T);
        boost::remove_edge(t, c0, T);
        BAG_(d, T) = BAG_(t, T);
    }

    treedec::nice::nicify_joins(T, c1);

    if(BAG_(t, T) != BAG_(c1, T)){
        typename boost::graph_traits<T_t>::vertex_descriptor d = boost::add_vertex(T);
        boost::add_edge(d, c1, T);
        boost::add_edge(t, d, T);
        boost::remove_edge(t, c1, T);
        BAG_(d, T) = BAG_(t, T);
    }
}

//Ensure that all nodes' bags are either a subset or a superset of their successors'.
//Complexity: Linear in the number of vertices of T.
template <class T_t>
void nicify_diffs(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor t, bool empty_leafs, bool cleanup){
    typename boost::graph_traits<T_t>::adjacency_iterator c, c_end;
    typename boost::graph_traits<T_t>::vertex_descriptor c0, c1;

    boost::tie(c, c_end) = boost::adjacent_vertices(t, T);

    switch(boost::out_degree(t, T)){
        case 0:
            if(empty_leafs){
            //Ensures that leafs have empty bags.
                if(BAG_(t, T).size())
                    boost::add_edge(t, boost::add_vertex(T), T);
            }
            return;
        case 1:
            break;
        case 2:
            c0 = *c++;
            c1 = *c;

            treedec::nice::nicify_diffs(T, c0, empty_leafs, cleanup);
            treedec::nice::nicify_diffs(T, c1, empty_leafs, cleanup);
            return;
        default:
            //An error occured.
            return;
    }

    c0 = *c;
    treedec::nice::nicify_diffs(T, c0, empty_leafs, cleanup);

    if(cleanup){
        // Redundant bags are isolated, and thus marked for later removal.
        if (T[t].bag == T[c0].bag){
            T[c0].bag.clear();
            boost::remove_edge(t, c0, T);
            typename boost::graph_traits<T_t>::adjacency_iterator c, c_end;
            for(boost::tie(c, c_end) = adjacent_vertices(c0, T); c != c_end; ++c){
                boost::add_edge(t, *c, T);
                boost::remove_edge(c0, *c, T);
            }
        }
    }

    if(std::includes(BAG_(t, T).begin(), BAG_(t, T).end(),
                     BAG_(c0, T).begin(), BAG_(c0, T).end())
    || std::includes(BAG_(c0, T).begin(), BAG_(c0, T).end(),
                     BAG_(t, T).begin(), BAG_(t, T).end()))
    {
        return;
    }

    typename boost::graph_traits<T_t>::vertex_descriptor d = boost::add_vertex(T);

    boost::add_edge(d, c0, T);
    boost::add_edge(t, d, T);
    boost::remove_edge(t, c0, T);

    std::set_intersection(BAG_(t, T).begin(), BAG_(t, T).end(),
                          BAG_(c0, T).begin(), BAG_(c0, T).end(),
                          std::inserter(BAG_(d, T), BAG_(d, T).begin()));
}

//Ensure that all bag sizes of adjacent bags differ by at most one.
template <class T_t>
void nicify_diffs_more(T_t &T, typename boost::graph_traits<T_t>::vertex_descriptor t){
    typename boost::graph_traits<T_t>::adjacency_iterator c, c_end;
    typename boost::graph_traits<T_t>::vertex_descriptor c0, c1;

    boost::tie(c, c_end) = boost::adjacent_vertices(t, T);

    switch(boost::out_degree(t, T)){
        case 0:
            if(BAG_(t, T).size() > 1){
                typename boost::graph_traits<T_t>::vertex_descriptor d = boost::add_vertex(T);
                BAG_(d, T) = BAG_(t, T);
                BAG_(d, T).erase(BAG_(d, T).begin());
                boost::add_edge(t, d, T);

                treedec::nice::nicify_diffs_more(T, t);
            }
            return;
        case 1:
            break;
        case 2:
            c0 = *c++;
            c1 = *c;

            treedec::nice::nicify_diffs_more(T, c0);
            treedec::nice::nicify_diffs_more(T, c1);
            return;
        default:
            //an error occured
            return;
    }

    c0 = *c;

    size_t c0_size, t_size;
    t_size = BAG_(t, T).size();
    c0_size = BAG_(c0, T).size();

    if(t_size <= c0_size + 1 && t_size + 1 >= c0_size){
        treedec::nice::nicify_diffs_more(T, c0);
        return;
    }

    typename boost::graph_traits<T_t>::vertex_descriptor d = boost::add_vertex(T);
    boost::add_edge(d, c0, T);
    boost::add_edge(t, d, T);
    boost::remove_edge(t, c0, T);

    BAG_(d, T) = BAG_(t_size > c0_size ? t : c0, T);
    std::set<unsigned int>::iterator i;

    for(i = BAG_(d, T).begin();
        BAG_(t_size < c0_size ? t : c0, T).find(*i)
     != BAG_(t_size < c0_size ? t : c0, T).end(); ++i);

    BAG_(d, T).erase(i);

    treedec::nice::nicify_diffs_more(T, t);
}

//Transform a tree decomposition into a nice tree decomposition.
template <class T_t>
void nicify(T_t &T, bool empty_leafs=false, bool cleanup=false){ //TODO: test empty_leafs=true and cleanup=true
    typename boost::graph_traits<T_t>::vertex_descriptor t = treedec::nice::find_root(T);

    //Ensure we have an empty bag at the root.
    if(boost::get(treedec::bag_t(), T, t).size() > 0){
        typename boost::graph_traits<T_t>::vertex_descriptor d = t;
        t = boost::add_vertex(T);
        boost::add_edge(t, d, T);
    }

    treedec::nice::nicify_joins(T, t);
    treedec::nice::nicify_diffs(T, t, empty_leafs, cleanup);
    treedec::nice::nicify_diffs_more(T, t);

    if(cleanup){
        treedec::remove_isolated_vertices(T);
    }

    assert(boost::num_edges(T)+1 == boost::num_vertices(T));
}

//TODO: cleanup names
template<class T, template<class G_, class ...> class CFGT=treedec::algo::default_config>
class I_nicify{
public: // construct
    I_nicify(T &t, bool empty_leafs, bool cleanup=true) : _t(t), _empty_leafs(empty_leafs), _cleanup(cleanup) { //TODO should cleanup be always true?!
    }

public: // algo interface
    void do_it(){
        nicify(_t, _empty_leafs, _cleanup);
    }

private:
    T& _t;
    bool _empty_leafs;
    bool _cleanup;
}; // nicify


} //namespace nice

} //namespace treedec

#undef BAG_

#endif //TD_NICE_DECOMPOSITION

// vim:ts=8:sw=4:et


