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
// Offers miscellaneous algorithms about computing tree decompositions.
//
//
// These functions are most likely to be interesting for outside use:
//
// bool is_valid_decomposition(G_t const& G, T_t const& T)
// void trivial_decomposition(G_t &G, T_t &T)
// int get_width(T_t &T)
// float get_average_bag_size(T_t &T)
// unsigned int get_adhesion(T_t T)
// void make_small(T_t &T)
// void glue_decompositions(T_t &T1, T_t &T2)
//

#ifndef TD_MISC
#define TD_MISC

// here (for now)
#ifndef unreachable
#define unreachable()
#endif
#ifndef incomplete
#define incomplete()
#endif

#include <stack>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include "graph.hpp"
#include "simple_graph_algos.hpp"
#include "std.hpp"

namespace treedec{

// Find a root of an acyclic graph T.
// Complexity: Linear in the number of vertices of T.
template <typename T_t>
typename boost::graph_traits<T_t>::vertex_descriptor find_root(T_t &T){
    typename boost::graph_traits<T_t>::vertex_descriptor t = *(boost::vertices(T).first);
    typename boost::graph_traits<T_t>::in_edge_iterator e, e_end;
    std::vector<bool> visited(boost::num_vertices(T), false);

    visited[t] = true;

    for(boost::tie(e, e_end)=boost::in_edges(t, T); e!=e_end;
        boost::tie(e, e_end)=boost::in_edges(t, T)){
        if(!visited[t]){
            t = boost::source(*e, T);
            visited[t] = true;
        }
        else{ //T is undirected
            return t;
        }
    }

    return t;
}

template <typename T_t>
void postorder_traversal(T_t &T, std::stack<typename boost::graph_traits<T_t>::vertex_descriptor> &S){
    std::stack<typename boost::graph_traits<T_t>::vertex_descriptor> S_tmp;

    std::vector<bool> visited(boost::num_vertices(T), false);

    //The root can be chosen freely.
    typename boost::graph_traits<T_t>::vertex_descriptor root = treedec::find_root(T);
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

template <typename T_t>
bool validate_connectivity(T_t &T){
    //Compute a postorder traversal.
    std::stack<typename boost::graph_traits<T_t>::vertex_descriptor> S;
    treedec::postorder_traversal(T, S);

    std::vector<bool> visited(boost::num_vertices(T), false);
    typename noboost::treedec_traits<T_t>::bag_type forgotten;

    while(true){
        typename boost::graph_traits<T_t>::vertex_descriptor cur = S.top();
        S.pop();
        visited[cur] = true;

        //Get the parent of cur.
        typename boost::graph_traits<T_t>::vertex_descriptor parent;
        if(!S.empty()){
            typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(cur, T); nIt != nEnd; nIt++){
                if(!visited[*nIt]){
                    parent = *nIt;
                    break;
                }
            }
        }

        //Test if forgotten and noboost::bag(cur, T) have an entry in common.
        typedef typename noboost::treedec_traits<T_t>::bag_type::const_iterator const_iterator;
        const_iterator it1 = forgotten.begin();
        const_iterator it2 = noboost::bag(cur, T).begin();

        for(; it1 != forgotten.end() && it2 != noboost::bag(cur, T).end(); ){
            if(*it1 == *it2){
                //There are coded vertices, that are not connected in T.
                return false;
            }
            else if(*it1 < *it2){
                it1++;
            }
            else{
                it2++;
            }
        }

        if(S.empty()){
            return true;
        }

        std::set_difference(noboost::bag(cur, T).begin(),
                            noboost::bag(cur, T).end(),
                            noboost::bag(parent, T).begin(),
                            noboost::bag(parent, T).end(),
                            std::inserter(forgotten, forgotten.begin()));
    }
}


/* Checks if a tree decomposition is valid with respect to G.
 *
 *  0 = valid
 * -1 = (at least) not a tree
 * -2 = (at least) not all vertices covered (but a tree)
 * -3 = (at least) not all edges covered (but a tree and all vertices covered)
 * -4 = there exist vertices, coded in the bags of T, that are not connected in T
 *                           (but T is a tree and all edges/vertices are covered)
 * -5 = The empty graph must have a treedecomposition with one vertex and an empty bag
 */
template <typename G_t, typename T_t>
int is_valid_treedecomposition(G_t const& G, T_t const& T){
    if(boost::num_vertices(T) == 0){
        //The empty graph has a treedecomposition with 1 vertex and an empty bag.
        return -5;
    }

    //Checks if T is a tree.
    std::vector<int> component(boost::num_vertices(T));
    int num = boost::connected_components(T, &component[0]);
    if(num > 1 || boost::num_edges(T) > boost::num_vertices(T)-1){
        return -1;
    }

    //Checks if exactly the vertices of G are covered.
    typename noboost::treedec_traits<T_t>::bag_type coded_vertices;
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        coded_vertices.insert(noboost::bag(*tIt, T).begin(),
                              noboost::bag(*tIt, T).end());
    }

    typename noboost::treedec_traits<T_t>::bag_type vertices;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        typename noboost::treedec_traits<T_t>::bag_type::value_type v;
        v = *vIt;
        vertices.insert(v);
    }

    if(coded_vertices != vertices){
        return -2;
    }

    //Checks if all edges are covered.
    typename std::vector<typename noboost::treedec_traits<T_t>::bag_type> edges(boost::num_edges(G));
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
            if(*vIt >= *nIt){
                continue;
            }

            bool is_contained = false;
            for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
                typedef typename noboost::treedec_traits<T_t>::bag_type::value_type vd_type;
                vd_type v(*vIt);
                vd_type n(*nIt);

                if(noboost::bag(*tIt, T).find(v) != noboost::bag(*tIt, T).end() &&
                   noboost::bag(*tIt, T).find(n) != noboost::bag(*tIt, T).end()){
                    is_contained = true;
                    break;
                }
            }

            if(!is_contained){
                return -3; //Not all edges are covered.
            }
        }
    }

    if(!validate_connectivity(T)){
        return -4;
    }

    return 0;
}

template <typename G_t, typename T_t>
void trivial_decomposition(G_t &G, T_t &T){
    typename boost::graph_traits<T_t>::vertex_descriptor t = boost::add_vertex(T);

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        typename noboost::treedec_traits<T_t>::bag_type::value_type v = *vIt;
        noboost::bag(t, T).insert(v);
    }
}

template <typename T_t>
int get_width(T_t &T){
    int max = -1;
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        size_t bag_size = noboost::bag(*tIt, T).size();
        if((int)bag_size > max){
            max = (int)bag_size;
        }
    }

    return (max-1);
}

template <typename T_t>
float get_average_bag_size(T_t &T){
    float avg = 0.0;
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        avg += noboost::bag(*tIt, T).size();
    }

    return (boost::num_vertices(T) > 0)? avg/boost::num_vertices(T) : 0.0;
}

template <typename T_t>
void make_small(T_t &T){
    while(true){
        bool modified = false;
        std::vector<typename boost::graph_traits<T_t>::vertex_descriptor > N;
        typename boost::graph_traits<T_t>::vertex_descriptor child, parent;

        typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*tIt, T); nIt != nEnd; nIt++){
                if(*tIt == *nIt){
                    continue;
                }
                if(std::includes(noboost::bag(*nIt, T).begin(), noboost::bag(*nIt, T).end(),
                                 noboost::bag(*tIt, T).begin(), noboost::bag(*tIt, T).end())){
                    child = *tIt;
                    parent = *nIt;

                    N.resize(boost::degree(*tIt, T)-1);
                    unsigned int c = 0;
                    typename boost::graph_traits<T_t>::adjacency_iterator nIt2, nEnd2;
                    for(boost::tie(nIt2, nEnd2) = boost::adjacent_vertices(*tIt, T); nIt2 != nEnd2; nIt2++){
                        if(*nIt2 != parent){
                            N[c++] = *nIt2;
                        }
                    }

                    modified = true;
                    break;
                }
            }
            if(modified){
                for(unsigned int i = 0; i < N.size(); i++){
                    boost::add_edge(parent, N[i], T);
                }

                boost::clear_vertex(child, T);
                boost::remove_vertex(child, T);
                break;
            }
        }
        if(!modified){
            return;
        }
    }
}

//Glues a single bag with the current tree decomposition.
template <class bagtype, typename T_t>
void glue_bag(bagtype &bag, typename bagtype::value_type elim_vertex, T_t &T){
    typename boost::graph_traits<T_t>::vertex_descriptor t_dec_node;

    typename boost::graph_traits<T_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(T); vIt != vEnd; vIt++){
        if(std::includes(noboost::bag(*vIt, T).begin(),
                         noboost::bag(*vIt, T).end(),
                         bag.begin(), bag.end()))
        {
            if(noboost::bag(*vIt, T).find(elim_vertex) != noboost::bag(*vIt, T).end()){
                return;
            }

            t_dec_node = boost::add_vertex(T);
            noboost::bag(t_dec_node, T) = MOVE(bag);
            bag.clear();
            noboost::bag(t_dec_node, T).insert(elim_vertex);
            boost::add_edge(*vIt, t_dec_node, T);
            return;
        }
    }

    if(boost::num_vertices(T) > 0){
        boost::tie(vIt, vEnd) = boost::vertices(T);
    }

    t_dec_node = boost::add_vertex(T);
    noboost::bag(t_dec_node, T) = MOVE(bag);
    bag.clear();
    noboost::bag(t_dec_node, T).insert(elim_vertex);

    if(boost::num_vertices(T) > 1){
        boost::add_edge(*vIt, t_dec_node, T);
    }
}

//glues two "disjoint" decompositions (e.g. decompositions of two components of a graph)
template <typename T_t>
void glue_decompositions(T_t &T1, T_t &T2){
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;

    //Copy T2 to T1 and add an edge from root to an arbitrary vertex of T2.
    std::vector<typename boost::graph_traits<T_t>::vertex_descriptor> idxMap(boost::num_vertices(T2));
    std::map<typename boost::graph_traits<T_t>::vertex_descriptor, unsigned int> vertex_map;
    unsigned int id = 0;
    for(boost::tie(tIt, tEnd) = boost::vertices(T2); tIt != tEnd; tIt++){
        idxMap[id] = boost::add_vertex(T1);
        vertex_map.insert(std::pair<typename boost::graph_traits<T_t>::vertex_descriptor, unsigned int>(*tIt, id));
        noboost::bag(idxMap[id++], T1) = noboost::bag(*tIt, T2);
    }

    typename boost::graph_traits<T_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(T2); eIt != eEnd; eIt++){
        typename std::map<typename boost::graph_traits<T_t>::vertex_descriptor, unsigned int>::iterator v, w;
        v = vertex_map.find(boost::source(*eIt, T2));
        w = vertex_map.find(boost::target(*eIt, T2));

        boost::add_edge(idxMap[v->second], idxMap[w->second], T1);
    }

    typename boost::graph_traits<T_t>::vertex_iterator tIt2, tEnd2;

    boost::tie(tIt, tEnd) = boost::vertices(T1);

    boost::add_edge(*tIt, idxMap[0], T1);
}


template <typename G_t>
void map_descriptors(std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &S,
                     std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &S_,
                     G_t &H,
                     typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &vdMap)
{
    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt
                                 = S.begin(); sIt != S.end(); sIt++)
    {
        unsigned int pos = noboost::get_pos(*sIt, H);
        S_.insert(vdMap[pos]);
    }
}

template <typename G_t>
void map_descriptors_to_bags(std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &S,
                             typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &B)
{
    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt =
                  S.begin(); sIt != S.end(); sIt++)
    {
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type::value_type vd = *sIt;
        B.insert(vd);
    }
}

template <typename G_t, typename T_t>
void apply_map_on_treedec(T_t &T, G_t &G, typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &vdMap){
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        typename noboost::treedec_traits<T_t>::bag_type bag_old, bag_new;
        bag_old = noboost::bag(*tIt, T);
        for(typename noboost::treedec_traits<T_t>::bag_type::iterator sIt = bag_old.begin(); sIt != bag_old.end(); sIt++){
            unsigned int pos = noboost::get_pos(*sIt, G);
            bag_new.insert(vdMap[pos]);
        }
        noboost::bag(*tIt, T) = MOVE(bag_new);
    }
}


#ifndef TD_SUBSETS
#define TD_SUBSETS

//Collects all subsets of 'X' of size 'k' and stores it in 'subs'.
template <typename T>
void subsets(std::set<T> &X, int size, int k, int idx, std::vector<T> &sub, std::vector<std::set<T> > &subs){
    if(k==0){
        typename std::set<T> subS;
        for(unsigned int i = 0; i < sub.size(); i++){
            subS.insert(sub[i]);
        }
        subs.push_back(subS);
        return;
    }

    int i = idx;
    typename std::set<T>::iterator sIt = X.begin();
    std::advance(sIt, i);
    for(; i<size;i++){
        sub.push_back(*sIt);
        subsets(X,size,k-1,i+1,sub, subs);
        sub.pop_back();
        sIt++;
    }
}

#endif //TD_SUBSETS

} // namespace treedec

namespace noboost { // MOVE
template<class G, class CB>
void make_clique_and_hijack(
        typename boost::graph_traits<G>::vertex_descriptor c,
        G& g, CB* cb, typename outedge_set<G>::type& bag)
{ itested();
    if(cb){incomplete();
        // probably unneeded now.
    }
    noboost::make_clique(boost::adjacent_vertices(c, g), g);
    return hijack_neighborhood(c, g, bag);
}

} // noboost

#endif //TD_MISC

// vim:ts=8:sw=4:et
