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
// bool is_valid_decomposition(G_t &G, T_t T)
// void trivial_decomposition(G_t &G, T_t &T)
// int get_width(T_t &T)
// float get_average_bag_size(T_t &T)
// unsigned int get_adhesion(T_t T)
// void make_small(T_t &T)
// void glue_decompositions(T_t &T1, T_t &T2)
//

#ifndef TD_MISC
#define TD_MISC

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include "TD_simple_graph_algos.hpp"

namespace treedec{


/* checks if a tree decomposition is valid with respect to G
 *
 *  0 = valid
 * -1 = (at least) not a tree
 * -2 = (at least) not all vertices covered (but a tree)
 * -3 = (at least) not all edges covered (but a tree and all vertices covered)
 * -4 = there exist vertices, coded in the bags of T, that are not connected in T
 *                           (but T is a tree and all edges/vertices are covered)
 */
template <typename G_t, typename T_t>
int is_valid_treedecomposition(G_t G, T_t T){
    //checks if T is a tree
    std::vector<int> component(boost::num_vertices(T));
    int num = boost::connected_components(T, &component[0]);
    if(num > 1 || boost::num_edges(T) > boost::num_vertices(T)-1)
        return -1;

    //checks if exactly the vertices of G are covered
    std::set<unsigned int> coded_vertices;
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++)
        coded_vertices.insert(T[*tIt].bag.begin(), T[*tIt].bag.end());

    std::set<unsigned int> vertices;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        vertices.insert(G[*vIt].id);

    if(coded_vertices != vertices)
        return -2;

    //checks if all edges are covered
    std::vector<std::set<unsigned int> > edges;

    for (boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        std::set<unsigned int> edge;
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
            edge.insert(G[*vIt].id);
            edge.insert(G[*nIt].id);
            edges.push_back(edge);
            edge.clear();
        }
    }
    for(std::vector<std::set<unsigned int> >::iterator it = edges.begin(); it != edges.end(); it++){
        bool isSubset = false;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            if(std::includes(T[*tIt].bag.begin(), T[*tIt].bag.end(), it->begin(), it->end())){
                isSubset = true;
                break;
            }
        }
        if(!isSubset)
            //not all edges covered
            return -3;
    }
    std::set<unsigned int> forgotten;

    while(true){
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            unsigned int degree = boost::out_degree(*tIt, T);
            if(degree < 2){
                std::set<unsigned int> intersection;
                std::set_intersection(forgotten.begin(), forgotten.end(), T[*tIt].bag.begin(), T[*tIt].bag.end(), std::inserter(intersection, intersection.begin()));
                if(!intersection.empty()){
                    //there are coded vertices, that are not connected in T
                    return -4;
                }

                if(degree == 1){
                    typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                    typename boost::graph_traits<T_t>::vertex_descriptor parent;
                    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*tIt, T); nIt != nEnd; nIt++)
                        parent = *nIt;

                    std::set_difference(T[*tIt].bag.begin(), T[*tIt].bag.end(), T[parent].bag.begin(), T[parent].bag.end(), std::inserter(forgotten, forgotten.begin()));
                }

                boost::clear_vertex(*tIt, T);
                boost::remove_vertex(*tIt, T);
                break;
            }
        }

        if(boost::num_vertices(T) == 0)
            return 0;
    }
}

template <typename G_t, typename T_t>
void trivial_decomposition(G_t &G, T_t &T){
    std::set<unsigned int> bag;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        bag.insert(G[*vIt].id);

    typename boost::graph_traits<T_t>::vertex_descriptor t;
    t = boost::add_vertex(T);
    T[t].bag = bag;
}

template <typename T_t>
int get_width(T_t &T){
    int max = -1;
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++)
        max = ((int)T[*tIt].bag.size() > max)? (int)T[*tIt].bag.size() : max;

    return (max-1);
}

template <typename T_t>
float get_average_bag_size(T_t &T){
    float avg = 0.0;
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++)
        avg += T[*tIt].bag.size();

    return avg/boost::num_vertices(T);
}

template <typename T_t>
unsigned int get_adhesion(T_t T){
    unsigned int max = 0;
    while(boost::num_vertices(T) > 1){
        typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            if(T[*tIt].bag.size() == 0 || boost::out_degree(*tIt, T) == 0){
                boost::clear_vertex(*tIt, T);
                boost::remove_vertex(*tIt, T);
                continue;
            }
            if(boost::out_degree(*tIt, T) == 1){
                typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                typename boost::graph_traits<T_t>::vertex_descriptor parent;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*tIt, T); nIt != nEnd; nIt++)
                    parent = *nIt;

                std::set<unsigned int> intersection;
                std::set_intersection(T[*tIt].bag.begin(), T[*tIt].bag.end(), T[parent].bag.begin(), T[parent].bag.end(), std::inserter(intersection, intersection.begin()));

                max = (intersection.size() > max)? intersection.size() : max;

                boost::clear_vertex(*tIt, T);
                boost::remove_vertex(*tIt, T);
                break;
            }
        }
    }
    return max;
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
                if(*tIt == *nIt)
                    continue;
                if(std::includes(T[*nIt].bag.begin(), T[*nIt].bag.end(), T[*tIt].bag.begin(), T[*tIt].bag.end())){
                    child = *tIt;
                    parent = *nIt;

                    typename boost::graph_traits<T_t>::adjacency_iterator nIt2, nEnd2;
                    for(boost::tie(nIt2, nEnd2) = boost::adjacent_vertices(*tIt, T); nIt2 != nEnd2; nIt2++){
                        if(*nIt2 != parent)
                            N.push_back(*nIt2);
                    }

                    modified = true;
                    break;
                }
            }
            if(modified){
                for(unsigned int i = 0; i < N.size(); i++)
                    boost::add_edge(parent, N[i], T);

                boost::clear_vertex(child, T);
                boost::remove_vertex(child, T);
                break;
            }
        }
        if(!modified)
            return;
    }
}

//glues a single bag with the current tree decomposition
template <typename T_t>
void glue_bag(std::set<unsigned int> &bag, unsigned int elim_vertex, T_t &T){
    typename boost::graph_traits<T_t>::vertex_descriptor t_dec_node;

    typename boost::graph_traits<T_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(T); vIt != vEnd; vIt++){
        if(std::includes(T[*vIt].bag.begin(), T[*vIt].bag.end(), bag.begin(), bag.end())){
            if(T[*vIt].bag.find(elim_vertex) != T[*vIt].bag.end())
                return;

            t_dec_node = boost::add_vertex(T);
            bag.insert(elim_vertex);
            T[t_dec_node].bag = bag;
            boost::add_edge(*vIt, t_dec_node, T);
            bag.clear();
            return;
        }
    }

    if(boost::num_vertices(T) > 0)
        boost::tie(vIt, vEnd) = boost::vertices(T);

    t_dec_node = boost::add_vertex(T);
    bag.insert(elim_vertex);
    T[t_dec_node].bag = bag;
    bag.clear();

    if(boost::num_vertices(T) > 1)
        boost::add_edge(*vIt, t_dec_node, T);
}

//glues two "disjoint" decompositions (e.g. decompositions of two components of a graph)
template <typename T_t>
void glue_decompositions(T_t &T1, T_t &T2){
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;

    //copy T2 to T1 and add an edge from root to an arbitrary vertex of T2
    std::vector<typename boost::graph_traits<T_t>::vertex_descriptor> idxMap(boost::num_vertices(T2));
    std::map<typename boost::graph_traits<T_t>::vertex_descriptor, unsigned int> vertex_map;
    unsigned int id = 0;
    for(boost::tie(tIt, tEnd) = boost::vertices(T2); tIt != tEnd; tIt++){
        idxMap[id] = boost::add_vertex(T1);
        vertex_map.insert(std::pair<typename boost::graph_traits<T_t>::vertex_descriptor, unsigned int>(*tIt, id));
        T1[idxMap[id++]].bag = T2[*tIt].bag;
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
void reorder_ids_graph(G_t &G, std::vector<unsigned int> &id_map){
    id_map.resize(boost::num_vertices(G));
    unsigned int k = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        id_map[k] = G[*vIt].id;
        G[*vIt].id = k++;
    }
}

template <typename G_t>
void undo_reorder_ids_graph(G_t &G, std::vector<unsigned int> &id_map){
    unsigned int k = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        G[*vIt].id = id_map[k++];
}

template <typename G_t>
void descriptor_bag_to_id_bag(G_t &G, std::set<unsigned int> &id_bag,
        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &desc_bag)
{
    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt =
            desc_bag.begin(); sIt != desc_bag.end(); sIt++){
        unsigned id = noboost::get_id(G,*sIt);
        id_bag.insert(id);
    }
}

template <typename G_t>
void id_bag_to_descriptor_bag(G_t &G, typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &desc_bag, std::set<unsigned int> &id_bag){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(id_bag.find(G[*vIt].id) != id_bag.end())
            desc_bag.insert(*vIt);
    }
}

template <typename T_t>
void reorder_ids_decomposition(T_t &T, std::vector<unsigned int> &id_map){
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        std::set<unsigned int> new_bag;
        for(std::set<unsigned int>::iterator sIt = T[*tIt].bag.begin(); sIt != T[*tIt].bag.end(); sIt++)
            new_bag.insert(id_map[*sIt]);

        T[*tIt].bag = new_bag;
    }
}

template <typename G_t>
void remove_isolated_vertices(G_t &H, G_t &G){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    unsigned int max = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        max = (G[*vIt].id > max)? G[*vIt].id : max;

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap(max+1);
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) > 0){
            idxMap[G[*vIt].id] = boost::add_vertex(H);
            H[idxMap[G[*vIt].id]].id = G[*vIt].id;
        }
    }
    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++)
        boost::add_edge(idxMap[G[boost::source(*eIt, G)].id], idxMap[G[boost::target(*eIt, G)].id], H);
}

#ifndef TD_SUBSETS
#define TD_SUBSETS

//collects all subsets of X of size k in subs (todo: replacement by enumeration in hunt())
static void subsets(std::set<unsigned int> &X, int size, int k, unsigned int idx, std::vector<unsigned int> &sub, std::vector<std::set<unsigned int> > &subs){
    if(k==0){
        std::set<unsigned int> subS;
        for(unsigned int i = 0; i < sub.size(); i++)
            subS.insert(sub[i]);
        subs.push_back(subS);
        return;
    }

    unsigned int i = idx;
    std::set<unsigned int>::iterator sIt = X.begin();
    std::advance(sIt, i);
    for(; i<X.size();i++){
        sub.push_back(*sIt);
        subsets(X,X.size(),k-1,i+1,sub, subs);
        sub.pop_back();
        sIt++;
    }
}

#endif

#ifndef TD_POWERSET
#define TD_POWERSET

inline void powerset(std::set<unsigned int> &X, std::vector<std::set<unsigned int> > &subs){
    std::vector<unsigned int> sub;
    for(unsigned int i = 0; i <=X.size(); i++){
        subsets(X, X.size(), i, 0, sub, subs);
    }
}

#endif

}

#endif
// vim:ts=8:sw=4:et
