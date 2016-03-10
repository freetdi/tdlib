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

#include <stack>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include "TD_simple_graph_algos.hpp"
#include "TD_noboost.hpp"
#include "TD_std.hpp"

#include <iostream>

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

/*
/ An iterative function to do post order traversal of a given binary tree
void postOrderIterative(struct Node* root)
{
    // Create two stacks
    struct Stack* s1 = createStack(MAX_SIZE);
    struct Stack* s2 = createStack(MAX_SIZE);
 
    // push root to first stack
    push(s1, root);
    struct Node* node;
 
    // Run while first stack is not empty
    while (!isEmpty(s1))
    {
        // Pop an item from s1 and push it to s2
        node = pop(s1);
        push(s2, node);
 
        // Push left and right children of removed item to s1
        if (node->left)
            push(s1, node->left);
        if (node->right)
            push(s1, node->right);
    }
 
    // Print all elements of second stack
    while (!isEmpty(s2))
    {
        node = pop(s2);
        printf("%d ", node->data);
    }
}
*/

template <typename T_t>
bool validate_connectivity(T_t &T, typename noboost::treedec_traits<T_t>::bag_type &forgotten){
    //Compute a postorder traversal.
    std::stack<typename boost::graph_traits<T_t>::vertex_descriptor> s1;
    std::stack<typename boost::graph_traits<T_t>::vertex_descriptor> s2;

    std::vector<bool> visited(boost::num_vertices(T), false);

    //The root can be chosen freely.
    typename boost::graph_traits<T_t>::vertex_descriptor root = *(boost::vertices(T).first);
    s1.push(root);
    visited[noboost::pos(root, T)] = true;

    while(!s1.empty()){
        typename boost::graph_traits<T_t>::vertex_descriptor v = s1.top();
        s1.pop();
        s2.push(v);

        typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, T); nIt != nEnd; nIt++){
            if(!visited[noboost::get_pos(*nIt, T)]){
                s1.push(*nIt);
                visited[noboost::pos(*nIt, T)] = true;
            }
        }
    }

    while(!s2.empty()){
        typename boost::graph_traits<T_t>::vertex_descriptor cur = s2.top();
        s2.pop();
        typename boost::graph_traits<T_t>::vertex_descriptor parent;
        if(!s2.empty()){ parent = s2.top(); }

        typename noboost::treedec_traits<T_t>::bag_type::iterator it1 = forgotten.begin();
        typename noboost::treedec_traits<T_t>::bag_type::iterator it2 = noboost::bag(T, cur).begin();

        //Test if forgotten and noboost::bag(T, cur) have an entry in common.
        for(; it1 != forgotten.end() && it2 != noboost::bag(T, cur).end(); ){
            if(*it1 == *it2){
                //There are coded vertices, that are not connected in T.
                return false;
            }
            else if(*it1 < *it2){ it1++; }
            else{ it2++; }
        }

        if(s2.empty()){
            return true;
        }

        std::set_difference(noboost::bag(T, cur).begin(),
                            noboost::bag(T, cur).end(),
                            noboost::bag(T, parent).begin(),
                            noboost::bag(T, parent).end(),
                            std::inserter(forgotten, forgotten.begin()));

    }
}

template <typename G_t, typename T_t>
int is_valid_treedecomposition(G_t G, T_t T){
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
        coded_vertices.insert(noboost::bag(T, *tIt).begin(),
                              noboost::bag(T, *tIt).end());
    }

    typename noboost::treedec_traits<T_t>::bag_type vertices;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        vertices.insert((typename noboost::treedec_traits<T_t>::bag_type::value_type) *vIt);
    }

    if(coded_vertices != vertices){
        return -2;
    }

    //Checks if all edges are covered.
    typename std::vector<typename noboost::treedec_traits<T_t>::bag_type> edges(boost::num_edges(G));
    unsigned int i = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
            if(*vIt >= *nIt){
                continue;
            }

            bool is_contained = false;
            for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
                if(noboost::bag(T,*tIt).find((typename noboost::treedec_traits<T_t>::bag_type::value_type) *vIt)
                != noboost::bag(T,*tIt).end() &&
                   noboost::bag(T,*tIt).find((typename noboost::treedec_traits<T_t>::bag_type::value_type) *nIt)
                != noboost::bag(T,*tIt).end()){
                    is_contained = true;
                    break;
                }
            }

            if(!is_contained){
                return -3; //Not all edges are covered.
            }
        }
    }

    typename noboost::treedec_traits<T_t>::bag_type forgotten;
    if(!validate_connectivity(T, forgotten)){
        return -4;
    }

    return 0;
}

template <typename G_t, typename T_t>
void trivial_decomposition(G_t &G, T_t &T){
    typename boost::graph_traits<T_t>::vertex_descriptor t = boost::add_vertex(T);

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned id=noboost::get_id(G,*vIt);
        noboost::bag(T, t).insert(id);
    }
}

template <typename T_t>
int get_width(T_t &T){
    int max = -1;
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        size_t bag_size = noboost::bag(T, *tIt).size();
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
        avg += noboost::bag(T, *tIt).size();
    }

    return (boost::num_vertices(T) > 0)? avg/boost::num_vertices(T) : 0.0;
}

template <typename T_t>
unsigned int get_adhesion(T_t T){
    unsigned int max = 0;
    while(boost::num_vertices(T) > 1){
        typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            if(noboost::bag(T, *tIt).size() == 0 || boost::out_degree(*tIt, T) == 0){
                boost::clear_vertex(*tIt, T);
                boost::remove_vertex(*tIt, T);
                continue;
            }
            if(boost::out_degree(*tIt, T) == 1){
                typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                typename boost::graph_traits<T_t>::vertex_descriptor parent;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*tIt, T); nIt != nEnd; nIt++){
                    parent = *nIt;
                }

                std::set<unsigned int> intersection;
                std::set_intersection(noboost::bag(T, *tIt).begin(),
                                      noboost::bag(T, *tIt).end(),
                                      noboost::bag(T, parent).begin(),
                                      noboost::bag(T, parent).end(),
                                      std::inserter(intersection, intersection.begin()));

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
                if(*tIt == *nIt){
                    continue;
                }
                if(std::includes(noboost::bag(T, *nIt).begin(), noboost::bag(T, *nIt).end(),
                                 noboost::bag(T, *tIt).begin(), noboost::bag(T, *tIt).end())){
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
        if(std::includes(noboost::bag(T, *vIt).begin(),
                         noboost::bag(T, *vIt).end(),
                         bag.begin(), bag.end()))
        {
            if(noboost::bag(T, *vIt).find(elim_vertex) != noboost::bag(T, *vIt).end()){
                return;
            }

            t_dec_node = boost::add_vertex(T);
            noboost::bag(T,t_dec_node) = MOVE(bag);
            bag.clear();
            noboost::bag(T,t_dec_node).insert(elim_vertex);
            boost::add_edge(*vIt, t_dec_node, T);
            return;
        }
    }

    if(boost::num_vertices(T) > 0){
        boost::tie(vIt, vEnd) = boost::vertices(T);
    }

    t_dec_node = boost::add_vertex(T);
    noboost::bag(T,t_dec_node) = MOVE(bag);
    bag.clear();
    noboost::bag(T,t_dec_node).insert(elim_vertex);

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
        noboost::bag(T1, idxMap[id++]) = noboost::bag(T2, *tIt);
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
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        G[*vIt].id = id_map[k++];
    }
}

template <typename G_t>
void descriptor_bag_to_id_bag(G_t &G, std::set<unsigned int> &id_bag,
        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &desc_bag)
{
    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt =
         desc_bag.begin(); sIt != desc_bag.end(); sIt++)
    {
        unsigned id = noboost::get_id(G,*sIt);
        id_bag.insert(id);
    }
}

template <typename G_t>
void id_bag_to_descriptor_bag(G_t &G,
        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &desc_bag,
        std::set<unsigned int> &id_bag)
{
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned id=noboost::get_id(G,*vIt);
        if(id_bag.find(id) != id_bag.end()){
            desc_bag.insert(*vIt);
        }
    }
}

template <typename T_t>
void reorder_ids_decomposition(T_t &T, std::vector<unsigned int> &id_map){
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        std::set<unsigned int> new_bag;
        for(std::set<unsigned int>::iterator sIt = noboost::bag(T, *tIt).begin(); sIt != noboost::bag(T, *tIt).end(); sIt++){
            new_bag.insert(id_map[*sIt]);
        }

        noboost::bag(T, *tIt) = MOVE(new_bag);
    }
}

template <typename G_t>
void remove_isolated_vertices(G_t &H, G_t &G){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    unsigned int max = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned id=noboost::get_id(G,*vIt);
        max = (id > max)? id : max;
    }

    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap(max+1);
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) > 0){
            unsigned id=noboost::get_id(G,*vIt);
            idxMap[id] = boost::add_vertex(H);
            H[idxMap[id]].id = id;
        }
    }
    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt != eEnd; eIt++){
        unsigned sid=noboost::get_id(G, boost::source(*eIt, G));
        unsigned tid=noboost::get_id(G, boost::target(*eIt, G));
        boost::add_edge(idxMap[sid], idxMap[tid], H);
    }
}

#ifndef TD_SUBSETS
#define TD_SUBSETS

//collects all subsets of X of size k in subs (todo: replacement by enumeration in hunt())
static void subsets(std::set<unsigned int> &X, int /*size*/, int k, unsigned int idx,
        std::vector<unsigned int> &sub, std::vector<std::set<unsigned int> > &subs){
    if(k==0){
        std::set<unsigned int> subS;
        for(unsigned int i = 0; i < sub.size(); i++){
            subS.insert(sub[i]);
        }
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

#endif //TD_SUBSETS

#ifndef TD_POWERSET
#define TD_POWERSET

static void powerset(std::set<unsigned int> &X, std::vector<std::set<unsigned int> > &subs){
    std::vector<unsigned int> sub;
    for(unsigned int i = 0; i <=X.size(); i++){
        subsets(X, X.size(), i, 0, sub, subs);
    }
}

#endif //ifdef TD_POWERSET

} // namespace treedec

#endif //TD_MISC

// vim:ts=8:sw=4:et
