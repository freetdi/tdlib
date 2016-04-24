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
// Offers functionality to compute a tree decomposition of exact width.
//
// A tree decomposition is a graph that has a set of vertex indices as bundled property, e.g.:
//
// struct tree_dec_node
// {
//  std::set<unsigned int> bag;
// };
// typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, tree_dec_node> tree_dec_t;
//
// Vertices of the input graph have to provide the attribute 'id', e.g.:
//
// struct Vertex
// {
//  unsigned int id;
// };
// typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, Vertex> TD_graph_t;
//
//
//
// These functions are most likely to be interesting for outside use:
//
// void CR_dynamic_decomp(G_t &G, T_t &T, int lb)
// void CR_dynamic_decomp(G_t &G, T_t &T)
//

#ifndef TD_DYNAMICCR
#define TD_DYNAMICCR

#include <map>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>

#include "simple_graph_algos.hpp"
#include "misc.hpp"
#include "noboost.hpp"

namespace treedec{

//computes the robber space with respect to X and y and saves it in R
template <typename G_t>
void get_robber_components(G_t G, std::set<unsigned int> &X, std::vector<std::set<unsigned int> > &Rcomps){
    //G \ X
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(std::set<unsigned int>::iterator sIt = X.begin(); sIt != X.end(); sIt++){
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int pos = noboost::get_pos(*vIt, G);
            if(pos == *sIt){
                boost::clear_vertex(*vIt, G);
                boost::remove_vertex(*vIt, G);
                break;
            }
        }
    }

    get_components(G, Rcomps);
}


//computes the robber space with respect to X and y and saves it in R
static void get_robber_component(std::set<unsigned int> &X_prime, std::set<unsigned int> &R, std::vector<std::set<unsigned int> > &Rcomps){
    for(unsigned int i = 0; i < Rcomps.size(); i++){
        std::set<unsigned int> intersection;
        std::set_intersection(Rcomps[i].begin(), Rcomps[i].end(), X_prime.begin(), X_prime.end(),
                              std::inserter(intersection, intersection.begin()));
        if(!intersection.empty()){
            R.insert(Rcomps[i].begin(), Rcomps[i].end());
        }
    }
}

//checks if comp(G\X, y) = comp(G\ (X ^ X'), y) := R and X' ^ R != emptyset
template <typename G_t>
bool is_monotone_dynamicCR(G_t &G, std::set<unsigned int> &X, std::set<unsigned int> &X_prime, std::set<unsigned int> &oldR, std::set<unsigned int> &newR, std::vector<std::set<unsigned int> > &Rcomps){ 
    if(X == X_prime){
        return false;
    }

    //robber_space(G\X)
    std::set<unsigned int> R1;
    get_robber_component(X_prime, R1, Rcomps);

    //robber_space(G\(X ^ X'))
    std::set<unsigned int> is_X_X_prime;
    std::set_intersection(X.begin(), X.end(), X_prime.begin(), X_prime.end(), std::inserter(is_X_X_prime, is_X_X_prime.begin()));
    std::set<unsigned int> R2;
    std::vector<std::set<unsigned int> > Rcomps2;
    get_robber_components(G, is_X_X_prime, Rcomps2);
    get_robber_component(X_prime, R2, Rcomps2);

    if(R1 == R2 && ((oldR.size() == 0 && std::includes(X_prime.begin(), X_prime.end(), R1.begin(), R1.end()))
    || std::includes(oldR.begin(), oldR.end(), R1.begin(), R1.end()))){
        std::set_union(R1.begin(), R1.end(), X_prime.begin(), X_prime.end(), std::inserter(newR, newR.begin()));
        return true;
    }

    return false;
}

template <typename G_t>
bool make_layer(G_t &G, std::vector<std::vector<boost::tuple<std::set<unsigned int>, std::set<unsigned int>, std::vector<unsigned int> > > > &W, std::set<unsigned int> &vertices, unsigned int k, unsigned int idx){
    W.resize(idx+1);

    std::vector<unsigned int> sub;
    std::vector<std::set<unsigned int> > subs; //all possible X
    subsets(vertices, vertices.size(), k, 0, sub, subs);

    //in the last layer (idx == 0), we just insert all subsets, without computing anything (the robber is catched)
    if(idx != 0){
        for(unsigned int i = 0; i <  subs.size(); i++){
            std::vector<unsigned int> indices;

            std::vector<std::set<unsigned int> > Rcomps;
            get_robber_components(G, subs[i], Rcomps);
            std::set<unsigned int> newR;

            //search for all X', such that a turn is strict monotone
            for(unsigned int j = 0; j < W[idx-1].size(); j++){
                //if a turn X -> X' is monotone, save the index of the entry in the layer below for computing a tree decomposition
                //at the end
                if(is_monotone_dynamicCR(G, subs[i], W[idx-1][j].get<0>(), W[idx-1][j].get<1>(), newR, Rcomps)){
                    indices.push_back(j);
                }
            }
            //if there is a at least one monotone turn, including X, we have to add X to W
            if(indices.size() != 0){
                W[idx].push_back(boost::tuple<std::set<unsigned int>, std::set<unsigned int>, std::vector<unsigned int> >(subs[i], newR, indices));

                std::set<unsigned int> union_R_X;
                std::set_union(newR.begin(), newR.end(), subs[i].begin(), subs[i].end(), std::inserter(union_R_X, union_R_X.begin()));

                //test if a tree decomposition has been found
                if(union_R_X.size() == vertices.size()){
                    return true;
                }
            }
        }
    }

    //We catched the robber(y in X). To cut down memory usage, we save the empty robber space instead of all possible y in X (or just X).
    //With this modification, we additionally have to check "if(|R|== 0 && X' includes newR)" for monotonicity in is_monotone_dynamicCR

    //todo: copy the whole "graph" in W displaced by one layer, such we would not compute monotone turns several times
    for(unsigned int i = 0; i < subs.size(); i++){
        W[idx].push_back(boost::tuple<std::set<unsigned int>, std::set<unsigned int>,
                         std::vector<unsigned int> >(subs[i], std::set<unsigned int>(),
                         std::vector<unsigned int>()));
    }
    return false;
}


template <typename T_t>
void dynamicCR_glue_bags(T_t &T, std::set<unsigned int> bag1, std::set<unsigned int> &bag2){
    typename boost::graph_traits<T_t>::vertex_iterator vIt1, vIt2, vEnd;
    typename boost::graph_traits<T_t>::vertex_descriptor b1, b2;

    for(boost::tie(vIt1, vEnd) = boost::vertices(T); vIt1 != vEnd; vIt1++){
        if(noboost::bag(*vIt1, T) == bag1){
            b1 = *vIt1;
            break;
        }
    }

    for(boost::tie(vIt2, vEnd) = boost::vertices(T); vIt2 != vEnd; vIt2++){
        if(noboost::bag(*vIt2, T) == bag2){
            b2 = *vIt2;
            break;
        }
    }

    if(vIt1 != vEnd && vIt2 != vEnd){
        return;
    }

    if(vIt1 == vEnd){
        b1 = boost::add_vertex(T);
        noboost::bag(b1, T) = bag1;
    }

    if(vIt2 == vEnd){
        b2 = boost::add_vertex(T);
        noboost::bag(b2, T) = bag2;
    }

    boost::add_edge(b1, b2, T);
}

//follows the "pointers" stored in the third entry of W[layer][i] recursively
template <typename T_t>
void make_tree_decomposition(T_t &T, std::vector<std::vector<boost::tuple<std::set<unsigned int>, std::set<unsigned int>, std::vector<unsigned int> > > > &W, unsigned int layer, unsigned int idx){
    std::set<unsigned int> R;
    for(unsigned int i = 0; i < W[layer][idx].get<2>().size(); i++){
        if(std::includes(R.begin(), R.end(), W[layer-1][W[layer][idx].get<2>()[i]].get<0>().begin(), W[layer-1][W[layer][idx].get<2>()[i]].get<0>().end()))
            continue;

        for(std::set<unsigned int>::iterator sIt = W[layer-1][W[layer][idx].get<2>()[i]].get<0>().begin(); sIt != W[layer-1][W[layer][idx].get<2>()[i]].get<0>().end(); sIt++)
            R.insert(*sIt);
        dynamicCR_glue_bags(T, W[layer][idx].get<0>(), W[layer-1][W[layer][idx].get<2>()[i]].get<0>());
        make_tree_decomposition(T, W, layer-1, W[layer][idx].get<2>()[i]);
    }
}

template <typename G_t, typename T_t>
void CR_dynamic_decomp(G_t &G, T_t &T, int lb){
    std::set<unsigned int> vertices;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    if(boost::num_vertices(G) <= 1 || 2*boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        typename boost::graph_traits<T_t>::vertex_descriptor t = boost::add_vertex(T);
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int pos = noboost::get_pos(*vIt, G);
            noboost::bag(t, T).insert(pos);
        }

        return;
    }

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int pos = noboost::get_pos(*vIt, G);
        vertices.insert(pos);
    }

    std::vector<std::vector<boost::tuple<std::set<unsigned int>, std::set<unsigned int>, std::vector<unsigned int> > > > W(1);
    unsigned int k = (unsigned int)lb;

    while(true){
        k++;
        for(unsigned int i = 0; i < vertices.size(); i++){
            if(make_layer(G, W, vertices, k, i)){
                make_tree_decomposition(T, W, W.size()-1, W.back().size()-1);
                return;
            }
        }

        W.clear();
        W.resize(1);
    }
}


template <typename G_t, typename T_t>
static void CR_dynamic_decomp(G_t &G, T_t &T){
    int lb = -1;
    CR_dynamic_decomp(G, T, lb);
}

} //namespace treedec

#endif //ifdef TD_DYNAMICCR

// vim:ts=8:sw=4:et
