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
// Offers functionality to compute a tree decomposition of exact width.
//
// A tree decomposition is a graph that has a set of vertex indices as bundled property, e.g.:
//
// struct tree_dec_node
// {
//  std::set<unsigned int> bag;
// };
// typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, tree_dec_node> tree_dec_t;
//
// typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> TD_graph_t;
//
//

/*
  These functions are most likely to be interesting for outside use:

- void CR_dynamic_decomp(G_t &G, T_t &T, int lb)
- void CR_dynamic_decomp(G_t &G, T_t &T)
*/

#ifndef TD_DYNAMICCR
#define TD_DYNAMICCR

#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>

#include "simple_graph_algos.hpp"
#include "misc.hpp"
#include "graph.hpp"

namespace treedec{

//Computes the robber space-components with respect to X.
template <typename G_t>
void get_robber_components(G_t &G,
     typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &X,
     typename std::vector<typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type> &Rcomps)
{
    typedef typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type bag_type;

    //G \ X
    std::vector<bool> disabled(boost::num_vertices(G), false);
    for(typename bag_type::iterator sIt = X.begin(); sIt != X.end(); sIt++){
        unsigned int pos = noboost::get_pos(*sIt, G);
        disabled[pos] = true;
    }

    treedec::get_components_provided_map(G, Rcomps, disabled);
}


//Computes the robber space with respect to X.
template <typename G_t>
static void get_robber_component(
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &X_prime,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &R,
        typename std::vector<typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type> &Rcomps)
{
    for(unsigned int i = 0; i < Rcomps.size(); i++){
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type intersection;
        std::set_intersection(Rcomps[i].begin(), Rcomps[i].end(),
                              X_prime.begin(), X_prime.end(),
                              std::inserter(intersection, intersection.begin()));

        if(!intersection.empty()){
            R.insert(Rcomps[i].begin(), Rcomps[i].end());
        }
    }
}

//Checks if comp(G\X, y) = comp(G\ (X ^ X'), y) := R and X' ^ R != emptyset.
template <typename G_t>
bool is_monotone_dynamicCR(G_t &G,
          typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &X,
          typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &X_prime,
          typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &oldR,
          typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &newR,
          typename std::vector<typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type> &Rcomps)
{
    typedef typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type bag_type;

    if(X == X_prime){
        return false;
    }

    //robber_space(G\X)
    bag_type R1;
    get_robber_component<G_t>(X_prime, R1, Rcomps);



    //robber_space(G\(X ^ X'))
    bag_type IS_X_X_prime, R2;
    std::vector<bag_type> Rcomps2;
    std::set_intersection(X.begin(), X.end(),
                          X_prime.begin(), X_prime.end(),
                          std::inserter(IS_X_X_prime, IS_X_X_prime.begin()));

    get_robber_components(G, IS_X_X_prime, Rcomps2);



    get_robber_component<G_t>(X_prime, R2, Rcomps2);



    if(R1 == R2){
        if((oldR.size() == 0 && std::includes(X_prime.begin(), X_prime.end(), R1.begin(), R1.end()))
         || std::includes(oldR.begin(), oldR.end(), R1.begin(), R1.end()))
        {
            std::set_union(R1.begin(), R1.end(), X_prime.begin(), X_prime.end(), std::inserter(newR, newR.begin()));

            return true;
        }
    }

    return false;
}

template <typename G_t, typename W_t>
bool make_layer(G_t &G, W_t &W,
         typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type &vertices,
         unsigned int k, unsigned int idx)
{
    typedef typename noboost::treedec_chooser<G_t>::type T_t;
    typedef typename noboost::treedec_traits<T_t>::bag_type bag_type;
    typedef typename noboost::treedec_traits<T_t>::bag_type::value_type value_type;

    W.resize(idx+1);

    std::vector<value_type> sub;
    std::vector<bag_type> subs; //all possible X
    subsets(vertices, vertices.size(), k, 0, sub, subs);

    //in the last layer (idx == 0), we just insert all subsets, without
    //computing anything (the robber is caught)
    if(idx != 0){
        for(unsigned int i = 0; i <  subs.size(); i++){
            std::vector<unsigned int> indices;

            std::vector<bag_type> Rcomps;
            get_robber_components(G, subs[i], Rcomps);

            bag_type newR;

            bag_type forbidden;

            //search for all X', such that a turn is strict monotone
            for(unsigned int j = 0; j < W[idx-1].size(); j++){
                //if a turn X -> X' is monotone, save the index of the entry in
                //the layer below for computing a tree decomposition at the end

                // FIXME: don't create intersection. stop at first common elt.
                bag_type inters;
                std::set_intersection(subs[i].begin(), subs[i].end(), 
                                      W[idx-1][j].get<3>().begin(), W[idx-1][j].get<3>().end(), std::inserter(inters, inters.begin()));

                if(!inters.empty()){
                    continue;
                }

                if(is_monotone_dynamicCR(G, subs[i], W[idx-1][j].get<0>(), W[idx-1][j].get<1>(), newR, Rcomps)){
                    indices.push_back(j);
                    std::set_difference(W[idx-1][j].get<0>().begin(), W[idx-1][j].get<0>().end(), subs[i].begin(), subs[i].end(),
                                        std::inserter(forbidden, forbidden.begin()));
                }
            }
            //if there is at least one monotone turn including X, we have to add X to W
            if(indices.size() != 0){
                W[idx].push_back(boost::tuple<bag_type, bag_type, std::vector<unsigned int>, bag_type>(subs[i], newR, indices, forbidden));

                bag_type union_R_X;
                std::set_union(newR.begin(), newR.end(),
                               subs[i].begin(), subs[i].end(),
                               std::inserter(union_R_X, union_R_X.begin()));

                //test if a tree decomposition has been found
                if(union_R_X.size() == vertices.size()){
                    return true;
                }
            }
        }
    }

    //We caught the robber(y in X). To cut down memory usage, we save the empty
    //robber space instead of all possible y in X (or just X). With this
    //modification, we additionally have to check "if(|R|== 0 && X' includes
    //newR)" for monotonicity in is_monotone_dynamicCR

    //TODO: copy the whole "graph" in W displaced by one layer, such we would
    //not compute monotone turns several times
    for(unsigned int i = 0; i < subs.size(); i++){
        W[idx].push_back(boost::tuple<bag_type, bag_type, std::vector<unsigned int>, bag_type>
                  (subs[i], bag_type(), std::vector<unsigned int>(), bag_type()));
    }
    return false;
}

//Follows the "pointers" stored in the third entry of W[layer][i] recursively.
template <typename T_t, typename W_t, typename G_t>
void make_tree_decomposition(T_t &T, W_t &W, typename noboost::treedec_traits<T_t>::bag_type &R,
                             unsigned int layer, unsigned int idx, G_t &G)
{
    for(unsigned int i = 0; i < W[layer][idx].get<2>().size(); i++){
        if(std::includes(R.begin(), R.end(),
           W[layer-1][W[layer][idx].get<2>()[i]].get<0>().begin(),
           W[layer-1][W[layer][idx].get<2>()[i]].get<0>().end()))
        {
            continue;
        }

        R.insert(W[layer-1][W[layer][idx].get<2>()[i]].get<0>().begin(), W[layer-1][W[layer][idx].get<2>()[i]].get<0>().end());

        glue_two_bags(T, W[layer][idx].get<0>(), W[layer-1][W[layer][idx].get<2>()[i]].get<0>());
        make_tree_decomposition(T, W, R, layer-1, W[layer][idx].get<2>()[i], G);
    }
}

template <typename G_t, typename T_t>
void CR_dynamic_decomp(G_t &G, T_t &T, int lb){
    typedef typename noboost::treedec_traits<T_t>::bag_type bag_type;

    bag_type vertices;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    typename boost::graph_traits<G_t>::vertices_size_type nv=boost::num_vertices(G);
    typename boost::graph_traits<G_t>::edges_size_type ne=boost::num_edges(G);

    if( nv<=1 || 2*ne==nv*(nv-1u) ){
        typename boost::graph_traits<T_t>::vertex_descriptor t = boost::add_vertex(T);
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            noboost::bag(t, T).insert(*vIt);
        }
        return;
    }

    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        vertices.insert(*vIt);
    }

    std::vector<std::vector<
       boost::tuple<bag_type, bag_type,
                    std::vector<unsigned int>, bag_type> > > W(1);

    for(unsigned int k = (unsigned int)lb; k < vertices.size(); k++){
        for(unsigned int i = 0; i < vertices.size(); i++){
            if(make_layer(G, W, vertices, k, i)){
                typename noboost::treedec_traits<T_t>::bag_type R;
                R.insert(W[W.size()-1][W.back().size()-1].get<0>().begin(), W[W.size()-1][W.back().size()-1].get<0>().end());
                make_tree_decomposition(T, W, R, W.size()-1, W.back().size()-1, G);
                return;
            }
        }

        W.clear();
        W.resize(1);
    }
}


template <typename G_t, typename T_t>
static void CR_dynamic_decomp(G_t &G, T_t &T){
    int lb = 0;
    CR_dynamic_decomp(G, T, lb);
}

} //namespace treedec

#endif //ifdef TD_DYNAMICCR

// vim:ts=8:sw=4:et
