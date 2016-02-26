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
// typedef boost::adjacency_list<boost::setS, boost::listS, boost::undirectedS, Vertex> TD_graph_t;
//
//


#ifndef TD_APPROXIMATIOM
#define TD_APPROXIMATION

#include <set>
#include "TD_approximate_seperator.hpp"
#include "TD_misc.hpp"
#include "TD_simple_graph_algos.hpp"

namespace treedec{

//the main procedure of the log-Seperator algorithm
template <typename G_t, typename T_t>
typename boost::graph_traits<T_t>::vertex_descriptor logApproximation(G_t &G, T_t &T, typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &Z, typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &W){

    typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> union_Z_W;
    std::set_union(W.begin(), W.end(), Z.begin(), Z.end(), std::inserter(union_Z_W, union_Z_W.begin()));

    //return a tree decomposition with one single node, containing (Z v W)
    if(3*Z.size() <= W.size()){
        std::set<unsigned int> bag;
        treedec::descriptor_bag_to_id_bag(G, bag, union_Z_W);
        typename boost::graph_traits<T_t>::vertex_descriptor t_dec_node = boost::add_vertex(T);
        T[t_dec_node].bag = bag;

        return t_dec_node;
    }

    typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> S, S_;

    std::vector<bool> disabled(boost::num_vertices(G), true);
    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = union_Z_W.begin(); sIt != union_Z_W.end(); sIt++)
        disabled[G[*sIt].id] = false;

    //find a 2/3-vertex seperator S of W in G[Z v W]
    std::vector<bool> disabled1(disabled);
    approximate_vertex_seperator(G, disabled1, W, S);

    //find a 2/3-vertex seperator S_ of (Z v W) in G[Z v W]
    std::vector<bool> disabled2(disabled);
    approximate_vertex_seperator(G, disabled, union_Z_W, S_);


    //compute the connected components of G[Z v W - (S v S_)]
    std::vector<bool> disabled3(disabled);
    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = S.begin(); sIt != S.end(); sIt++)
        disabled[G[*sIt].id] = true;
    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = S_.begin(); sIt != S_.end(); sIt++)
        disabled[G[*sIt].id] = true;

    std::vector<typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > C;
    get_components_provided_map(G, C, disabled3);

    //create a new root node with a bag, containing vertices (W v S v S')
    typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> union_S_S_, union_W_S_S_;

    std::set_union(S.begin(), S.end(), S_.begin(), S_.end(), std::inserter(union_S_S_, union_S_S_.begin()));
    std::set_union(W.begin(), W.end(), union_S_S_.begin(), union_S_S_.end(), std::inserter(union_W_S_S_, union_W_S_S_.begin()));

    typename boost::graph_traits<T_t>::vertex_descriptor t_dec_node = boost::add_vertex(T);
    std::set<unsigned int> root_bag;
    treedec::descriptor_bag_to_id_bag(G, root_bag, union_W_S_S_);
    T[t_dec_node].bag = root_bag;

    //let G_1,...,G_|C| be the connected components of G[Z v W - (S v S')]
    //call logApproximation with Z = Zi and W = (Wi v S v S_)
    //add an edge form the root of the current call to the roots of the recursive calls
    for(unsigned int i = 0; i < C.size(); i++){
        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> Zi, Wi, union_Wi_S_S_;

        //Zi = (V(Gi) ^ Z)
        std::set_intersection(C[i].begin(), C[i].end(), Z.begin(), Z.end(), std::inserter(Zi, Zi.begin()));

        //Wi = (V(Gi) ^ W)
        std::set_intersection(C[i].begin(), C[i].end(), W.begin(), W.end(), std::inserter(Wi, Wi.begin()));

        typename boost::graph_traits<T_t>::vertex_descriptor t_dec_rec_node = logApproximation(G, T, Zi, union_Wi_S_S_);
        boost::add_edge(t_dec_node, t_dec_rec_node, T);
    }
}

template <typename G_t, typename T_t>
void logApproximation(G_t &G, T_t &T){
    typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> V, emptyset;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        V.insert(*vIt);

    logApproximation(G, T, V, emptyset);
}

} //namespace treedec

#endif //ifdef TD_APPROXIMATIOM

// vim:ts=8:sw=4:et
