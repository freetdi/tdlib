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
// Offers some recommended combinations of the algorithms
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
//
// These functions are most likely to be interesting for outside use:
//
// void preprocessing_MD(G_t &G, T_t &T)
// void preprocessing_MD(G_t &G, T_t &T, int &low)
// void preprocessing_FI_TM(G_t &G, T_t &T)
// void preprocessing_FI_TM(G_t &G, T_t &T, int &low)
//

#ifndef TD_COMBINATIONS
#define TD_COMBINATIONS

#include <set>
#include <boost/graph/adjacency_list.hpp>
#include "TD_preprocessing.hpp"
#include "TD_elimination_orderings.hpp"
#include "TD_postprocessing.hpp"

namespace treedec{

//recursively applies preprocessing rules and glues corresponding bags with current tree decomposition
//this version applies the minDegree-heuristic on not fully preprocessable graph instances
template <typename G_t, typename T_t>
void preprocessing_MD(G_t &G, T_t &T, int &low){
    std::set<unsigned int> bag;
    unsigned int preprocessed_node;
    if(Islet(G, bag, preprocessed_node, low) || Twig(G, bag, preprocessed_node, low) || Series(G, bag, preprocessed_node, low) || 
       Triangle(G, bag, preprocessed_node, low) || Buddy(G, bag, preprocessed_node, low) || Cube(G, bag, preprocessed_node, low) ||
       Simplicial(G, bag, preprocessed_node, low) || AlmostSimplicial(G, bag, preprocessed_node, low)){

        preprocessing_MD(G, T, low);
        _glue_bag_preprocessing(bag, preprocessed_node, T);
        return;
    }
    if(boost::num_vertices(G) != 0){
        treedec::minDegree_decomp(G, T);
        low = (low > 4)? low : 4;
    }
}


//recursively applies preprocessing rules and glues corresponding bags with current tree decomposition
//this version applies the fillIn-heuristic on not fully preprocessable graph instances
template <typename G_t, typename T_t>
void preprocessing_FI_TM(G_t &G, T_t &T, int &low){
    std::set<unsigned int> bag;
    unsigned int preprocessed_node;
    if(Islet(G, bag, preprocessed_node, low) || Twig(G, bag, preprocessed_node, low) || Series(G, bag, preprocessed_node, low) || 
       Triangle(G, bag, preprocessed_node, low) || Buddy(G, bag, preprocessed_node, low) || Cube(G, bag, preprocessed_node, low) ||
       Simplicial(G, bag, preprocessed_node, low) || AlmostSimplicial(G, bag, preprocessed_node, low)){

        preprocessing_FI_TM(G, T, low);
        _glue_bag_preprocessing(bag, preprocessed_node, T);       
        return;
    }
    if(boost::num_vertices(G) != 0){
        std::vector<unsigned int> old_elim_ordering2, new_elim_ordering2;
        treedec::fillIn_ordering(G, old_elim_ordering2);
        treedec::minimalChordal(G, old_elim_ordering2, new_elim_ordering2);
        treedec::ordering_to_treedec(G, new_elim_ordering2, T);
        low = (low > 4)? low : 4;
    }
}

template <typename G_t, typename T_t>
void preprocessing_MD(G_t &G, T_t &T){
    int low = -1;
    preprocessing_MD(G, T, low);
}

template <typename G_t, typename T_t>
void preprocessing_FI_TM(G_t &G, T_t &T){
    int low = -1;
    preprocessing_FI_TM(G, T, low);
}

}

#endif

