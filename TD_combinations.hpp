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
// void preprocessing_FI(G_t &G, T_t &T)
// void preprocessing_FI(G_t &G, T_t &T, int &low)
// void preprocessing_FI_TM(G_t &G, T_t &T)
// void preprocessing_FI_TM(G_t &G, T_t &T, int &low)
// void FI_TM(G_t &G, T_t &T)
// void FI_TM(G_t &G, T_t &T, int &low)
// void exact_decomposition_cutset(G_t &G, T_t &T)
// void exact_decomposition_cutset(G_t &G, T_t &T, int low)
// void exact_decomposition_dynamic(G_t &G, T_t &T)
// void exact_decomposition_dynamic(G_t &G, T_t &T, int low)
//

#ifndef TD_COMBINATIONS
#define TD_COMBINATIONS

#include <set>
#include <boost/graph/adjacency_list.hpp>
#include "TD_preprocessing.hpp"
#include "TD_lower_bounds.hpp"
#include "TD_elimination_orderings.hpp"
#include "TD_postprocessing.hpp"
#include "TD_dynamicCR.hpp"
#include "TD_exact_cutset.hpp"

namespace treedec{

//recursively applies preprocessing rules and glues corresponding bags with current tree decomposition
//this version applies the minDegree-heuristic on not fully preprocessable graph instances
template <typename G_t, typename T_t>
void preprocessing_MD(G_t &G, T_t &T, int &low){
    std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
    treedec::preprocessing(G, bags, low);
    if(boost::num_edges(G) != 0)
        treedec::minDegree_decomp(G, T);
    treedec::preprocessing_glue_bags(bags, T);
}

//recursively applies preprocessing rules and glues corresponding bags with current tree decomposition
//this version applies the minDegree-heuristic on not fully preprocessable graph instances
template <typename G_t, typename T_t>
void preprocessing_FI(G_t &G, T_t &T, int &low){
    std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
    treedec::preprocessing(G, bags, low);
    treedec::fillIn_decomp(G, T);
    treedec::preprocessing_glue_bags(bags, T);
}

//recursively applies preprocessing rules and glues corresponding bags with current tree decomposition
//this version applies the fillIn-heuristic followed by triangulation minimization on not fully preprocessable graph instances
template <typename G_t, typename T_t>
void preprocessing_FI_TM(G_t &G, T_t &T, int &low){
    std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
    treedec::preprocessing(G, bags, low);
    treedec::remove_isolated_vertices(G);
    if(boost::num_edges(G) != 0){
        std::vector<unsigned int> old_elim_ordering, new_elim_ordering;
        treedec::fillIn_ordering(G, old_elim_ordering);
        treedec::minimalChordal(G, old_elim_ordering, new_elim_ordering);
        treedec::ordering_to_treedec(G, new_elim_ordering, T);
    }
    treedec::preprocessing_glue_bags(bags, T); 
}

//this version applies the fillIn-heuristic followed by triangulation minimization
template <typename G_t, typename T_t>
void FI_TM(G_t &G, T_t &T, int &low){
    std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
    treedec::Islet(G, bags, low);
    if(boost::num_edges(G) != 0){
        treedec::remove_isolated_vertices(G);
        std::vector<unsigned int> old_elim_ordering, new_elim_ordering;
        treedec::fillIn_ordering(G, old_elim_ordering);
        treedec::minimalChordal(G, old_elim_ordering, new_elim_ordering);
        treedec::ordering_to_treedec(G, new_elim_ordering, T);
    }
    treedec::preprocessing_glue_bags(bags, T); 
}

template <typename G_t, typename T_t>
void exact_decomposition_dynamic(G_t &G, T_t &T, int lb){
    //preprocessing
    int low = -1;
    std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
    treedec::preprocessing(G, bags, low);
    if(boost::num_vertices(G) == 0){
        treedec::preprocessing_glue_bags(bags, T);
        return;
    }

    lb = (low > lb)? low : lb;

    //compute a treedecomposition for each connected component of G and glue the decompositions together
    std::vector<std::set<unsigned int> > components;
    get_components(G, components);

    if(components.size() == 1){
        //reorder ids
        std::vector<unsigned int> id_map;
        treedec::reorder_ids_graph(G, id_map);

        treedec::CR_dynamic_decomp(G, T, lb);

        treedec::reorder_ids_decomposition(T, id_map);

        treedec::preprocessing_glue_bags(bags, T);
        return;
    }

    typename boost::graph_traits<T_t>::vertex_descriptor root = boost::add_vertex(T);

    for(unsigned int i = 0; i < components.size(); i++){
        //ignore isolated vertices
        if(components[i].size() == 1)
            continue;

        G_t G_;
        induced_subgraph(G_, G, components[i]);
        
        //reorder ids
        std::vector<unsigned int> id_map;
        treedec::reorder_ids_graph(G_, id_map);

        T_t T_;

        treedec::CR_dynamic_decomp(G_, T_, lb);

        treedec::reorder_ids_decomposition(T_, id_map);

        treedec::glue_decompositions(T, T_);
    }
    
    treedec::preprocessing_glue_bags(bags, T);
}

template <typename G_t, typename T_t>
void exact_decomposition_cutset(G_t &G, T_t &T, int lb){
    //preprocessing
    int low = -1;

    std::vector<boost::tuple<unsigned int, std::set<unsigned int> > > bags;
    treedec::preprocessing(G, bags, low);

    if(boost::num_edges(G) == 0){
        treedec::preprocessing_glue_bags(bags, T);
        return;
    }

    int lb_deltaC = treedec::lb::deltaC_least_c(G);

    lb = (low > lb)? low : lb;
    lb = (lb_deltaC > lb)? lb_deltaC : lb;

    //compute a treedecomposition for each connected component of G and glue the decompositions together
    std::vector<std::set<unsigned int> > components;
    get_components(G, components);

    typename boost::graph_traits<T_t>::vertex_descriptor root = boost::add_vertex(T);

    for(unsigned int i = 0; i < components.size(); i++){
        //ignore isolated vertices
        if(components[i].size() == 1)
            continue;

        G_t G_;
        induced_subgraph(G_, G, components[i]);

        //reorder ids
        std::vector<unsigned int> id_map;
        treedec::reorder_ids_graph(G_, id_map);

        T_t T_;

        treedec::exact_cutset(G_, T_, lb);

        treedec::reorder_ids_decomposition(T_, id_map);

        treedec::glue_decompositions(T, T_);
    }
    
    treedec::preprocessing_glue_bags(bags, T);
}

template <typename G_t, typename T_t>
void exact_decomposition_chordal(G_t &G, T_t &T){
    std::vector<unsigned int> elim_ordering;
    treedec::LEX_M_minimal_ordering(G, elim_ordering);
    treedec::ordering_to_treedec(G, elim_ordering, T);
}

template <typename G_t, typename T_t>
void preprocessing_MD(G_t &G, T_t &T){
    int low = -1;
    preprocessing_MD(G, T, low);
}

template <typename G_t, typename T_t>
void preprocessing_FI(G_t &G, T_t &T){
    int low = -1;
    preprocessing_FI(G, T, low);
}

template <typename G_t, typename T_t>
void preprocessing_FI_TM(G_t &G, T_t &T){
    int low = -1;
    preprocessing_FI_TM(G, T, low);
}

template <typename G_t, typename T_t>
void FI_TM(G_t &G, T_t &T){
    int low = -1;
    FI_TM(G, T, low);
}

template <typename G_t, typename T_t>
void exact_decomposition_dynamic(G_t &G, T_t &T){
    int lb = -1;
    exact_decomposition_dynamic(G, T, lb);
}

template <typename G_t, typename T_t>
void exact_decomposition_cutset(G_t &G, T_t &T){
    int lb = -1;
    exact_decomposition_cutset(G, T, lb);
}

}

#endif

