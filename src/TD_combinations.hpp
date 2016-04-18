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

/*
 Offers some recommended combinations of the algorithms.

 These functions are most likely to be interesting for outside use:

 - void PP_MD(G_t &G, T_t &T)
 - void PP_MD(G_t &G, T_t &T, int &low)
 - void PP_FI(G_t &G, T_t &T)
 - void PP_FI(G_t &G, T_t &T, int &low)
 - void PP_FI_TM(G_t &G, T_t &T)
 - void PP_FI_TM(G_t &G, T_t &T, int &low)
 - void FI_TM(G_t &G, T_t &T)
 - void FI_TM(G_t &G, T_t &T, int &low)
 - void exact_decomposition_cutset(G_t &G, T_t &T)
 - void exact_decomposition_cutset(G_t &G, T_t &T, int low)
 - void exact_decomposition_cutset_decision(G_t &G, T_t &T, int k)
 - void exact_decomposition_dynamic(G_t &G, T_t &T)
 - void exact_decomposition_dynamic(G_t &G, T_t &T, int low)
 - void seperator_algorithm_MSVS(G_t &G, T_t &T)
 - void seperator_algorithm_TM(G_t &G, T_t &T)
 - void MSVS_trivial(G_t &G, T_t &T)
*/

#ifndef TD_COMBINATIONS
#define TD_COMBINATIONS

#include <set>
#include <vector>
#include <boost/graph/adjacency_list.hpp>

#include "TD_preprocessing.hpp"
#include "TD_lower_bounds.hpp"
#include "TD_elimination_orderings.hpp"
#include "TD_postprocessing.hpp"
#include "TD_dynamicCR.hpp"
#include "TD_exact_cutset.hpp"
#include "TD_seperator_algorithm.hpp"
#include "TD_misc.hpp"

namespace treedec{

//Recursively applies preprocessing rules and glues corresponding bags with current tree decomposition
//this version applies the minDegree-heuristic on not fully preprocessable graph instances.
template <typename G_t, typename T_t>
void PP_MD(G_t &G, T_t &T, int &low){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > bags;

    treedec::preprocessing(G, bags, low);
    if(boost::num_edges(G) > 0){
        impl::minDegree_decomp(G, &T);
    }
    treedec::preprocessing_glue_bags(bags, T);
}

//Recursively applies preprocessing rules and glues corresponding bags with current tree decomposition
//this version applies the minDegree-heuristic on not fully preprocessable graph instances.
template <typename G_t, typename T_t>
void PP_FI(G_t &G, T_t &T, int &low){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > bags;

    treedec::preprocessing(G, bags, low);
    if(boost::num_edges(G) > 0){
        treedec::impl::fillIn_decomp(G, &T);
    }
    treedec::preprocessing_glue_bags(bags, T);
}


//Recursively applies preprocessing rules and glues corresponding bags with current tree decomposition.
//This version applies the fillIn-heuristic followed by triangulation minimization on not fully preprocessable graph instances.
template <typename G_t, typename T_t>
void PP_FI_TM(G_t &G, T_t &T, int &low){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > bags;

    treedec::preprocessing(G, bags, low);

    if(boost::num_edges(G) > 0){
        typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> old_elim_ordering;
        typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> new_elim_ordering;

        treedec::fillIn_ordering(G, old_elim_ordering, true); //true = ignore isolated vertices
#ifndef NDEBUG
        for( auto i : old_elim_ordering){
            assert(noboost::is_valid(i,G));
        }
#endif
        treedec::minimalChordal(G, old_elim_ordering, new_elim_ordering);
        treedec::ordering_to_treedec(G, new_elim_ordering, T, true); //true = ignore isolated vertices
    }

    treedec::preprocessing_glue_bags(bags, T);
}

//This version applies the fillIn-heuristic followed by triangulation minimization on the input graph.
template <typename G_t, typename T_t>
void FI_TM(G_t &G, T_t &T){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> old_elim_ordering;
    typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> new_elim_ordering;
    treedec::fillIn_ordering(G, old_elim_ordering);
    treedec::minimalChordal(G, old_elim_ordering, new_elim_ordering);
    treedec::ordering_to_treedec(G, new_elim_ordering, T);
}


template <typename G_t, typename T_t>
void exact_decomposition_cutset(G_t &G, T_t &T, int lb){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    //Preprocessing.
    int low = -1;

    std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > bags;
    treedec::preprocessing(G, bags, low);

    if(boost::num_edges(G) == 0){
        treedec::preprocessing_glue_bags(bags, T);
        return;
    }

    //Lower bound on the treewidth of the reduced instance of G.
    G_t H(G);
    int lb_deltaC = treedec::lb::deltaC_least_c(H);

    lb = (low > lb)? low : lb;
    lb = (lb_deltaC > lb)? lb_deltaC : lb;

    //Compute a treedecomposition for each connected component of G and glue the decompositions together.
    std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > components;
    treedec::get_components(G, components);

    typename boost::graph_traits<T_t>::vertex_descriptor root = boost::add_vertex(T);

    for(unsigned int i = 0; i < components.size(); i++){
        //Ignore isolated vertices (already included in 'bags').
        if(components[i].size() == 1){
            continue;
        }

        G_t G_;
        typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> vdMap;
        treedec::induced_subgraph(G_, G, components[i], vdMap);
        T_t T_;

        while(!treedec::exact_cutset(G_, T_, lb)){
            lb++;
        }

        treedec::apply_map_on_treedec(T_, G_, vdMap);

        treedec::glue_decompositions(T, T_);
    }

    treedec::preprocessing_glue_bags(bags, T);
}


template <typename G_t, typename T_t>
bool exact_decomposition_cutset_decision(G_t &G, T_t &T, int k){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        if(k >= -1){ return true; }
        else{ return false; }
    }

    //Preprocessing.
    int low = -1;

    std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > bags;
    treedec::preprocessing(G, bags, low);

    if(boost::num_edges(G) == 0){
        treedec::preprocessing_glue_bags(bags, T);
        if(low <= k){
            return true;
        }
        return false;
    }

    //Lower bound on the treewidth of the reduced instance of G.
    G_t H(G);
    int lb_deltaC = treedec::lb::deltaC_least_c(H);

    int lb = low;
    lb = (lb_deltaC > lb)? lb_deltaC : lb;

    if(lb > k){
        return false;
    }

    //Compute a treedecomposition for each connected component of G and glue the decompositions together.
    std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > components;
    treedec::get_components(G, components);

    typename boost::graph_traits<T_t>::vertex_descriptor root = boost::add_vertex(T);

    for(unsigned int i = 0; i < components.size(); i++){
        //Ignore isolated vertices (already included in 'bags').
        if(components[i].size() == 1){
            continue;
        }

        G_t G_;
        typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> vdMap;
        treedec::induced_subgraph(G_, G, components[i], vdMap);
        T_t T_;

        while(!treedec::exact_cutset(G, T_, lb)){
            lb++;
            if(lb > k){
                return false;
            }
        }
    }

    return true;
}

/*
template <typename G_t, typename T_t>
void exact_decomposition_dynamic(G_t &G, T_t &T, int lb){
    //preprocessing
    int low = -1;
    std::vector<boost::tuple<
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::vd_type,
        typename noboost::treedec_traits<typename noboost::treedec_chooser<G_t>::type>::bag_type
         > > bags;

    treedec::preprocessing(G, bags, low);
    if(boost::num_edges(G) == 0){
        treedec::preprocessing_glue_bags(bags, T);
        return;
    }

    lb = (low > lb)? low : lb;

    //Compute a treedecomposition for each connected component of G and glue the decompositions together.
    std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > components;
    treedec::get_components(G, components);

    if(components.size() == 1){
        treedec::CR_dynamic_decomp(G, T, lb);

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
*/

template <typename G_t, typename T_t>
void exact_decomposition_chordal(G_t &G, T_t &T){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> elim_ordering;
    treedec::LEX_M_minimal_ordering(G, elim_ordering);
    treedec::ordering_to_treedec(G, elim_ordering, T);
}

template <typename G_t, typename T_t>
void seperator_algorithm_MSVS(G_t &G, T_t &T){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    treedec::seperator_algorithm(G, T);
    treedec::MSVS(G, T);
}

template <typename G_t, typename T_t>
void seperator_algorithm_TM(G_t &G, T_t &T){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    treedec::seperator_algorithm(G, T);
    typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> old_elim_ordering;
    typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> new_elim_ordering;
    treedec::treedec_to_ordering<G_t, T_t>(T, old_elim_ordering);
    treedec::minimalChordal(G, old_elim_ordering, new_elim_ordering);
    T.clear();
    treedec::ordering_to_treedec(G, new_elim_ordering, T);
}

template <typename G_t, typename T_t>
void MSVS_trivial(G_t &G, T_t &T){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    treedec::trivial_decomposition(G, T);
    treedec::MSVS(G, T);
}

template <typename G_t, typename T_t>
void PP_MD(G_t &G, T_t &T){
    int low = -1;
    PP_MD(G, T, low);
}

template <typename G_t, typename T_t>
void PP_FI(G_t &G, T_t &T){
    int low = -1;
    PP_FI(G, T, low);
}

template <typename G_t, typename T_t>
void PP_FI_TM(G_t &G, T_t &T){
    int low = -1;
    PP_FI_TM(G, T, low);
}

/*
template <typename G_t, typename T_t>
void exact_decomposition_dynamic(G_t &G, T_t &T){
    int lb = -1;
    exact_decomposition_dynamic(G, T, lb);
}
*/

template <typename G_t, typename T_t>
void exact_decomposition_cutset(G_t &G, T_t &T){
    int lb = -1;
    exact_decomposition_cutset(G, T, lb);
}

} //namespace treedec

#endif //TD_COMBINATIONS

// vim:ts=8:sw=4:et
