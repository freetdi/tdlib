// Copyright (C) 2014-2017 Lukas Larisch
// Author: Lukas Larisch
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
 * Offers some recommended combinations of the algorithms.
 *
 * These functions are most likely to be interesting for outside use:
 *
 * - void PP_MD(G_t &G, T_t &T)
 * - void PP_MD(G_t &G, T_t &T, int &low)
 * - void PP_FI(G_t &G, T_t &T)
 * - void PP_FI(G_t &G, T_t &T, int &low)
 * - void PP_FI_TM(G_t &G, T_t &T)
 * - void PP_FI_TM(G_t &G, T_t &T, int &low)
 *
 * - void exact_decomposition_cutset(G_t &G, T_t &T)
 * - void exact_decomposition_cutset(G_t &G, T_t &T, int low)
 * - void exact_decomposition_cutset_decision(G_t &G, T_t &T, int k)
 * - void exact_decomposition_dynamic(G_t &G, T_t &T)
 * - void exact_decomposition_dynamic(G_t &G, T_t &T, int low)
 *
 * - void separator_algorithm_MSVS(G_t &G, T_t &T)
 * - void separator_algorithm_TM(G_t &G, T_t &T)
 * - void MSVS_trivial(G_t &G, T_t &T)
 *
 */

#ifndef TD_COMBINATIONS
#define TD_COMBINATIONS

#include <set>
#include <vector>
#include <boost/graph/adjacency_list.hpp>

#include "preprocessing.hpp"
#include "lower_bounds.hpp"
#include "elimination_orderings.hpp"
#include "postprocessing.hpp"
#include "dynamicCR.hpp"
#include "exact_cutset.hpp"
#include "separator_algorithm.hpp"
#include "misc.hpp"
#ifdef USE_GALA
#include "exact_ta.hpp"
#endif

namespace treedec{

namespace comb{


template<class G, template<class G_, class ...> class CFGT=algo::default_config>
class PP_MD {
public: // types

    typedef typename treedec::graph_traits<G>::treedec_type T;
public: // construct
    PP_MD(G& g) : _g(g){
    }

public: // random stuff, should be in algo. later
    void set_lower_bound(unsigned lb){
        _low_tw = lb-1;
    }
    unsigned lower_bound()const{
        return _low_tw + 1;
    }

    void do_it(){
        if(boost::num_vertices(_g) == 0){
            boost::add_vertex(_t);
            return;
        }

        // TODO: cleanup
        std::vector<boost::tuple<
            typename treedec_traits<typename treedec_chooser<G>::type>::vd_type,
                     typename treedec_traits<typename treedec_chooser<G>::type>::bag_type
                         > > bags;

        treedec::preprocessing(_g, bags, _low_tw);
        if(boost::num_edges(_g) > 0){
            treedec::minDegree_decomp(
                    _g, _t,
                    (typename std::vector<typename treedec_chooser<G>::value_type>*)NULL,
                    UINT_MAX, true); //ignore_isolated_vertices
        }
        treedec::glue_bags(bags, _t);
    }

    template<class TT>
    void get_tree_decomposition(TT& t) const{
        // todo: assemble td here.
        boost::copy_graph(_t, t);
    }

private:
    G& _g;
    T _t;
    int _low_tw;
};


// TODO: faster.
// TODO: more generic
template<class G_t, template<class G_, class ...> class CFGT=algo::default_config>
class PP_FI{
private:
    typedef typename treedec::graph_traits<G_t>::treedec_type T_t;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;

public: // construct
    PP_FI(G_t& g) : _g(g){
        _low_tw = -1;
    }

public: // random stuff
    void set_lower_bound(unsigned lb){
        _low_tw = lb-1;
    }
    unsigned lower_bound()const{
        return _low_tw + 1;
    }

public: // algo interface
    void do_it(){
        if(boost::num_vertices(_g) == 0){
            boost::add_vertex(_t);
            return;
        }

#ifndef NOBAGS
        std::vector<boost::tuple<
            typename treedec_traits<typename treedec_chooser<G_t>::type>::vd_type,
            typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type
             > > bags;

        treedec::preprocessing(_g, bags, _low_tw);
#else
        impl::preprocessing<G_t> A(_g);
        A.set_treewidth(low_tw, -1u);
        A.do_it();
        low_tw = A.get_treewidth();
        A.get_graph(_g);
#endif

        if(boost::num_edges(_g) > 0){
            unsigned low2=-1;
            treedec::impl::fillIn_decomp(_g, _t, low2, true); //ignore_isolated
            _low_tw = low2;
        }
#ifndef NOBAGS
        treedec::glue_bags(bags, _t);
#else
        skeleton<...> S(...)
        S.do_it();
#endif
    }

    template<class TT>
    void get_tree_decomposition(TT& t) const{
        // todo: assemble td here.
        boost::copy_graph(_t, t);
    }

private:
    G_t &_g;
    T_t _t;
    int _low_tw;
}; // PP_FI


// TODO: faster.
// TODO: more generic
template<class G, template<class G_, class ...> class CFGT=algo::default_config>
class PP_FI_TM{
private:
    typedef typename treedec::graph_traits<G>::treedec_type T;
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;

public: // construct
    PP_FI_TM(G& g) : _g(g){
        _low_tw = -1;
    }

public: // random stuff
    void set_lower_bound(unsigned lb){
        _low_tw = lb-1;
    }
    unsigned lower_bound()const{
        return _low_tw + 1;
    }

public: // algo interface
    void do_it(){

        if(boost::num_vertices(_g) == 0){
            boost::add_vertex(_t);
            return;
        }

        std::vector<boost::tuple<
            typename treedec_traits<typename treedec_chooser<G>::type>::vd_type,
                     typename treedec_traits<typename treedec_chooser<G>::type>::bag_type
                         > > bags;

        treedec::preprocessing(_g, bags, _low_tw);

        if(boost::num_edges(_g) > 0){
            typename std::vector<vertex_descriptor> old_elim_ordering;
            typename std::vector<vertex_descriptor> new_elim_ordering;

            G H(_g);
            //true = ignore isolated vertices
            treedec::fillIn_ordering(_g, old_elim_ordering, true);
            _g = H; // reset

            treedec::minimalChordal(_g, old_elim_ordering, new_elim_ordering);

            typename std::vector<vertex_descriptor>
                new_elim_ordering_(old_elim_ordering.size());

            unsigned c = 0;
            for(unsigned i = 0; i < new_elim_ordering.size(); i++){
                if(boost::out_degree(new_elim_ordering[i], _g) > 0){
                    new_elim_ordering_[c++] = new_elim_ordering[i];
                }
            }

            treedec::ordering_to_treedec(_g, new_elim_ordering_, _t);
        }

        treedec::glue_bags(bags, _t);
    }

    template<class TT>
    void get_tree_decomposition(TT& t) const{
        // todo: assemble td here.
        boost::copy_graph(_t, t);
    }

private:
    G& _g;
    T _t;
    int _low_tw;
}; // PP_FI_TM


template<class G, template<class G_, class ...> class CFGT=algo::default_config>
class FI_TM{
public: // types
    typedef typename treedec::graph_traits<G>::treedec_type T;
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
public: // construct
    FI_TM(G const& g) : _g(g) {
    }

public: // random stuff, should be in algo. later
    void set_lower_bound(unsigned lb){ untested(); //gets bs, works on tw, returns bs?
        _low_tw = lb-1;
    }
    unsigned lower_bound()const{ untested();
        return _low_tw + 1;
    }

public: // algo interface
    void do_it(){
        std::vector<unsigned int> old_elim_ordering;
        typename std::vector<vertex_descriptor> new_elim_ordering;
        treedec::fillIn_ordering(_g, old_elim_ordering);
        treedec::minimalChordal(_g, old_elim_ordering, new_elim_ordering);
        treedec::ordering_to_treedec(_g, new_elim_ordering, _t);
    }

private:
    G& _g;
    T _t;
    int _low_tw;
}; // FI_TM

#ifdef USE_GALA
namespace ex17choice{
  template<class G, template<class G_, class ...> class C=treedec::algo::default_config>
  using exact_ta_=treedec::exact_ta<G, C>;
}

template<class x, template<class G_, class ...> class C=treedec::algo::default_config>
using ex17=treedec::draft::exact_decomposition<x,
             C, ex17choice::exact_ta_>;
#endif

} // comb

//Recursively applies preprocessing rules and glues corresponding bags with
//current tree decomposition this version applies the minDegree-heuristic on
//not fully preprocessable graph instances.
template <typename G_t, typename T_t>
void PP_MD(G_t &G, T_t &T, int &low)
{
    treedec::comb::PP_MD<G_t> a(G);
    a.set_lower_bound(low+1);
    a.do_it();
    low=a.lower_bound()-1;
    a.get_tree_decomposition(T);
}

//Recursively applies preprocessing rules and glues corresponding bags with
//current tree decomposition this version applies the minDegree-heuristic on
//not fully preprocessable graph instances.
template <typename G_t, typename T_t>
void PP_FI(G_t &G, T_t &T, int &low_tw){
    treedec::comb::PP_FI<G_t> a(G);
    a.set_lower_bound(low_tw+1);
    a.do_it();
    low_tw=a.lower_bound()-1;
    a.get_tree_decomposition(T);
}


//Recursively applies preprocessing rules and glues corresponding bags with
//current tree decomposition. This version applies the fillIn-heuristic followed
//by triangulation minimization on not fully preprocessable graph instances.
template <typename G_t, typename T_t>
void PP_FI_TM(G_t &G, T_t &T, int &low){

    comb::PP_FI_TM<G_t> a(G);
    a.set_lower_bound(low+1);
    a.do_it();
    low=a.lower_bound()-1;
    a.get_tree_decomposition(T);
}

template <typename G_t, typename T_t>
void exact_decomposition_cutset(G_t &G, T_t &T, int lb)
{
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    //Preprocessing.
    int low = -1;

    std::vector<boost::tuple<
        typename treedec_traits<typename treedec_chooser<G_t>::type>::vd_type,
        typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type
         > > bags;

    treedec::preprocessing(G, bags, low);

    if(boost::num_edges(G) == 0){
        treedec::glue_bags(bags, T);
        return;
    }

    //Lower bound on the treewidth of the reduced instance of G.
    G_t H(G);
    int lb_deltaC = treedec::lb::deltaC_least_c(H);

    lb = (low > lb)? low : lb;
    lb = (lb_deltaC > lb)? lb_deltaC : lb;

    //Compute a treedecomposition for each connected component of G and glue the decompositions together.
    typedef std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > components_t;
    components_t components;
    treedec::get_components(G, components);

    // root
    boost::add_vertex(T);

    typename components_t::iterator i = components.begin();
    for(; i!=components.end(); ++i){
        //Ignore isolated vertices (already included in 'bags').
        if(i->size() == 1){
            continue;
        }

        G_t G_;
        typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> vdMap;
        treedec::induced_subgraph(G_, G, *i, vdMap);
        T_t T_;

        while(!treedec::exact_cutset(G_, T_, lb)){
            lb++;
        }

        treedec::apply_map_on_treedec(T_, G_, vdMap);

        treedec::glue_decompositions(T, T_);
    }

    treedec::glue_bags(bags, T);
}

//this is in SageMath
template <typename G_t, typename T_t>
bool exact_decomposition_cutset_decision(G_t &G, T_t &T, int k){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);

        if(k >= -1){
            return true;
        }
        else{
            return false;
        }
    }

    //Preprocessing.
    int low_tw = -1;

    std::vector<boost::tuple<
        typename treedec_traits<typename treedec_chooser<G_t>::type>::vd_type,
        typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type
         > > bags;

    treedec::preprocessing(G, bags, low_tw);

    if(boost::num_edges(G) == 0){
        treedec::glue_bags(bags, T);
        if(low_tw <= k){
            return true;
        }
        else{
            return false;
        }
    }

    //Lower bound on the treewidth of the reduced instance of G.
    G_t H(G);
    int lb_deltaC = treedec::lb::deltaC_least_c(H);

    int lb = low_tw;
    if(lb_deltaC > lb){
        lb = lb_deltaC;
    }

    if(lb > k){
        return false;
    }

    //Compute a treedecomposition for each connected component of G and glue the decompositions together.
    std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > components;
    treedec::get_components(G, components);

    // root
    boost::add_vertex(T);

    for(unsigned int i = 0; i < components.size(); i++){
        //Ignore isolated vertices (already included in 'bags').
        if(components[i].size() == 1){
            continue;
        }else{
        }

        G_t G_;
        typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> vdMap;
        treedec::induced_subgraph(G_, G, components[i], vdMap);
        T_t T_;

        while(!treedec::exact_cutset(G_, T_, lb)){
            lb++;
            if(lb > k){
                return false;
            }
        }
    }

    return true;
}

template <typename G_t, typename T_t>
void exact_decomposition_dynamic(G_t &G, T_t &T, int lb){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    //preprocessing
    int low = -1;
    std::vector<boost::tuple<
        typename treedec_traits<typename treedec_chooser<G_t>::type>::vd_type,
        typename treedec_traits<typename treedec_chooser<G_t>::type>::bag_type
         > > bags;

    treedec::preprocessing(G, bags, low);
    if(boost::num_edges(G) == 0){
        treedec::glue_bags(bags, T);
        return;
    }

/*  what is this?!
    if(low > lb){ untested();
    }else{ untested();
    }
*/

    //Compute a treedecomposition for each connected component of G and glue the decompositions together.
    std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > components;
    treedec::get_components(G, components);

    if(components.size() == 1){
        treedec::CR_dynamic_decomp(G, T, lb);

        treedec::glue_bags(bags, T);
        return;
    }

    // root
    boost::add_vertex(T);

    for(unsigned int i = 0; i < components.size(); i++){
        //Ignore isolated vertices (already included in 'bags').
        if(components[i].size() == 1){
            continue;
        }

        G_t G_;
        typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> vdMap;
        treedec::induced_subgraph(G_, G, components[i], vdMap);
        T_t T_;

        treedec::CR_dynamic_decomp(G_, T_, lb);

        treedec::apply_map_on_treedec(T_, G_, vdMap);

        treedec::glue_decompositions(T, T_);
    }

    treedec::glue_bags(bags, T);
}


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
void separator_algorithm_MSVS(G_t &G, T_t &T){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    treedec::separator_algorithm(G, T);
    treedec::MSVS(G, T);
}

template <typename G_t, typename T_t>
void separator_algorithm_TM(G_t &G, T_t &T){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    treedec::separator_algorithm(G, T);
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

template <typename G_t, typename T_t>
void MD_MSVS(G_t &G, T_t &T){
    G_t H(G);
    treedec::minDegree_decomp(G, T);
    treedec::MSVS(H, T);
}

template <typename G_t, typename T_t>
void FI_MSVS(G_t &G, T_t &T){
    int low = -1;
    G_t H(G);
    treedec::fillIn_decomp(G, T);
    treedec::MSVS(H, T);
}

} //namespace treedec

#endif //TD_COMBINATIONS

// vim:ts=8:sw=4:et
