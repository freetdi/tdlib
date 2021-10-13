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
// Foundation, 51 Franklin Street - Suite 500, Boston, MA 02110-1335, USA.
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

#ifndef TREEDEC_COMBINATIONS_HPP
#define TREEDEC_COMBINATIONS_HPP

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
#include "util.hpp"
#ifdef HAVE_GALA_GRAPH_H
# include "exact_ta.hpp" // needs gala/cbset.h
#endif
#include "treedec_copy.hpp"

#include <boost/graph/copy.hpp>
#include <boost/graph/graph_utility.hpp>

namespace treedec{

namespace comb{


template<class G, template<class G_, class ...> class CFGT=algo::default_config>
class PP_MD {
public: // types

    typedef typename treedec::graph_traits<G>::treedec_type T;
public: // construct
    PP_MD(G& g) :
        _g(g),
        _low_tw(-1)
    {
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
        }else{
        }

        // TODO: cleanup
        std::vector<boost::tuple<
            typename treedec_traits<typename treedec_chooser<G>::type>::vd_type,
                     typename treedec_traits<typename treedec_chooser<G>::type>::bag_type
                         > > bags;

        treedec::preprocessing(_g, bags, _low_tw);
        if(boost::num_edges(_g) > 0){

            // HACK. old mindegree does not work on bidirectional graphs
        boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> ug;
        boost::copy_graph(_g, ug,
				 boost::vertex_copy(hack::forgetprop()).
				 edge_copy(hack::forgetprop()));

            treedec::minDegree_decomp(
                    ug, _t,
                    (typename std::vector<typename treedec_chooser<G>::value_type>*)NULL,
                    UINT_MAX, true); //ignore_isolated_vertices
        }else{
        }
        treedec::glue_bags(bags, _t);
    }

    template<class TT>
    void get_tree_decomposition(TT& t) const{
        // todo: assemble td here.
#if 0
        // boost::copy_graph(_t, t); // FIXME
#else
        treedec::obsolete_copy_treedec(_t, t);
#endif
    }
    template<class O>
    void get_elimination_ordering(O&) const{
        incomplete();
    }

private:
    G& _g;
    // T _t; // FIXME. does not work yet
    boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS,
        boost::property<treedec::bag_t, std::set<unsigned> > > _t; // BUG

    int _low_tw;
}; // PP_MD


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
            trace2("calling fillIn", low2, num_vertices(_g));
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
}; // PP_FI_TM (old)


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

#ifdef HAVE_GALA_GRAPH_H
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

namespace pending {

// TODO: faster.
// TODO: more generic
template<class G, template<class G_, class ...> class CFGT=algo::default_config>
class PP_FI : public algo::draft::algo1 {
public:
private:
    typedef typename treedec::graph_traits<G>::treedec_type T;
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
public: // construct
    PP_FI(G& g) :
       algo::draft::algo1("pp_fi"),
       _g(&g),
       _low_tw(-1){
    }
    template<class G_in>
    PP_FI(G_in const& g) :
       algo::draft::algo1("pp_fi"),
       _g(new G),
       _low_tw(-1),
       _own_g(true){
        boost::copy_graph(g, *_g);
    }
    ~PP_FI() {
        if(_own_g){
            delete _g;
        }else{
        }
    }
public: // random stuff
    void set_lower_bound(unsigned lb){
        _low_tw = lb-1;
    }
    unsigned lower_bound()const{
        return _low_tw + 1;
    }
    unsigned bagsize()const{
        // incomplete
        return get_bagsize(_t);
    }
public: // algo interface
    void do_it(){
        // incomplete(); // use comp::
        if(boost::num_vertices(g()) == 0){
            boost::add_vertex(_t);
            return;
        }else{

            // BUG, somehow need to cast CFGT to ppconfig
            // "message" is getting lost here, need pp_cfg+CFGT
            impl::preprocessing<G> A(g());
            A.set_treewidth(_low_tw, -1u);
            A.do_it();

            trace0("doing the rest");

            // this is the rest from exact_base
            A.template do_the_rest<T, treedec::impl::fillIn > (_t);
        }
    }

    template<class TT>
    void get_tree_decomposition(TT& t) const{ untested();
        boost::copy_graph(_t, t);
    }
    template<class TT>
    void get_tree_decomposition(TT& t) {
        boost::copy_graph(_t, t);
    }
    template<class O>
    void get_elimination_ordering(O&) const{
        incomplete();
    }

private:
    G& g(){return *_g;};

private:
    G* _g{nullptr};
    T _t;
    // std::vector<vertex_descriptor> _o;
    int _low_tw;
    bool _own_g{false};
}; // PPFI

struct do_nothing {
    template <typename T1, typename T2>
    void operator()(const T1& , T2& ) const {
    }
};

// pending
template<class G, template<class G_, class ...> class CFGT=algo::default_config>
class PP_FI_TM : public algo::draft::algo1 {
private:
    typedef typename treedec::graph_traits<G>::treedec_type T;
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;

public: // construct
    PP_FI_TM(G& g)
      : algo::draft::algo1("pp_fi_tm"),
       _g(&g),
       _low_tw(-1){
    }
    template<class G_in>
    PP_FI_TM(G_in const& g)
      : algo::draft::algo1("pp_fi_tm"),
       _g(new G),
       _own_g(true),
       _low_tw(-1) {
        boost::copy_graph(g, *_g);
    }
    ~PP_FI_TM() {
        if(_own_g){
            delete _g;
        }else{
        }
    }

public: // random stuff
    void set_lower_bound(unsigned lb){ untested();
        _low_tw = lb-1;
    }
    unsigned lower_bound()const{ untested();
        return _low_tw + 1;
    }
    unsigned bagsize()const{
        return get_bagsize(_t);
    }
    template<class O>
    void get_elimination_ordering(O&) const{
        incomplete();
    }

public: // algo interface
    void do_it(){

        if(boost::num_vertices(g()) == 0){
            boost::add_vertex(_t);
            return;
        }

        std::vector<boost::tuple<
            typename treedec_traits<typename treedec_chooser<G>::type>::vd_type,
                     std::vector<vertex_descriptor> > > bags;

#if 0
        treedec::preprocessing(g(), bags, _low_tw);
#else
            impl::preprocessing<G> A(g());
            A.set_treewidth(_low_tw, -1u);
            A.do_it();
            A.get_bags(bags);
            A.get_graph(g());

            trace1("done PP in PPFITM", boost::num_edges(g()));
#endif

        for(auto const& x : bags){
            auto& B=boost::get<1>(x);
            trace1("B", B.size());
        }

        if(boost::num_edges(g()) > 0){
//            typename std::vector<vertex_descriptor> old_elim_ordering;
            typename std::vector<vertex_descriptor> new_elim_ordering;

            // BUG this must work in ordering_to_treedec
            boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS> H;
            boost::copy_graph(g(), H,
				 boost::vertex_copy(hack::forgetprop()).
				 edge_copy(hack::forgetprop()));

            trace0("backedup to H");
            trace0("doit");
            treedec::impl::fillIn<G> a(g());
            a.set_ignore_isolated();
            a.do_it();
            trace0("done fi");
            auto& old_elim_ordering = a.get_elimination_ordering();

            trace2("mC", boost::num_edges(H), boost::num_vertices(g()));

            // changes H, but doesn't matter
            treedec::minimalChordal(H, old_elim_ordering, new_elim_ordering);
            trace3("", old_elim_ordering.size(), new_elim_ordering.size(), boost::num_vertices(g()));
            assert(is_vertex_permutation(new_elim_ordering, g()));


            typename std::vector<vertex_descriptor>
            
            new_elim_ordering_(boost::num_vertices(H));

            unsigned c = 0;
            for(auto n : new_elim_ordering) {
                trace4("eo", c, n, degree(n, g()), boost::out_degree(n, g()));
                if(degree(n, g()) > 0){
                    assert(c<new_elim_ordering_.size());
                    new_elim_ordering_[c++] = n;
                }else{
                }
            }
            new_elim_ordering_.resize(c);

            assert(is_vertex_permutation(new_elim_ordering, H));

            assert(is_permutation(new_elim_ordering));

            trace3("to_treedec", c, new_elim_ordering.size(), boost::num_vertices(H));
            trace3("", boost::num_edges(H), new_elim_ordering.size(), boost::num_vertices(H));
            treedec::ordering_to_treedec(H, new_elim_ordering_, _t);
            trace0("PPFITM ordered_to_treedec");
            // boost::print_graph(_t);
        }

        trace0("gluing");
        treedec::glue_bags(bags, _t);
    }

    template<class TT>
    void get_tree_decomposition(TT& t){
        // incomplete();
        treedec::obsolete_copy_treedec(_t, t);
    }
    template<class TT>
    void get_tree_decomposition(TT& t) const{ untested();
        // todo: assemble td here.
#if 0
        // boost::copy_graph(_t, t); // FIXME
#else
        treedec::obsolete_copy_treedec(_t, t);
#endif
    }

private:
    G& g() { assert(_g); return *_g; }
private:
    G* _g{nullptr};
    bool _own_g{false};
    // T _t; // FIXME. does not work yet

    // TD_tree_dec_t _t; // BUG
    // needs to work with ordering_to_treedec... (which?!)
    boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
        boost::property<treedec::bag_t, std::set<unsigned> > > _t; // BUG

    int _low_tw;
}; // PP_FI_TM

} // pending

//Recursively applies preprocessing rules and glues corresponding bags with
//current tree decomposition this version applies the minDegree-heuristic on
//not fully preprocessable graph instances.
template <typename G_t, typename T_t>
void PP_FI(G_t &G, T_t &T, int &low_tw)
{
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
void PP_FI_TM(G_t &G, T_t &T, int &low)
{
    comb::PP_FI_TM<G_t> a(G);
    a.set_lower_bound(low+1);
    a.do_it();
    low=a.lower_bound()-1;
    a.get_tree_decomposition(T);
}

#ifdef HAVE_GALA_GRAPH_H
template <typename G_t, typename T_t>
void exact_decomposition_ex17(G_t &G, T_t &T, int lb_tw)
{
    using draft::exact_decomposition;
    auto alg=exact_decomposition<G_t, algo::default_config, exact_ta>(G);
                                   // really^?
    return alg.try_it(T, lb_tw+1);
}
#endif

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

#endif //TREEDEC_COMBINATIONS_HPP

// vim:ts=8:sw=4:et
