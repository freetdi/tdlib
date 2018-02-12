// Lukas Larisch, 2014 - 2017
// Felix Salfelder, 2016 - 2017
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
 * Offers functionality to compute lower bounds on the tree width of a graph.
 *
 * Provides following functions (namespace treedec::lb):
 *
 * - int delta(G_t &G)
 * - int delta2(G_t &G)
 * - int gamma(G_t &G)
 * - int deltaD(G_t G)
 * - int delta2D(G_t &G)
 * - int gammaD_left(G_t &G)
 * - int gammaD_right(G_t &G)
 * - int gammaD_min_e(G_t &G)
 * - int deltaC_min_d(G_t &G)
 * - int deltaC_max_d(G_t &G)
 * - int deltaC_least_c(G_t &G)
 *
 * - void k_neighbour_improved_graph(G_t &G, unsigned int k)
 * - int LBN_deltaD(G_t &G)
 * - int LBN_deltaC(G_t &G)
 * - int LBNC_deltaD(G_t &G)
 * - int LBNC_deltaC(G_t &G)
 * - void k_path_improved_graph(G_t &G, unsigned int k)
 * - int LBP_deltaD(G_t &G)
 * - int LBP_deltaC(G_t &G)
 * - int LBPC_deltaD(G_t &G)
 * - int LBPC_deltaC(G_t &G)
 *
 * - int MCS(G_t &G)
 * - int MCSC(G_t &G)
 *
 * - int relation_edges_vertices(G_t &G)
 *
 * The main idea of all algorithms included in this file is this one:
 *
 *    For any tree decomposition T of width k, there is a 'small tree decomposition' T' (no bag is subset of another bag) of width k.
 *    If a graph G is a complete graph, than all trees of tree decompositions of G consists of an isolated vertex t with B(t) = V(G).
 *    In this case, the minimal degree of a vertex in G matches the treewidth of G. If G is not a complete graph, than all trees of
 *    tree decompositions of G of minimal width have at least two vertices and hence at least two leafs. Let t be a leaf of a small tree
 *    decomposition T' of G of minimal width and let t' be adjacent with t in T'. Then B(t) \ B(t') is not empty and all neighbours
 *    of a vertex v in (B(t) \ B(t')) are contained in B(t). Let d be the degree of v in G. The minimal degree d_min of a vertex in G is
 *    a lower bound of the treewidth of G, since d_min <= d <= tw(G) (algorithm delta). The same argument holds for the second minimal
 *    degree of vertices in in G (algorithm delta2) and some variation, which is introduced in gamma. The minimal degree-method also
 *    holds for subgraphs and minors of G. The degeneracy of G is the maximum over all smallest degrees of vertices in the graphs
 *    obtained by successivly removing a vertex of minimal degree. The algorithms ...D apply the algorithms delta, delta2 and gamma
 *    on all (degeneracy-)subgraphs of G. The algorithms ...C apply the algorithms delta, delta2 and gamma on some heuristically choosen
 *    minors of G.
 *
 * For more information, see e.g.:
 *
 *     Hans L. Bodlaender, Arie M.C.A. Koster, Treewidth computations II. Lower bounds, Information and Computation,
 *     Volume 209, Issue 7, July 2011, Pages 1103-1119
 *
 */

#ifndef TREEDEC_LOWER_BOUNDS_HPP
#define TREEDEC_LOWER_BOUNDS_HPP

#include <set>
#include <vector>
#include <climits>

#include <boost/graph/adjacency_list.hpp>
#include <utility>

#include "degree.hpp"
#include "graph.hpp"
#include "misc.hpp"
#include "network_flow.hpp"
#include "simple_graph_algos.hpp"
#include "algo.hpp"
#include "trace.hpp"
#include "degree_config.hpp"

#include "impl/greedy_base.hpp"

namespace treedec{

namespace lb{


/* DEGREE BASED */

namespace impl{

//Smallest vertex-degree in G.
template <typename G_t>
class delta : public treedec::algo::draft::algo1{
public:
    delta(const G_t &G) : algo1("lb::delta"), _g(G), _lb(0){}

    void do_it(){ untested();
        timer_on();

        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(_g); vIt != vEnd; vIt++){ untested();
            unsigned degree = boost::out_degree(*vIt, _g);
            _lb = (degree < _lb)? degree : _lb;
        }

        timer_off();
    }

    unsigned lower_bound_bagsize(){ untested();
        return _lb+1u;
    }

private:
    const G_t& _g;
    unsigned _lb;
};

} //namespace impl


//Smallest vertex-degree in G.
template <typename G_t>
int delta(const G_t &G){ untested();
    if(boost::num_vertices(G) == 0){ untested();
        return -1;
    }else{ untested();
    }

    impl::delta<G_t> delta(G);
    delta.do_it();
    return (int)delta.lower_bound_bagsize()-1;
}


namespace impl{

//Second smallest vertex-degree in G.
template <typename G_t>
class delta2 : public treedec::algo::draft::algo1{
public:
    delta2(const G_t &G) : algo1("lb::delta2"), _g(G), _lb(0){}

    void do_it(){ untested();
        timer_on();

        unsigned min = boost::num_vertices(_g);
        unsigned snd = min;

        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(_g); vIt != vEnd; vIt++){ untested();
            unsigned int degree = boost::out_degree(*vIt, _g);
            if(degree <= min){ untested();
                snd = min;
                min = degree;
            }else{ untested();
            }
            if(degree > min && degree < snd){ untested();
                snd = degree;
            }else{ untested();
            }
        }

        _lb = snd;

        timer_off();
    }

    unsigned lower_bound_bagsize(){ untested();
        return _lb+1u;
    }

private:
    const G_t& _g;
    unsigned _lb;
};

} //namespace impl


//Second smallest vertex-degree in G.
template <typename G_t>
int delta2(const G_t &G){ untested();
    if(boost::num_vertices(G) == 0){ untested();
        return -1;
    }else if(boost::num_vertices(G) == 1){ untested();
        return 0;
    }else{ untested();
    }

    impl::delta2<G_t> delta2(G);
    delta2.do_it();
    return (int)delta2.lower_bound_bagsize()-1;
}


namespace impl{

template <typename G_t>
class gamma : public treedec::algo::draft::algo1{
public:
    gamma(const G_t &G) : algo1("lb::gamma"), _g(G), _lb(0){}

    void do_it(){ untested();
        timer_on();

        //Sort the vertices of G according to rising degree -> degree_sequence.
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> degree_sequence;
        make_degree_sequence(_g, degree_sequence);

        //Take the degree of the right vertex in the first not-edge.
        for(unsigned int i = 0; i < boost::num_vertices(_g); i++){ untested();
            for(unsigned int j = 0; j < i; j++){ untested();
                if(!boost::edge(degree_sequence[i], degree_sequence[j], _g).second){ untested();
                    _lb = boost::out_degree(degree_sequence[i], _g);
                    timer_off();
                    return;
                }else{ untested();
                }
            }
        }
        unreachable();
    }

    unsigned lower_bound_bagsize(){ untested();
        return _lb+1u;
    }

private:
    const G_t& _g;
    unsigned _lb;
};

} //namespace impl


template <typename G_t>
int gamma(const G_t &G)
{ untested();
    unsigned int V = boost::num_vertices(G);
    unsigned int E = boost::num_edges(G);

    if(V == 0){ untested();
        return -1;
    }else if(E == 0){ untested();
        return 0;
    }else if(2*E+1 == V*(V-1u)){ untested();
        return V-1;
    }else{ untested();
    }

    impl::gamma<G_t> gamma(G);
    gamma.do_it();
    return (int)gamma.lower_bound_bagsize()-1;
}


namespace impl{

template <typename G_t>
class deltaD : public treedec::algo::draft::algo1{
public:
    deltaD(G_t &G) : algo1("lb::deltaD"), _g(G), _lb(0){}

    void do_it(){ untested();
        timer_on();

        unsigned int maxmin = 0;
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex =
                                                *(boost::vertices(_g).second);

        while(true){ untested();
            unsigned int min_degree = boost::num_vertices(_g);
            for(boost::tie(vIt, vEnd) = boost::vertices(_g); vIt != vEnd; vIt++){ untested();
                unsigned int degree = boost::out_degree(*vIt, _g);
                if(degree < min_degree && degree > 0){ untested();
                    min_degree = degree;
                    min_vertex = *vIt;
                }
            }

            if(min_degree == boost::num_vertices(_g)){ untested();
                timer_off();
                _lb = maxmin;
            }else{ untested();
            }

            maxmin = (maxmin>min_degree)? maxmin : min_degree;

            boost::clear_vertex(min_vertex, _g);
        }
    }

    unsigned lower_bound_bagsize(){ untested();
        return _lb+1u;
    }

private:
    G_t& _g;
    unsigned _lb;
};

} //namespace impl


template <typename G_t>
int deltaD(G_t& G)
{ untested();
    if(boost::num_vertices(G) == 0){ untested();
        return -1;
    }else{ untested();
    }

    impl::deltaD<G_t> deltaD(G);
    deltaD.do_it();
    return (int)deltaD.lower_bound_bagsize()-1;
}

namespace impl{

template <typename G_t>
class delta2D : public treedec::algo::draft::algo1{
public:
    delta2D(const G_t &G) : algo1("lb::delta2D"), _g(G), _lb(0){}

    void do_it(){ untested();
        timer_on();

        G_t H;
        boost::copy_graph(_g, H);

        typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> assumed_minimal;

        typename boost::graph_traits<G_t>::vertex_iterator hIt, hEnd;
        for(boost::tie(hIt, hEnd) = boost::vertices(H); hIt != hEnd; hIt++){ untested();
            assumed_minimal.push_back(*hIt);
        }

        unsigned int min_degree, maxmin = 0;
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex = *(boost::vertices(H).first);

        for(unsigned int i = 0; i < assumed_minimal.size(); i++){ untested();
            while(boost::num_edges(H) > 0){ untested();
                min_degree = boost::num_vertices(H);

                for(boost::tie(hIt, hEnd) = boost::vertices(H); hIt != hEnd; hIt++){ untested();
                    if(*hIt == assumed_minimal[i]){ untested();
                        continue;
                    }else{ untested();
                    }
                    unsigned int degree = boost::out_degree(*hIt, H);
                    if(degree < min_degree && degree > 0){ untested();
                        min_degree = degree;
                        min_vertex = *hIt;
                    }else{ untested();
                    }
                }
                if(min_degree == boost::num_vertices(H)){ untested();
                    break;
                }else{ untested();
                }

                maxmin = (maxmin>min_degree)? maxmin : min_degree;
                boost::clear_vertex(min_vertex, H);
            }
            H.clear();
            boost::copy_graph(_g, H);
        }

        _lb = maxmin;

        timer_off();
    }

    unsigned lower_bound_bagsize(){ untested();
        return _lb+1u;
    }

private:
    const G_t& _g;
    unsigned _lb;
};

} //namespace impl

//Assume each vertex as the minimal one and apply deltaD.
template <typename G_t>
int delta2D(const G_t &G){ untested();
    if(boost::num_vertices(G) == 0){ untested();
        return -1;
    }else{ untested();
    }

    impl::delta2D<G_t> delta2D(G);
    delta2D.do_it();
    return (int)delta2D.lower_bound_bagsize()-1;
}

namespace impl{

template <typename G_t>
class gammaD_left : public treedec::algo::draft::algo1{
public:
    gammaD_left(G_t &G) : algo1("lb::gammaD_left"), _g(G), _lb(0){}

    void gammaD_left_recursion(){ untested();
        if(boost::num_edges(_g) == 0){ untested();
            return;
        }else{ untested();
        }

        //Sort the vertices of G according to rising degree -> degree_sequence.
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> degree_sequence;
        make_degree_sequence(_g, degree_sequence);

        for(unsigned int i = 0; i < degree_sequence.size(); i++){ untested();
            for(unsigned int j = 0; j < i; j++){ untested();
                if(boost::edge(degree_sequence[i], degree_sequence[j], _g).second){ untested();
                    continue;
                }else{ untested();
                }

                //gammaD-left heuristic
                unsigned int degree = boost::out_degree(degree_sequence[i], _g);
                for(unsigned int k = 0; k < i; k++){ untested();
                    boost::clear_vertex(degree_sequence[k], _g);
                }

                _lb = (degree > _lb)? degree : _lb;

                gammaD_left_recursion();
                return;
            }
        }
    }

    void do_it(){ untested();
        timer_on();
        gammaD_left_recursion();
        timer_off();
    }

    unsigned lower_bound_bagsize(){ untested();
        return _lb+1u;
    }

private:
    G_t& _g;
    unsigned _lb;
};

} //namespace impl


template <typename G_t>
int gammaD_left(G_t& G)
{ untested();
    unsigned int V = boost::num_vertices(G);
    unsigned int E = boost::num_edges(G);

    if(V == 0){ untested();
        return -1;
    }else if(E == 0){ untested();
        return 0;
    }else if(2*E == V*(V-1u)){ untested();
        return V-1;
    }else{ untested();
        impl::gammaD_left<G_t> gammaD_left(G);
        gammaD_left.do_it();
        return (int)gammaD_left.lower_bound_bagsize()-1;
    }
}

namespace impl{

template <typename G_t>
class gammaD_right : public treedec::algo::draft::algo1{
public:
    gammaD_right(G_t &G) : algo1("lb::gammaD_right"), _g(G), _lb(0){}

        void gammaD_right_recursion(){ untested();
        if(boost::num_edges(_g) == 0){ untested();
            return;
        }else{ untested();
        }

        //Sort the vertices of G according to rising degree -> degree_sequence.
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> degree_sequence;
        make_degree_sequence(_g, degree_sequence);

        for(unsigned int i = 0; i < degree_sequence.size(); i++){ untested();
            for(unsigned int j = 0; j < i; j++){ untested();
                if(boost::edge(degree_sequence[i], degree_sequence[j], _g).second){ untested();
                    continue;
                }

                //gammaD-right heuristic
                unsigned int degree = boost::out_degree(degree_sequence[i], _g);
                boost::clear_vertex(degree_sequence[i], _g);

                _lb = (degree > _lb)? degree : _lb;

                gammaD_right_recursion();
                return;
            }
        }
    }

    void do_it(){ untested();
        timer_on();
        gammaD_right_recursion();
        timer_off();
    }

    unsigned lower_bound_bagsize(){ untested();
        return _lb+1u;
    }

private:
    G_t& _g;
    unsigned _lb;
};

} //namespace impl


template <typename G_t>
int gammaD_right(G_t& G)
{ untested();
    unsigned int V = boost::num_vertices(G);
    unsigned int E = boost::num_edges(G);

    if(V == 0){ untested();
        return -1;
    }
    else if(E == 0){ untested();
        return 0;
    }
    else if(2*E == V*(V-1u)){ untested();
        return V-1;
    }
    else{ untested();
        impl::gammaD_right<G_t> gammaD_right(G);
        gammaD_right.do_it();
        return (int)gammaD_right.lower_bound_bagsize()-1;
    }
}

namespace impl{

template <typename G_t>
class gammaD_min_e : public treedec::algo::draft::algo1{
public:
    gammaD_min_e(G_t &G) : algo1("lb::gammaD_min_e"), _g(G), _lb(0){}

    void gammaD_min_e_recursion(){ untested();
        if(boost::num_edges(_g) == 0){ untested();
            return;
        }else{ untested();
        }

        //Sort the vertices of G according to rising degree -> degree_sequence.
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> degree_sequence;
        make_degree_sequence(_g, degree_sequence);

        for(unsigned int i = 0; i < degree_sequence.size(); i++){ untested();
            for(unsigned int j = 0; j < i; j++){ untested();
                if(boost::edge(degree_sequence[i], degree_sequence[j], _g).second){ untested();
                    continue;
                }else{ untested();
                }

                //gammaD-min-e heuristic
                unsigned int degree_right = boost::out_degree(degree_sequence[i], _g);
                unsigned int degree_left = 0;
                for(unsigned int k = 0; k < i; k++){ untested();
                    degree_left += boost::out_degree(degree_sequence[k], _g);
                }

                if(degree_left < degree_right){ untested();
                    for(unsigned int t = 0; t < i; t++){ untested();
                        boost::clear_vertex(degree_sequence[t], _g);
                    }
                }else{ untested();
                    boost::clear_vertex(degree_sequence[i], _g);
                }

                _lb = (degree_right > _lb)? degree_right : _lb;

                gammaD_min_e_recursion();
                return;
            }
        }
    }

    void do_it(){ untested();
        timer_on();
        gammaD_min_e_recursion();
        timer_off();
    }

    unsigned lower_bound_bagsize(){ untested();
        return _lb+1u;
    }

private:
    G_t& _g;
    unsigned _lb;
};

} //namespace impl


template <typename G_t>
int gammaD_min_e(G_t& G)
{ untested();
    unsigned int V = boost::num_vertices(G);
    unsigned int E = boost::num_edges(G);

    if(V == 0){ untested();
        return -1;
    }
    else if(E == 0){ untested();
        return 0;
    }
    else if(2*E+1 == V*(V-1u)){ untested();
        return V-1;
    }
    else{ untested();
        impl::gammaD_min_e<G_t> gammaD_min_e(G);
        gammaD_min_e.do_it();
        return (int)gammaD_min_e.lower_bound_bagsize()-1;
    }
}


namespace impl{

template <typename G_t>
class deltaC_min_d : public treedec::algo::draft::algo1{
public:
    deltaC_min_d(G_t &G) : algo1("lb::deltaC_min_d"), _g(G), _lb(0){}

    void do_it(){
        timer_on();

        while(boost::num_edges(_g) > 0){
            //Search a minimum-degree-vertex.
            typename boost::graph_traits<G_t>::vertex_descriptor min_vertex
                          = get_min_degree_vertex(_g, true); //ignore isolated vertices

            _lb = (_lb>boost::out_degree(min_vertex, _g))? _lb : boost::out_degree(min_vertex, _g);

            //min_d heuristic: Search a neighbour of min_vertex with minimal degree.
            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            unsigned int min_degree_w = boost::num_vertices(_g);

            boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, _g);
            typename boost::graph_traits<G_t>::vertex_descriptor w = *nIt;
            for(; nIt != nEnd; nIt++){
                unsigned int degree = boost::out_degree(*nIt, _g);
                if(degree <= min_degree_w){
                    min_degree_w = degree;
                    w = *nIt;
                }
            }

            contract_edge(min_vertex, w, _g);
        }

        timer_off();
    }

    unsigned lower_bound_bagsize(){
        return _lb+1u;
    }

private:
    G_t& _g;
    unsigned _lb;
};

} //namespace impl

template <typename G_t>
int deltaC_min_d(G_t& G)
{
    unsigned int V = boost::num_vertices(G);
    unsigned int E = boost::num_edges(G);

    if(V == 0){
        return -1;
    }
    else if(E == 0){
        return 0;
    }
    else if(2*E+1 == V*(V-1u)){ untested();
        return V-1;
    }
    else{
        impl::deltaC_min_d<G_t> deltaC_min_d(G);
        deltaC_min_d.do_it();
        return (int)deltaC_min_d.lower_bound_bagsize()-1;
    }
}


namespace impl{

template <typename G_t>
class deltaC_max_d : public treedec::algo::draft::algo1{
public:
    deltaC_max_d(G_t &G) : algo1("lb::deltaC_max_d"), _g(G), _lb(0){}

    void do_it(){
        timer_on();

        while(boost::num_edges(_g) > 0){
            //Search a minimum-degree-vertex.
            typename boost::graph_traits<G_t>::vertex_descriptor min_vertex
                          = get_min_degree_vertex(_g, true); //ignore isolated vertices

            _lb = (_lb>boost::out_degree(min_vertex, _g))? _lb : boost::out_degree(min_vertex, _g);

            //min_d heuristic: Search a neighbour of min_vertex with minimal degree.
            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            unsigned int max_degree = 0;

            boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, _g);
            typename boost::graph_traits<G_t>::vertex_descriptor w = *nIt;
            for(; nIt != nEnd; nIt++){
                unsigned int degree = boost::out_degree(*nIt, _g);
                if(degree > max_degree){
                    max_degree = degree;
                    w = *nIt;
                }
            }

            contract_edge(min_vertex, w, _g);
        }

        timer_off();
    }

    unsigned lower_bound_bagsize(){
        return _lb+1u;
    }

private:
    G_t& _g;
    unsigned _lb;
};

} //namespace impl


template <typename G_t>
int deltaC_max_d(G_t& G)
{
    unsigned int V = boost::num_vertices(G);
    unsigned int E = boost::num_edges(G);

    if(V == 0){
        return -1;
    }
    else if(E == 0){
        return 0;
    }
    else if(2*E == V*(V-1u)){
        return V-1u;
    }
    else{
        impl::deltaC_max_d<G_t> deltaC_max_d(G);
        deltaC_max_d.do_it();
        return (int)deltaC_max_d.lower_bound_bagsize()-1;
    }
}

template<typename G_t>
struct degree_decrease
   : public vertex_callback<typename boost::graph_traits<G_t>::vertex_descriptor>{ //
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename misc::DEGS<G_t>::bag_type degbag;
    typedef typename deg_chooser<G_t>::type degs_type;

    degree_decrease(degs_type* d, G_t*g) :
        _degs(d), G(g){}

    void operator()(vertex_descriptor v){ untested();
        size_t degree = boost::out_degree(v, *G);
        if(degree==0){ untested();
            // unreachable
            // unconnected nodes are unreachable throught adj iterator
            assert(false);
        }else if(degree==1){ untested();
            // unreachable
            // a degree one node does not change its degree during collapse
            assert(false);
        }else{ untested();
            _degs->unlink(v);
            _degs->reg(v, degree-1);
        }
    }
private:
    degs_type*_degs;
    G_t* G;
};

namespace impl{

// TODO: move to impl, then use from lb/order
template <typename G_t,
          template<class G, class...> class CFGT=algo::default_config>
class deltaC_least_c
  : private treedec::impl::greedy_base<
                 G_t,
                 std::vector< typename boost::graph_traits<G_t>::vertex_descriptor >,
                 CFGT>
{ //
public:
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
private:
    typedef treedec::impl::greedy_base<
                 G_t, std::vector<vertex_descriptor>, CFGT> baseclass;

    using typename baseclass::graph_type;
    // typedef typename deg_chooser<typename baseclass::graph_type>::type degs_type;
    typedef DEGS<graph_type, ::treedec::degs::mapped_config> degs_type;

    using typename baseclass::marker_type;
    using baseclass::_g;
    using baseclass::_marker;
    using baseclass::_numbering;
    using baseclass::_degreemap;
    using baseclass::_subgraph;
    using baseclass::_num_edges;
public:

    deltaC_least_c(G_t &G)
      : baseclass(G, -1u),
//        _g(G), // baseclass _g?
        _lb_tw(0)
    {
        trace2("deltacleastc", boost::num_vertices(G), boost::num_edges(G));
    }

    void do_it(){

#ifndef NDEBUG
        auto p=boost::vertices(_g);
        for(; p.first!=p.second; ++p.first){
            assert(_degreemap[*p.first] == boost::out_degree(*p.first, _g));
        }
#endif

        degs_type degs(_g, baseclass::_degreemap);
//        degree_decrease<graph_type> cb(&degs, &_g);

        unsigned int min_ntd = 2;

        while(_num_edges){
            //Search a minimum-degree-vertex.
            if(min_ntd>1){
                --min_ntd;
            }else{
            }

            std::pair<vertex_descriptor, vertices_size_type> min_pair;
            min_pair = degs.pick_min(min_ntd);
            min_ntd = min_pair.second;
            trace2("dclc", min_pair.first, min_ntd);

            if(_lb_tw < min_ntd){
                _lb_tw = min_ntd;
            }else{
            }

            vertex_descriptor min_vertex;
            min_vertex = min_pair.first;
            assert(_degreemap[min_vertex]);
            assert(_degreemap[min_vertex]==min_pair.second);

            //least-c heuristic: search the neighbour of min_vertex such that
            //contracting {min_vertex, w} removes the least edges
            vertex_descriptor w = get_least_common_vertex(
                    min_vertex, _marker, _subgraph);

            //Contract the edge between min_vertex into w.
            //Clear min_vertex and rearrange degs through callback.
            contract_edge(min_vertex, w, degs);
        }
    }

    unsigned lower_bound_bagsize(){
        return _lb_tw+1u;
    }
private:
    // sort of "eliminate v".
    template<class D> // TODO...
    void contract_edge(vertex_descriptor v,
                       vertex_descriptor target, D& _degs) {
        _numbering.put(v);
        _degs.unlink(v);

        trace3("elim", v, _num_edges, _degreemap[v]);
//        _numbering.increment();
#ifndef NDEBUG
        // precondition: the neighs of v are marked!
        {
            auto p=adjacent_vertices(v, _subgraph);
            for(; p.first!=p.second; ++p.first){
                assert(_marker.is_marked(*p.first));
            }
        }
#endif
        auto p=adjacent_vertices(target, _subgraph);
        for(; p.first!=p.second; ++p.first){
            _marker.unmark(*p.first);
            assert(*p.first!=v);
        }
        _marker.unmark(target);

        // TODO/LATER: modify baseclass::subgraph ...
        auto q=boost::adjacent_vertices(v, _subgraph);
        for(; q.first!=q.second; ++q.first){
            trace3("n", v, target, *q.first);
            if(_marker.is_marked(*q.first)){
                assert(*q.first!=target);
            }else{
            }

            if(*q.first==target){
                --_num_edges;
                assert(_degreemap[*q.first]);
                --_degreemap[*q.first];
            }else if(_marker.is_marked(*q.first)){
                // a neigh of v not connected to target.
                // "move" edge.
                assert(!boost::edge(target, *q.first, _g).second);
                assert(!boost::edge(*q.first, target, _g).second);
                treedec::add_edge(target, *q.first, _g);
                assert(boost::edge(target, *q.first, _g).second);
                assert(boost::edge(*q.first, target, _g).second);
                // ++_degreemap[*q.first];
                ++_degreemap[target];
                // ++_num_edges;
            }else{
                // this one has been connected to both. now only one.
                // tell degs...
                --_num_edges;
                assert(_degreemap[*q.first]);
                --_degreemap[*q.first];
                assert(_degreemap[*q.first]);
                _degs.update(*q.first);
            }


        }
        _degs.update(target);

//        boost::clear_vertex(v, g);
    }

private: // overrides
    bool next(vertex_descriptor &) { incomplete(); return false;}
    void eliminate(vertex_descriptor) { incomplete(); }
private:
    unsigned _lb_tw;
};

} //namespace impl



template <typename G_t>
int deltaC_least_c(G_t& G)
{
    auto V=boost::num_vertices(G);
    auto E=treedec::num_edges(G);

    if(V == 0){
        return -1;
    }
    else if(E == 0){
        return 0;
    }
    else if(2*E == V*(V-1u)){
        return V-1u;
    }
    else{
        impl::deltaC_least_c<G_t> deltaC_least_c(G);
        deltaC_least_c.do_it();
        return (int)deltaC_least_c.lower_bound_bagsize()-1;
    }
}


/* IMPROVED GRAPHS */

template <typename G_t>
void k_neighbour_improved_graph(G_t &G, unsigned int k){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> edges_to_add;

    typename boost::graph_traits<G_t>::vertex_iterator vIt1, vIt2, vEnd;
    for(boost::tie(vIt1, vEnd) = boost::vertices(G); vIt1 != vEnd; vIt1++){
        vIt2 = vIt1;
        vIt2++;
        for(; vIt2 != vEnd; vIt2++){
            if(!boost::edge(*vIt1, *vIt2, G).second){
                std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N1, N2;
                typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt1, G); nIt != nEnd; nIt++){
                    N1.insert(*nIt);
                }
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt2, G); nIt != nEnd; nIt++){
                    N2.insert(*nIt);
                }
                std::set<typename boost::graph_traits<G_t>::vertex_descriptor> intersection;

                std::set_intersection(N1.begin(), N1.end(), N2.begin(), N2.end(),
                                      std::inserter(intersection, intersection.begin()));

                if(intersection.size() >= k){
                    edges_to_add.push_back(*vIt1);
                    edges_to_add.push_back(*vIt2);
                }
            }
        }
    }

    for(unsigned int i = 0; i < edges_to_add.size(); ){
        boost::add_edge(edges_to_add[i], edges_to_add[i+1], G);
        ++i;
        ++i;
    }
}


namespace impl{

template <typename G_t, typename CFG_t>
class LB_improved_base : public treedec::algo::draft::algo1{
public:
    LB_improved_base(G_t &G) : algo1(CFG_t::name()), _g(G), _lb(0){}

    void do_it(){
        timer_on();

        G_t H(_g);
        int lb = CFG_t::lb_algo(H);

        while(true){
            G_t H;
            boost::copy_graph(_g, H);
            CFG_t::improvement_algo(H, lb+1);

            int new_lb = CFG_t::lb_algo(H);
            if(new_lb > lb){
                lb++;
            }else{
                break;
            }
        }

        _lb = lb;

        timer_off();
    }


    unsigned lower_bound_bagsize(){
        return _lb+1u;
    }

private:
    G_t& _g;
    unsigned _lb;
};

} //namespace impl



template <typename G_t>
struct CFG_LBN_deltaD{
    static int lb_algo(G_t &H){ untested();
        impl::deltaD<G_t> deltaD(H);
        deltaD.do_it();
        return (int)deltaD.lower_bound_bagsize() - 1;
    }

    static void improvement_algo(G_t &H, unsigned k){ untested();
        treedec::lb::k_neighbour_improved_graph(H, k);
    }

    static const std::string name(){ untested();
        return "lb::LBN_deltaD";
    }
};

template <typename G_t>
int LBN_deltaD(G_t &G)
{ untested();
    unsigned int V = boost::num_vertices(G);
    unsigned int E = boost::num_edges(G);

    if(V == 0){ untested();
        return -1;
    }else if(E == 0){ untested();
        return 0;
    }else if(2*E == V*(V-1u)){ untested();
        return V-1u;
    }else{ untested();
    }

    impl::LB_improved_base<G_t, CFG_LBN_deltaD<G_t> > LBN_deltaD(G);
    LBN_deltaD.do_it();
    return (int)LBN_deltaD.lower_bound_bagsize()-1;
}

template <typename G_t>
int LBN_deltaD(G_t const&G)
{ untested();
    G_t H(G);
    return LBN_deltaD(H);
}

template <typename G_t>
struct CFG_LBN_deltaC{
    static int lb_algo(G_t &H){
        impl::deltaC_least_c<G_t> deltaC(H);
        deltaC.do_it();
        return (int)deltaC.lower_bound_bagsize() - 1;
    }

    static void improvement_algo(G_t &H, unsigned k){
        treedec::lb::k_neighbour_improved_graph(H, k);
    }

    static const std::string name(){
        return "lb::LBN_deltaC";
    }
};

template <typename G_t>
int LBN_deltaC(G_t &G){
    unsigned int V = boost::num_vertices(G);
    unsigned int E = boost::num_edges(G);

    if(V == 0){
        return -1;
    }
    else if(E == 0){
        return 0;
    }
    else if(2*E == V*(V-1u)){
        return V-1u;
    }

    impl::LB_improved_base<G_t, CFG_LBN_deltaC<G_t> > LBN_deltaC(G);
    LBN_deltaC.do_it();
    return (int)LBN_deltaC.lower_bound_bagsize()-1;
}


namespace impl{

template <typename G_t, typename CFG_t>
class LB_improved_contraction_base : public treedec::algo::draft::algo1{
public:
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
private:
    typedef treedec::draft::sMARKER<vertices_size_type, vertices_size_type> marker_type;
public:
    LB_improved_contraction_base(G_t &G)
      : algo1(CFG_t::name()), _g(G), _lb(0),
        _marker(boost::num_vertices(G))
    {
    }

    void do_it(){
        timer_on();

        G_t H(_g);
        int lb = CFG_t::lb_algo(H);

        while(true){
            G_t H;
            boost::copy_graph(_g, H);
            CFG_t::improvement_algo(H, lb+1);

            int new_lb=0;

            while(boost::num_edges(H) > 0){
                new_lb = CFG_t::lb_algo(H);
                if(new_lb > lb){
                    break;
                }else{
                }

                auto min_vertex=get_min_degree_vertex(H, true); //ignore isolated vertices
                auto w=get_least_common_vertex(min_vertex, _marker, H);

                contract_edge(min_vertex, w, H);

                CFG_t::improvement_algo(H, lb+1);
            }
            if(new_lb > lb){
                lb++;
            }
            else{
                break;
            }
        }

        _lb = lb;

        timer_off();
    }


    unsigned lower_bound_bagsize(){
        return _lb+1u;
    }

private:
    G_t& _g;
    unsigned _lb;
    marker_type _marker;
};

} //namespace impl

template <typename G_t>
struct CFG_LBNC_deltaD{
    static int lb_algo(G_t &H){ untested();
        impl::deltaD<G_t> deltaD(H);
        deltaD.do_it();
        return (int)deltaD.lower_bound_bagsize() - 1;
    }

    static void improvement_algo(G_t &H, unsigned k){ untested();
        treedec::lb::k_neighbour_improved_graph(H, k);
    }

    static const std::string name(){ untested();
        return "lb::LBNC_deltaD";
    }
};

template <typename G_t>
int LBNC_deltaD(G_t &G)
{ untested();
    unsigned int V = boost::num_vertices(G);
    unsigned int E = boost::num_edges(G);

    if(V == 0){ untested();
        return -1;
    }else if(E == 0){ untested();
        return 0;
    }else if(2*E == V*(V-1u)){ untested();
        return V-1u;
    }else{ untested();
    }

    impl::LB_improved_contraction_base<G_t, CFG_LBNC_deltaD<G_t> > LBNC_deltaD(G);
    LBNC_deltaD.do_it();
    return (int)LBNC_deltaD.lower_bound_bagsize()-1;
}

template <typename G_t>
int LBNC_deltaD(G_t const&G)
{ untested();
    G_t H(G);
    return LBNC_deltaD(H);
}

template <typename G_t>
struct CFG_LBNC_deltaC{
    static int lb_algo(G_t &H){
        impl::deltaC_least_c<G_t> deltaC(H);
        deltaC.do_it();
        return (int)deltaC.lower_bound_bagsize() - 1;
    }

    static void improvement_algo(G_t &H, unsigned k){
        treedec::lb::k_neighbour_improved_graph(H, k);
    }

    static const std::string name(){
        return "lb::LBNC_deltaC";
    }
};

template <typename G_t>
int LBNC_deltaC(G_t &G){
    unsigned int V = boost::num_vertices(G);
    unsigned int E = boost::num_edges(G);

    if(V == 0){
        return -1;
    }
    else if(E == 0){
        return 0;
    }
    else if(2*E == V*(V-1u)){
        return V-1u;
    }

    impl::LB_improved_contraction_base<G_t, CFG_LBNC_deltaC<G_t> > LBNC_deltaC(G);
    LBNC_deltaC.do_it();
    return (int)LBNC_deltaC.lower_bound_bagsize()-1;
}


template <typename G_t>
void k_path_improved_graph(G_t &G, unsigned int k){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> edges_to_add;

    typename boost::graph_traits<G_t>::vertex_iterator vIt1, vIt2, vEnd;
    for(boost::tie(vIt1, vEnd) = boost::vertices(G); vIt1 != vEnd; vIt1++){
        vIt2 = vIt1;
        vIt2++;
        for(; vIt2 != vEnd; vIt2++){
            if(!boost::edge(*vIt1, *vIt2, G).second){
                typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> X, Y, S;

                typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt1, G); nIt != nEnd; nIt++){
                    X.insert(*nIt);
                }
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt2, G); nIt != nEnd; nIt++){
                    Y.insert(*nIt);
                }

                std::vector<BOOL> disabled(boost::num_vertices(G), false);
                auto pos1 = boost::get(boost::vertex_index, G, *vIt1);
                auto pos2 = boost::get(boost::vertex_index, G, *vIt2);

                unsigned num_dis=0;
                if(!disabled[pos1]) ++num_dis;
                if(!disabled[pos2]) ++num_dis;

                disabled[pos1] = true;
                disabled[pos2] = true;

                treedec::seperate_vertices(G, disabled, num_dis, X, Y, S);

                if(S.size() >= k){
                    edges_to_add.push_back(*vIt1);
                    edges_to_add.push_back(*vIt2);
                }
            }
        }
    }

    for(unsigned int i = 0; i < edges_to_add.size(); ){
        boost::add_edge(edges_to_add[i], edges_to_add[i+1], G);
        ++i;
        ++i;
    }
}

template <typename G_t>
struct CFG_LBP_deltaD{
    static int lb_algo(G_t &H){ untested();
        impl::deltaD<G_t> deltaD(H);
        deltaD.do_it();
        return (int)deltaD.lower_bound_bagsize() - 1;
    }

    static void improvement_algo(G_t &H, unsigned k){ untested();
        treedec::lb::k_path_improved_graph(H, k);
    }

    static const std::string name(){ untested();
        return "lb::LBP_deltaD";
    }

};

template <typename G_t>
int LBP_deltaD(G_t &G)
{ untested();
    unsigned int V = boost::num_vertices(G);
    unsigned int E = boost::num_edges(G);

    if(V == 0){ untested();
        return -1;
    }else if(E == 0){ untested();
        return 0;
    }else if(2*E == V*(V-1u)){ untested();
        return V-1u;
    }else{ untested();
    }

    impl::LB_improved_base<G_t, CFG_LBP_deltaD<G_t> > LBP_deltaD(G);
    LBP_deltaD.do_it();
    return (int)LBP_deltaD.lower_bound_bagsize()-1;
}

template <typename G_t>
int LBP_deltaD(const G_t &G)
{ untested();
    G_t H(G);
    return LBNC_deltaD(G);
}


template <typename G_t>
struct CFG_LBP_deltaC{
    static int lb_algo(G_t &H){
        impl::deltaC_least_c<G_t> deltaC(H);
        deltaC.do_it();
        return (int)deltaC.lower_bound_bagsize() - 1;
    }

    static void improvement_algo(G_t &H, unsigned k){
        treedec::lb::k_path_improved_graph(H, k);
    }

    static const std::string name(){
        return "lb::LBP_deltaC";
    }
};

template <typename G_t>
int LBP_deltaC(G_t &G){
    unsigned int V = boost::num_vertices(G);
    unsigned int E = boost::num_edges(G);

    if(V == 0){
        return -1;
    }
    else if(E == 0){
        return 0;
    }
    else if(2*E == V*(V-1u)){
        return V-1u;
    }

    impl::LB_improved_base<G_t, CFG_LBP_deltaC<G_t> > LBP_deltaC(G);
    LBP_deltaC.do_it();
    return (int)LBP_deltaC.lower_bound_bagsize()-1;
}

template <typename G_t>
struct CFG_LBPC_deltaD{
    static int lb_algo(G_t &H){ untested();
        impl::deltaD<G_t> deltaD(H);
        deltaD.do_it();
        return (int)deltaD.lower_bound_bagsize() - 1;
    }

    static void improvement_algo(G_t &H, unsigned k){ untested();
        treedec::lb::k_path_improved_graph(H, k);
    }

    static const std::string name(){ untested();
        return "lb::LBPC_deltaD";
    }
};

template <typename G_t>
int LBPC_deltaD(G_t &G)
{ untested();
    unsigned int V = boost::num_vertices(G);
    unsigned int E = boost::num_edges(G);

    if(V == 0){ untested();
        return -1;
    }else if(E == 0){ untested();
        return 0;
    }else if(2*E == V*(V-1u)){ untested();
        return V-1u;
    }

    impl::LB_improved_contraction_base<G_t, CFG_LBPC_deltaD<G_t> > LBPC_deltaD(G);
    LBPC_deltaD.do_it();
    return (int)LBPC_deltaD.lower_bound_bagsize()-1;
}

template <typename G_t>
int LBPC_deltaD(const G_t &G)
{ untested();
    G_t H(G); // does this work with all graphs?
    return LBPC_deltaD(H);
}

template <typename G_t>
struct CFG_LBPC_deltaC{
    static int lb_algo(G_t &H){
        impl::deltaC_least_c<G_t> deltaC(H);
        deltaC.do_it();
        return (int)deltaC.lower_bound_bagsize() - 1;
    }

    static void improvement_algo(G_t &H, unsigned k){
        treedec::lb::k_path_improved_graph(H, k);
    }

    static const std::string name(){
        return "lb::LBPC_deltaC";
    }
};

template <typename G_t>
int LBPC_deltaC(G_t &G){
    unsigned int V = boost::num_vertices(G);
    unsigned int E = boost::num_edges(G);

    if(V == 0){
        return -1;
    }
    else if(E == 0){
        return 0;
    }
    else if(2*E == V*(V-1u)){
        return V-1u;
    }else{
    }

    impl::LB_improved_contraction_base<G_t, CFG_LBPC_deltaC<G_t> > LBPC_deltaC(G);
    LBPC_deltaC.do_it();
    return (int)LBPC_deltaC.lower_bound_bagsize()-1;
}


/* Maximum Cardinality Search-based */

namespace detail{

//Applies a maximum cardinality search on G and returns the vertex descriptors of the maximum visited
//degree-vertex and the last visited vertex.
template <typename G_t>
typename std::pair<int, typename boost::graph_traits<G_t>::vertex_descriptor> MCS(G_t &G){ untested();
    std::vector<int> visited_degree(boost::num_vertices(G), 0);

    int max = -1;
    typename boost::graph_traits<G_t>::vertex_descriptor max_vertex=*(boost::vertices(G).first);

    for(unsigned int i = 0; i < boost::num_vertices(G); i++){ untested();
        int cur_max = -1;
        unsigned int cur_pos = 0;
        for(unsigned int j = 0; j < visited_degree.size(); j++){ untested();
            if(visited_degree[j] > cur_max){ untested();
                cur_max = visited_degree[j];
                cur_pos = j;
            }else{ untested();
            }
        }

        typename boost::graph_traits<G_t>::vertex_iterator cur_vertex_it = boost::vertices(G).first;
        std::advance(cur_vertex_it, cur_pos);

        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*cur_vertex_it, G); nIt != nEnd; nIt++){ untested();
            auto pos = boost::get(boost::vertex_index, G, *nIt);
            if(visited_degree[pos] > -1){ untested();
                visited_degree[pos]++;
            }else{ untested();
            }
        }

        if(visited_degree[cur_pos] > max){ untested();
            max = visited_degree[cur_pos];
            max_vertex = *cur_vertex_it;
        }else{ untested();
        }

        visited_degree[cur_pos] = -1;
    }

    return std::make_pair(max, max_vertex);
}
} //namespace detail (for MCS)



//Returns the maximum visited degree in a maximum cardinality search.
template <typename G_t>
int MCS(G_t &G){ untested();
    return treedec::lb::detail::MCS(G).first;
}


//This is the contraction version of the maximal cardinality search lower bound algorithm.
//Heuristic: max_mcs.
template <typename G_t>
int MCSC(G_t& G)
{ untested();
    int max = -1;

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    while(boost::num_edges(G) > 0){ untested();
        typename std::pair<int, typename boost::graph_traits<G_t>::vertex_descriptor> result = treedec::lb::detail::MCS(G);
        max = (result.first > max)? result.first : max;
        typename boost::graph_traits<G_t>::vertex_descriptor v = result.second;

        auto w=get_least_common_vertex(v, G);

        contract_edge(w, v, G);
    }

    return max;
}


// why int?
// why not average_degree?
template <typename G_t>
int relation_edges_vertices(G_t &G){ untested();
    if(boost::num_vertices(G) == 0){ untested();
        return -1;
    }else{ untested();
    }
    return (int)(2*boost::num_edges(G)/boost::num_vertices(G));
}

} //namespace lb

} //namespace treedec

#endif //TREEDEC_LOWER_BOUNDS_HPP

// vim:ts=8:sw=4:et
