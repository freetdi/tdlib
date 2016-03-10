// Lukas Larisch, 2014 - 2016
// Felix Salfelder, 2016
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

#ifndef TD_LOWER_BOUNDS
#define TD_LOWER_BOUNDS

#include <set>
#include <vector>
#include <climits>

#include <boost/graph/adjacency_list.hpp>
#include <utility>

#include "TD_simple_graph_algos.hpp"
#include "TD_NetworkFlow.hpp"
#include "TD_misc.hpp"
#include "TD_noboost.hpp"

namespace treedec{

namespace lb{


/* DEGREE BASED */

//Smallest vertex-degree in G.
template <typename G_t>
int delta(const G_t &G){
    if(boost::num_vertices(G) == 0){
        return -1;
    }

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    unsigned int min = boost::num_vertices(G);
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::degree(*vIt, G);
        min = (degree < min)? degree : min;
    }
    return (int)min;
}

//Second smallest vertex-degree in G.
template <typename G_t>
int delta2(const G_t &G){
    if(boost::num_vertices(G) == 0){
        return -1;
    }
    else if(boost::num_vertices(G) == 1){
        return 0;
    }

    unsigned int min = boost::num_vertices(G);
    unsigned int snd = min;

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::degree(*vIt, G);
        if(degree <= min){
            snd = min;
            min = degree;
        }
        if(degree > min && degree < snd){
            snd = degree;
        }
    }
    return (int)snd;
}

template <typename G_t>
int gamma(const G_t &G)
{
    BOOST_AUTO(V, boost::num_vertices(G));
    BOOST_AUTO(E, boost::num_edges(G));

    if(V == 0){
        return -1;
    }else if(E == 0){
        return 0;
    }else if(2*E+1 == V*(V-1u)) {
        return V-1;
    }

    //Sort the vertices of G according to rising degree -> degree_sequence.
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> degree_sequence;
    noboost::make_degree_sequence(G, degree_sequence);

    //Take the degree of the right vertex in the first not-edge.
    for(unsigned int i = 0; i < boost::num_vertices(G); i++){
        for(unsigned int j = 0; j < i; j++){
            if(!boost::edge(degree_sequence[i], degree_sequence[j], G).second){
                unsigned int deg = boost::degree(degree_sequence[i], G);
                return (int)deg;
            }
        }
    }
    unreachable();
    return 0;
}


namespace detail{

template <typename G_t>
int deltaD(G_t &G){
    unsigned int maxmin = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::vertex_descriptor min_vertex = *(boost::vertices(G).second);

    while(true){
        unsigned int min_degree = boost::num_vertices(G);
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int degree = boost::degree(*vIt, G);
            if(degree < min_degree && degree > 0){
                min_degree = degree;
                min_vertex = *vIt;
            }
        }

        if(min_degree == boost::num_vertices(G)){
            return (int)maxmin;
        }

        maxmin = (maxmin>min_degree)? maxmin : min_degree;

        boost::clear_vertex(min_vertex, G);
    }
}

} //namespace detail (for deltaD)

template <typename G_t>
int deltaD(G_t& G)
{
    if(boost::num_vertices(G) == 0){
        return -1;
    }

    return treedec::lb::detail::deltaD(G);
}

//Assume each vertex as the minimal one and apply deltaD.
template <typename G_t>
int delta2D(const G_t &G){
    if(boost::num_vertices(G) == 0){
        return -1;
    }

    G_t H;
    boost::copy_graph(G, H);

    typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> assumed_minimal;

    typename boost::graph_traits<G_t>::vertex_iterator hIt, hEnd;
    for(boost::tie(hIt, hEnd) = boost::vertices(H); hIt != hEnd; hIt++){
        assumed_minimal.push_back(*hIt);
    }

    unsigned int min_degree, maxmin = 0;
    typename boost::graph_traits<G_t>::vertex_descriptor min_vertex = *(boost::vertices(H).first);

    for(unsigned int i = 0; i < assumed_minimal.size(); i++){
        while(boost::num_edges(H) > 0){
            min_degree = boost::num_vertices(H);

            for(boost::tie(hIt, hEnd) = boost::vertices(H); hIt != hEnd; hIt++){
                if(*hIt == assumed_minimal[i]){
                    continue;
                }
                unsigned int degree = boost::degree(*hIt, H);
                if(degree < min_degree && degree > 0){
                    min_degree = degree;
                    min_vertex = *hIt;
                }
            }
            if(min_degree == boost::num_vertices(H)){
                break;
            }

            maxmin = (maxmin>min_degree)? maxmin : min_degree;
            boost::clear_vertex(min_vertex, H);
        }
        H.clear();
        boost::copy_graph(G, H);
    }
    return (int)maxmin;
}

namespace detail{

template <typename G_t>
void gammaD_left(G_t &G, unsigned int &lb){
    if(boost::num_edges(G) == 0){
        return;
    }

    //Sort the vertices of G according to rising degree -> degree_sequence.
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> degree_sequence;
    noboost::make_degree_sequence(G, degree_sequence);

    for(unsigned int i = 0; i < degree_sequence.size(); i++){
        for(unsigned int j = 0; j < i; j++){
            if(boost::edge(degree_sequence[i], degree_sequence[j], G).second){
                continue;
            }

            //gammaD-left heuristic
            unsigned int degree = boost::degree(degree_sequence[i], G);
            for(unsigned int k = 0; k < i; k++){
                boost::clear_vertex(degree_sequence[k], G);
            }

            lb = (degree > lb)? degree : lb;

            gammaD_left(G, lb);
            return;
        }
    }
}

} //namespace detail (for gammaD_left)

template <typename G_t>
int gammaD_left(G_t& G)
{
    BOOST_AUTO(V, boost::num_vertices(G));
    BOOST_AUTO(E, boost::num_edges(G));

    if(V == 0){
        return -1;
    }else if(E == 0){
        return 0;
    }else if(2*E == V*(V-1u)) {
        return V-1;
    }else{
        unsigned int lb = 0;
        detail::gammaD_left(G, lb);
        return lb;
    }
}

namespace detail{

template <typename G_t>
void gammaD_right(G_t &G, unsigned int &lb){
    if(boost::num_edges(G) == 0){
        return;
    }

    //Sort the vertices of G according to rising degree -> degree_sequence.
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> degree_sequence;
    noboost::make_degree_sequence(G, degree_sequence);

    for(unsigned int i = 0; i < degree_sequence.size(); i++){
        for(unsigned int j = 0; j < i; j++){
            if(boost::edge(degree_sequence[i], degree_sequence[j], G).second){
                continue;
            }

            //gammaD-right heuristic
            unsigned int degree = boost::degree(degree_sequence[i], G);
            boost::clear_vertex(degree_sequence[i], G);

            lb = (degree > lb)? degree : lb;

            gammaD_right(G, lb);
            return;
        }
    }
}

} //namespace detail (for gammaD_right)


template <typename G_t>
int gammaD_right(G_t& G)
{
    BOOST_AUTO(V, boost::num_vertices(G));
    BOOST_AUTO(E, boost::num_edges(G));

    if(V == 0){
        return -1;
    }else if(E == 0){
        return 0;
    }else if(2*E == V*(V-1u)) {
        return V-1;
    }else{
        unsigned int lb = 0;
        detail::gammaD_right(G, lb);
        return lb;
    }
}

namespace detail{

template <typename G_t>
void gammaD_min_e(G_t &G, unsigned int &lb){
    if(boost::num_edges(G) == 0){
        return;
    }

    //Sort the vertices of G according to rising degree -> degree_sequence.
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> degree_sequence;
    noboost::make_degree_sequence(G, degree_sequence);

    for(unsigned int i = 0; i < degree_sequence.size(); i++){
        for(unsigned int j = 0; j < i; j++){
            if(boost::edge(degree_sequence[i], degree_sequence[j], G).second){
                continue;
            }

            //gammaD-min-e heuristic
            unsigned int degree_right = boost::degree(degree_sequence[i], G);
            unsigned int degree_left = 0;
            for(unsigned int k = 0; k < i; k++){
                degree_left += boost::degree(degree_sequence[k], G);
            }

            if(degree_left < degree_right){
                for(unsigned int t = 0; t < i; t++){
                    boost::clear_vertex(degree_sequence[t], G);
                }
            }
            else{
                boost::clear_vertex(degree_sequence[i], G);
            }

            lb = (degree_right > lb)? degree_right : lb;

            gammaD_min_e(G, lb);
            return;
        }
    }
}

} //namespace detail (for gammaD_min_e)

template <typename G_t>
int gammaD_min_e(G_t& G)
{
    BOOST_AUTO(V, boost::num_vertices(G));
    BOOST_AUTO(E, boost::num_edges(G));

    if(V == 0){
        return -1;
    }else if(E == 0){
        return 0;
    }else if(2*E+1 == V*(V-1u)) {
        return V-1;
    }else{
        unsigned int lb = 0;
        detail::gammaD_min_e(G, lb);
        return lb;
    }
}

namespace detail{

template <typename G_t>
int deltaC_min_d(G_t &G){
    unsigned int lb = 0;

    while(boost::num_edges(G) > 0){
        //Search a minimum-degree-vertex.
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex
                      = noboost::get_min_degree_vertex(G, true); //ignore isolated vertices

        lb = (lb>boost::degree(min_vertex, G))? lb : boost::degree(min_vertex, G);

        //min_d heuristic: Search a neighbour of min_vertex with minimal degree.
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        unsigned int min_degree_w = boost::num_vertices(G);

        boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, G);
        typename boost::graph_traits<G_t>::vertex_descriptor w = *nIt;
        for(; nIt != nEnd; nIt++){
            unsigned int degree = boost::degree(*nIt, G);
            if(degree <= min_degree_w){
                min_degree_w = degree;
                w = *nIt;
            }
        }

        noboost::contract_edge(min_vertex, w, G, false);
    }
    return (int)lb;
}

} //namespace detail (for deltaC_min_d)


template <typename G_t>
int deltaC_min_d(G_t& G)
{
    BOOST_AUTO(V, boost::num_vertices(G));
    BOOST_AUTO(E, boost::num_edges(G));

    if(V == 0){
        return -1;
    }else if(E == 0){
        return 0;
    }else if(2*E+1 == V*(V-1u)) {
        return V-1;
    }else{
        return detail::deltaC_min_d(G);
    }
}

namespace detail{

template <typename G_t>
int deltaC_max_d(G_t &G){
    unsigned int lb = 0;

    while(boost::num_edges(G) > 0){
        //Search a minimum-degree-vertex.
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex
                      = noboost::get_min_degree_vertex(G, true); //ignore isolated vertices

        lb = (lb>boost::degree(min_vertex, G))? lb : boost::degree(min_vertex, G);

        //max_d heuristic: Search the neighbour of min_vertex with maximal degree.
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        unsigned int max_degree = 0;

        boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, G);
        typename boost::graph_traits<G_t>::vertex_descriptor w = *nIt;
        for(; nIt != nEnd; nIt++){
            unsigned int degree = boost::degree(*nIt, G);
            if(degree > max_degree){
                max_degree = degree;
                w = *nIt;
            }
        }

        noboost::contract_edge(min_vertex, w, G, false);
    }
    return (int)lb;
}

} //namespace detail (for deltaC_max_d)

template <typename G_t>
int deltaC_max_d(G_t& G)
{
    BOOST_AUTO(V, boost::num_vertices(G));
    BOOST_AUTO(E, boost::num_edges(G));

    if(V == 0){
        return -1;
    }else if(E == 0){
        return 0;
    }else if(2*E == V*(V-1u)) {
        return V-1u;
    }else{
        return detail::deltaC_max_d(G);
    }
}

template<typename G_t>
struct degree_decrease : public noboost::vertex_callback<G_t>{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    degree_decrease(std::vector<std::set<vertex_descriptor> >*d, G_t*g) :
        degs(d), G(g){}

    void operator()(vertex_descriptor v){
        size_t degree = boost::degree(v, *G);
        if(degree==0){
            // unreachable
            // unconnected nodes are unreachable throught adj iterator
            assert(false);
        }else if(degree==1){
            // unreachable
            // a degree one node does not change its degree during collapse
            assert(false);
        }else{
            size_t found=(*degs)[degree].erase(v);
            assert(found); // sanity check on degs.
            (void) found;
            bool done=(*degs)[degree-1].insert(v).second;
            assert(done);
            (void) done;
        }
    }
    private:
        std::vector<std::set<vertex_descriptor> >*degs;
        G_t* G;
};

namespace detail{

template <typename G_t>
int deltaC_least_c(G_t &G)
{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;

    unsigned int lb = 0;
    misc::DEGS<G_t> degs(G);
    degree_decrease<G_t> cb(&degs._degs, &G);

    unsigned int min_degree = 1;

    while(boost::num_edges(G) > 0){
        //Search a minimum-degree-vertex.
        if(degs[min_degree].empty()){
            for(min_degree = 1; min_degree < degs.size(); min_degree++){
                if(!degs[min_degree].empty()){
                    break;
                }
            }
        }

        if(lb <= min_degree){
            lb = min_degree;
        }

        vertex_descriptor min_vertex;
        min_vertex = *degs[min_degree].begin();

        //least-c heuristic: search the neighbour of min_vertex such that
        //contracting {min_vertex, w} removes the least edges
        typename boost::graph_traits<G_t>::vertex_descriptor w = noboost::get_least_common_vertex(min_vertex, G);

        size_t outdegw = boost::degree(w, G);
        size_t outdegmin = boost::degree(min_vertex, G);
        assert(degs[outdegw].find(w) != degs[outdegw].end());
        assert(degs[outdegmin].find(min_vertex) != degs[outdegmin].end());

        degs[outdegw].erase(w);
        degs[outdegmin].erase(min_vertex);

        //Contract the edge between min_vertex into w.
        //Clear min_vertex and rearrange degs through callback.
        noboost::contract_edge(min_vertex, w, G, false, &cb);

        assert(0==boost::degree(min_vertex, G));
        assert(boost::degree(min_vertex, G)==0);

        outdegw = boost::degree(w, G);
        degs[outdegw].insert(w);
    }
    return (int)lb;
}

} //namespace detail (for deltaC_least_c)


template <typename G_t>
int deltaC_least_c(G_t& G)
{
    if(boost::num_vertices(G) == 0){
        return -1;
    }
    else if(boost::num_edges(G) == 0){
        return 0;
    }
    else if(2*boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1u)){
        return boost::num_vertices(G)-1;
    }

    return treedec::lb::detail::deltaC_least_c(G);
}


/* IMPROVED GRAPHS */

template <typename G_t>
void k_neighbour_improved_graph(G_t &G, unsigned int k){
    std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > edges_to_add;

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

                std::set_intersection(N1.begin(), N1.end(), N2.begin(), N2.end(), std::inserter(intersection, intersection.begin()));

                if(intersection.size() >= k){
                    std::set<typename boost::graph_traits<G_t>::vertex_descriptor> edge;
                    edge.insert(*vIt1); edge.insert(*vIt2);
                    edges_to_add.push_back(edge);
                }
            }
        }
    }

    for(unsigned int i = 0; i < edges_to_add.size(); i++){
        boost::add_edge(*(edges_to_add[i].begin()), *(++edges_to_add[i].begin()), G);
    }
}


template <typename G_t>
int LBN_deltaD(const G_t &G){
    if(boost::num_vertices(G) == 0){
        return -1;
    }
    else if(boost::num_edges(G) == 0){
        return 0;
    }
    else if(2*boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return boost::num_vertices(G)-1;
    }

    int lb = treedec::lb::deltaD(G);

    while(true){
        G_t H;
        boost::copy_graph(G, H);
        treedec::lb::k_neighbour_improved_graph(H, lb+1);

        if(treedec::lb::detail::deltaD(H) > lb){
            lb++;
        }
        else{
            break;
        }
    }
    return lb;
}

template <typename G_t>
int LBN_deltaC(G_t &G){
    if(boost::num_vertices(G) == 0){
        return -1;
    }
    else if(boost::num_edges(G) == 0){
        return 0;
    }
    else if(2*boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return boost::num_vertices(G)-1;
    }

    int lb = treedec::lb::deltaC_least_c(G);

    while(true){
        G_t H;
        boost::copy_graph(G, H);
        treedec::lb::k_neighbour_improved_graph(H, lb+1);

        if(treedec::lb::detail::deltaC_least_c(H) > lb){
            lb++;
        }
        else{
            break;
        }
    }
    return lb;
}

template <typename G_t>
int LBNC_deltaD(const G_t &G){
    if(boost::num_vertices(G) == 0){
        return -1;
    }
    else if(boost::num_edges(G) == 0){
        return 0;
    }
    else if(2*boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return boost::num_vertices(G)-1;
    }

    int lb = deltaD(G);

    while(true){
        G_t H;
        boost::copy_graph(G, H);
        treedec::lb::k_neighbour_improved_graph(H, lb+1);

        int new_lb;

        while(boost::num_edges(H) > 0){
            new_lb = treedec::lb::deltaD(H);
            if(new_lb > lb){
                break;
            }

            typename boost::graph_traits<G_t>::vertex_descriptor min_vertex
                      = noboost::get_min_degree_vertex(H, true); //ignore isolated vertices

            typename boost::graph_traits<G_t>::vertex_descriptor w = noboost::get_least_common_vertex(min_vertex, H);

            noboost::contract_edge(min_vertex, w, H, false);

            treedec::lb::k_neighbour_improved_graph(H, lb+1);
        }
        if(new_lb > lb){
            lb++;
        }else{
            break;
        }
    }
    return lb;
}

template <typename G_t>
int LBNC_deltaC(G_t &G){
    if(boost::num_vertices(G) == 0){
        return -1;
    }
    else if(boost::num_edges(G) == 0){
        return 0;
    }
    else if(2*boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return boost::num_vertices(G)-1;
    }

    int lb = treedec::lb::deltaC_least_c(G);

    while(true){
        G_t H;
        boost::copy_graph(G, H);
        treedec::lb::k_neighbour_improved_graph(H, lb+1);

        int new_lb = 0;

        while(boost::num_edges(H) > 0){
            new_lb = deltaC_least_c(H);
            if(new_lb > lb){
                break;
            }

            typename boost::graph_traits<G_t>::vertex_descriptor min_vertex
                      = noboost::get_min_degree_vertex(H, true); //ignore isolated vertices

            typename boost::graph_traits<G_t>::vertex_descriptor w = noboost::get_least_common_vertex(min_vertex, H);

            noboost::contract_edge(min_vertex, w, H, false);

            treedec::lb::k_neighbour_improved_graph(H, lb+1);
        }
        if(new_lb > lb){
            lb++;
        }
        else{
            break;
        }
    }
    return lb;
}


template <typename G_t>
void k_path_improved_graph(G_t &G, unsigned int k){
    std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > edges_to_add;

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

                std::vector<bool> disabled(boost::num_vertices(G), false);
                unsigned int pos1 = noboost::get_pos(*vIt1, G);
                unsigned int pos2 = noboost::get_pos(*vIt2, G);
                disabled[pos1] = true;
                disabled[pos2] = true;

                treedec::seperate_vertices(G, disabled, X, Y, S);

                if(S.size() >= k){
                    std::set<typename boost::graph_traits<G_t>::vertex_descriptor> edge;
                    edge.insert(*vIt1); edge.insert(*vIt2);
                    edges_to_add.push_back(edge);
                }
            }
        }
    }

    for(unsigned int i = 0; i < edges_to_add.size(); i++){
        boost::add_edge(*(edges_to_add[i].begin()), *(++edges_to_add[i].begin()), G);
    }
}

template <typename G_t>
int LBP_deltaD(const G_t &G){
    if(boost::num_vertices(G) == 0){
        return -1;
    }
    else if(boost::num_edges(G) == 0){
        return 0;
    }
    else if(2*boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return boost::num_vertices(G)-1;
    }

    int lb = treedec::lb::deltaD(G);

    while(true){
        G_t H;
        boost::copy_graph(G, H);

        treedec::lb::k_path_improved_graph(H, lb+1);

        if(treedec::lb::detail::deltaD(H) > lb){
            lb++;
        }
        else{
            break;
        }
    }
    return lb;
}

template <typename G_t>
int LBP_deltaC(G_t &G){
    if(boost::num_vertices(G) == 0){
        return -1;
    }
    else if(boost::num_edges(G) == 0){
        return 0;
    }
    else if(2*boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return boost::num_vertices(G)-1;
    }

    int lb = treedec::lb::deltaC_least_c(G);

    while(true){
        G_t H;
        boost::copy_graph(G, H);

        treedec::lb::k_path_improved_graph(H, lb+1);

        if(treedec::lb::detail::deltaC_least_c(H) > lb){
            lb++;
        }
        else{
            break;
        }
    }
    return lb;
}

template <typename G_t>
int LBPC_deltaD(const G_t &G){
    if(boost::num_vertices(G) == 0){
        return -1;
    }
    else if(boost::num_edges(G) == 0){
        return 0;
    }
    else if(2*boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return boost::num_vertices(G)-1;
    }

    int lb = treedec::lb::deltaD(G);

    while(true){
        G_t H;
        boost::copy_graph(G, H);

        treedec::lb::k_path_improved_graph(H, lb+1);

        int new_lb;

        while(boost::num_edges(H) > 0){
            new_lb = treedec::lb::deltaD(H);
            if(new_lb > lb){
                break;
            }

            typename boost::graph_traits<G_t>::vertex_descriptor min_vertex
                      = noboost::get_min_degree_vertex(H, true); //ignore isolated vertices

            typename boost::graph_traits<G_t>::vertex_descriptor w = noboost::get_least_common_vertex(min_vertex, H);

            noboost::contract_edge(min_vertex, w, H, false);

            treedec::lb::k_path_improved_graph(H, lb+1);

        }
        if(new_lb > lb){
            lb++;
        }
        else{
            break;
        }
    }
    return lb;
}

template <typename G_t>
int LBPC_deltaC(G_t &G){
    if(boost::num_vertices(G) == 0){
        return -1;
    }
    else if(boost::num_edges(G) == 0){
        return 0;
    }
    else if(2*boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return boost::num_vertices(G)-1;
    }

    int lb = treedec::lb::deltaC_least_c(G);

    while(true){
        G_t H;
        boost::copy_graph(G, H);

        treedec::lb::k_path_improved_graph(H, lb+1);

        int new_lb = 0;

        while(boost::num_edges(H) > 0){
            new_lb = treedec::lb::deltaC_least_c(H);

            if(new_lb > lb){
                break;
            }

            typename boost::graph_traits<G_t>::vertex_descriptor min_vertex
                      = noboost::get_min_degree_vertex(H, true); //ignore isolated vertices

            typename boost::graph_traits<G_t>::vertex_descriptor w = noboost::get_least_common_vertex(min_vertex, H);

            noboost::contract_edge(min_vertex, w, H, false);

            treedec::lb::k_path_improved_graph(H, lb+1);
        }
        if(new_lb > lb){
            lb++;
        }
        else{
            break;
        }
    }
    return lb;
}


/* Maximum Cardinality Search-based */

namespace detail{

//Applies a maximum cardinality search on G and returns the vertex descriptors of the maximum visited
//degree-vertex and the last visited vertex.
template <typename G_t>
typename std::pair<int, typename boost::graph_traits<G_t>::vertex_descriptor> MCS(G_t &G){
    std::vector<int> visited_degree(boost::num_vertices(G), 0);

    int max = -1;
    BOOST_AUTO(vertices, boost::vertices(G));
    typename boost::graph_traits<G_t>::vertex_descriptor max_vertex=*vertices.first;

    for(unsigned int i = 0; i < boost::num_vertices(G); i++){
        int cur_max = -1;
        unsigned int cur_pos = 0;
        for(unsigned int j = 0; j < visited_degree.size(); j++){
            if(visited_degree[j] > cur_max){
                cur_max = visited_degree[j];
                cur_pos = j;
            }
        }

        typename boost::graph_traits<G_t>::vertex_iterator cur_vertex_it = vertices.first;
        std::advance(cur_vertex_it, cur_pos);

        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*cur_vertex_it, G); nIt != nEnd; nIt++){
            unsigned int pos = noboost::get_pos(*nIt, G);
            if(visited_degree[pos] > -1){
                visited_degree[pos]++;
            }
        }

        if(visited_degree[cur_pos] > max){
            max = visited_degree[cur_pos];
            max_vertex = *cur_vertex_it;
        }

        visited_degree[cur_pos] = -1;
    }

    return std::make_pair(max, max_vertex);
}
} //namespace detail (for MCS)



//Returns the maximum visited degree in a maximum cardinality search.
template <typename G_t>
int MCS(G_t &G){
    return treedec::lb::detail::MCS(G).first;
}


//This is the contraction version of the maximal cardinality search lower bound algorithm.
//Heuristic: max_mcs.
template <typename G_t>
int MCSC(G_t& G)
{
    int max = -1;

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    while(boost::num_edges(G) > 0){
        typename std::pair<int, typename boost::graph_traits<G_t>::vertex_descriptor> result = treedec::lb::detail::MCS(G);
        max = (result.first > max)? result.first : max;
        typename boost::graph_traits<G_t>::vertex_descriptor v = result.second;

        typename boost::graph_traits<G_t>::vertex_descriptor w = noboost::get_least_common_vertex(v, G);

        noboost::contract_edge(w, v, G);
    }

    return max;
}


template <typename G_t>
int relation_edges_vertices(G_t &G){
    if(boost::num_vertices(G) == 0){
        return -1;
    }
    return (int)(2*boost::num_edges(G)/boost::num_vertices(G));
}

} //namespace lb

} //namespace treedec

#endif //TD_LOWER_BOUNDS

// vim:ts=8:sw=4:et
