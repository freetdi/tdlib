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
// Offers functionality to compute tree decompositions of graphs
// according to various heuristic and functions, that convert
// tree decompositions to elimination orderings and vice versa.
// Also the LEX-M algorithm is included in this header
//
//

/*
 These functions are most likely to be interesting for outside use:

 - void minDegree_decomp(G_t G&, T_t &T)
 - void fillIn_decomp(G_t G&, T_t &T)
 - void minDegree_ordering(G_t G,
         typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering)
 - void fillIn_ordering(G_t G,
         typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering)
 - void ordering_to_treedec(G_t G,
         typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering, T_t &T)
 - void treedec_to_ordering<G_t, T_t>(T_t T,
         typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering)
 - void LEX_M_minimal_ordering(G_t &G,
         typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering)
*/

#ifndef TD_ELIMINATION_ORDERING
#define TD_ELIMINATION_ORDERING

#include <cmath>
#include <climits>

#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "TD_preprocessing.hpp"
#include "TD_simple_graph_algos.hpp"
#include "TD_misc.hpp"
#include "TD_noboost.hpp"
#include "TD_std.hpp"
//#include "TD_trace.hpp"

namespace treedec{

namespace detail{

//DRAFT (no useful interface).
template<typename G>
struct degree_mod : public noboost::vertex_callback<G>{
    typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
    degree_mod(misc::DEGS<G>* d, G* g) : _degs(d), _g(g){}

    //Reinsert with degree-1.
    void operator()(vertex_descriptor v){
        unsigned deg=boost::degree(v,*_g);
        (void)deg;
        assert(deg);
        bool done=(*_degs)[deg-1].insert(v).second;
        (void)done;
        assert(done);
    }
    //private: not yet.
//    BUG:: hardcoded type
        misc::DEGS<G>* _degs;
    private:
        degree_mod(const degree_mod&){}
        G* _g;
};

} //namespace detail

// register a 1-neigborhood to DEGS
template<class U, class G, class B, class D>
void redegree(U, G& g, B& neighborhood, D& d)
{
    BOOST_AUTO( i, neighborhood.begin() );
    BOOST_AUTO( e, neighborhood.end() );

    for(; i!=e ; ++i){
        size_t deg=boost::degree(*i,g);
        // assert(deg>0);
#ifdef NDEBUG
        if(deg)
#endif
        {
        d.reg(*i,deg);
        }
    }
}

namespace impl {

//Construct a tree decomposition T of G using the elimination ordering
//obtained by the minimum-degree heuristic. Ignore isolated vertices.
template <typename G_t, typename T_t/*=typename noboost::treedec_chooser<G_t>::type c++11)*/>
size_t /*FIXME*/ minDegree_decomp(G_t &G, T_t *T/*=NULL (need c++11)*/)
{
    typedef typename noboost::treedec_chooser<G_t>::value_type my_vd;
    typedef typename boost::graph_traits<G_t>::vertex_iterator vertex_iterator;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::adjacency_iterator adjacency_iterator;
    typedef typename misc::deg_chooser<G_t>::type degs_type;
    std::vector<typename noboost::treedec_traits<T_t>::bag_type> bags(boost::num_vertices(G));
    std::vector<my_vd> elim_vertices(boost::num_vertices(G));

    degs_type degs(G);

#ifndef NDEBUG
    noboost::check(G);
#endif

    unsigned int i = 0;
    unsigned min_ntd = 1; // minimum nontrivial vertex degree
    unsigned num_vert = boost::num_vertices(G);
    while(boost::num_edges(G) > 0){
        assert(min_ntd != num_vert);
        //Search a minimum degree vertex.
        vertex_iterator vIt, vEnd;
        boost::tie(vIt, vEnd) = boost::vertices(G);

        // recompute ntd can only increase from here
        vertex_descriptor c;
        if(min_ntd>1){
            --min_ntd;
        }
        boost::tie(c, min_ntd) = degs.pick_min(min_ntd, num_vert);
        assert(min_ntd == boost::degree(c,G));

        adjacency_iterator I, E;
        for(boost::tie(I, E) = boost::adjacent_vertices(c, G); I!=E; ++I){
                assert(*I!=c); // no self loops...
                vertex_descriptor i=*I;
                assert(boost::out_degree(*I,G)>0);
                degs.unlink(i);
        }

        noboost::make_clique_and_hijack(c, G, (detail::degree_mod<G_t>*)NULL, bags[i]);
#ifndef NDEBUG // safety net.
        noboost::check(G);
#endif
        redegree(NULL, G, bags[i], degs);

        degs.unlink(c, min_ntd);
        assert(boost::out_degree(c, G)==0);

        // this is an unsigned for balu...
        elim_vertices[i] = noboost::get_vd(G, c);

        if(min_ntd>1){
            // if a neigbor of c was already connected to the others,
            // min_ntd possibly decreases by one.
            --min_ntd;
        }
        assert(boost::out_degree(c, G)==0);

#ifndef NDEBUG //redundant, for checking only
        degs.reg(c,0);
#endif
        ++i; // number of nodes in tree decomposition tree
    }
    assert(boost::num_edges(G)==0);


    if(T){ untested();
        for(; i > 0; i--){
            typename noboost::treedec_chooser<G_t>::value_type e=elim_vertices[i-1];
            glue_bag(bags[i-1], e, *T);
        }
        return 0;
    }else{
        incomplete();
        return 17;
    }
}

} // impl

// FIXME: remove (underscore prefix: private)
template <typename G_t, typename T_t>
void _minDegree_decomp(G_t &G, T_t& T)
{
    impl::minDegree_decomp(G, &T);
}

//Constructs a tree decomposition from the elimination ordering obtained by the
//minimum-degree heuristic.
template <typename G_t, typename T_t>
void minDegree_decomp(G_t &G, T_t &T){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    std::vector< boost::tuple<typename noboost::treedec_traits<T_t>::bag_type::value_type,
                              typename noboost::treedec_traits<T_t>::bag_type> > bags;

    Islet(G, bags);
    _minDegree_decomp(G, T);
    treedec::preprocessing_glue_bags(bags, T);
}


//Constructs a tree decomposition from the elimination ordering obtained by the
//fill-in heuristic. Ignores isolated vertices.
template <typename G_t, typename T_t>
void _fillIn_decomp(G_t &G, T_t &T){
    std::vector<typename noboost::treedec_traits<T_t>::bag_type> bags(boost::num_vertices(G));
    std::vector<typename noboost::treedec_traits<T_t>::bag_type::value_type> elim_vertices(boost::num_vertices(G));

    unsigned int i = 0;
    while(boost::num_edges(G) > 0){
        //Search a vertex v such that least edges are missing for making the neighbourhood of v a clique.
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        boost::tie(vIt, vEnd) = boost::vertices(G);
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex = *vIt;

        unsigned int min_fill = UINT_MAX;
        for(; vIt != vEnd; vIt++){
            if(boost::out_degree(*vIt, G) == 0){
                continue;
            }

            unsigned int current_fill = 0;

            typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
            for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(*vIt, G); nIt1 != nEnd; nIt1++){
                nIt2 = nIt1;
                nIt2++;
                for(; nIt2 != nEnd; nIt2++){
                    if(!boost::edge(*nIt1, *nIt2, G).second){
                        current_fill++;
                    }
                }
            }

            if(current_fill < min_fill){
                min_fill = current_fill;
                min_vertex = *vIt;
                if(current_fill == 0){
                    break;
                }
            }
        }

        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, G); nIt != nEnd; nIt++){
            bags[i].insert((typename noboost::treedec_traits<T_t>::bag_type::value_type) *nIt);
        }

        noboost::make_clique(boost::adjacent_vertices(min_vertex, G), G); //replace this with make_clique_and_hijack?!

        elim_vertices[i++] = min_vertex;

        boost::clear_vertex(min_vertex, G);
    }

    for(; i > 0; i--){
        treedec::glue_bag(bags[i-1], elim_vertices[i-1], T);
    }
}

//Constructs a tree decomposition from the elimination ordering obtained by the
//fill-in heuristic.
template <typename G_t, typename T_t>
void fillIn_decomp(G_t &G, T_t &T){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    std::vector< boost::tuple<typename noboost::treedec_traits<T_t>::bag_type::value_type,
                              typename noboost::treedec_traits<T_t>::bag_type> > bags;

    treedec::Islet(G, bags);
    treedec::_fillIn_decomp(G, T);
    treedec::preprocessing_glue_bags(bags, T);
}


//Computes an elimination ordering according to the minDegree heuristic (version used for postprocessing algorithms).
template<typename G_t>
void _minDegree_ordering(G_t G, std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering,
       std::vector<bool> &visited)
{
    unsigned int i = 0;
    while(true){
        //Search a minimum degree vertex.
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;

        unsigned int min_degree = UINT_MAX;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int pos = noboost::get_pos(*vIt, G);
            if(visited[pos]){
                continue;
            }

            unsigned int degree = boost::out_degree(*vIt, G);
            if(degree < min_degree){
                min_degree  = degree;
                min_vertex = *vIt;
            }
        }

        if(min_degree == UINT_MAX){
            return;
        }

        noboost::make_clique(boost::adjacent_vertices(min_vertex, G), G);

        elim_ordering[i++] = min_vertex;
        visited[noboost::get_pos(min_vertex, G)] = true;

        boost::clear_vertex(min_vertex, G);
    }
}

//Computes an elimination ordering according to minDegree heuristic.
template<typename G_t>
void minDegree_ordering(G_t G,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering,
      bool ignore_isolated_vertices=false)
{
    std::vector<bool> visited(boost::num_vertices(G), false);
    unsigned int n_ = 0;

    //Mark isolated vertices as already visited.
    if(ignore_isolated_vertices){
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            if(boost::degree(*vIt, G) == 0){
                unsigned int pos = noboost::get_pos(*vIt, G);
                visited[pos] = true;
            }
            else{ n_++; }
        }
    }
    else{ n_ = boost::num_vertices(G); }

    elim_ordering.resize(n_);

    _minDegree_ordering(G, elim_ordering, visited);
}


//Computes an elimination ordering according to fillIn heuristic (version used for postprocessing algorithms).
template<typename G_t>
void _fillIn_ordering(G_t &G,
       std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering,
       std::vector<bool> &visited)
{
    unsigned int i = 0;
    while(true){
        //Search a vertex v such that least edges are missing for making the neighbourhood of v a clique.
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;

        unsigned int min_fill = UINT_MAX;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int pos = noboost::get_pos(*vIt, G);
            if(visited[pos]){
                continue;
            }

            unsigned int current_fill = 0;

            for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(*vIt, G); nIt1 != nEnd; nIt1++){
                nIt2 = nIt1;
                nIt2++;
                for(; nIt2 != nEnd; nIt2++){
                    if(!boost::edge(*nIt1, *nIt2, G).second){
                        current_fill++;
                    }
                }
            }

            if(current_fill < min_fill){
                min_fill = current_fill;
                min_vertex = *vIt;
                if(current_fill == 0){
                    break;
                }
            }
        }

        if(min_fill == UINT_MAX){
            return;
        }

        noboost::make_clique(boost::adjacent_vertices(min_vertex, G), G);

        elim_ordering[i++] = min_vertex;
        visited[noboost::get_pos(min_vertex, G)] = true;

        boost::clear_vertex(min_vertex, G);
    }
}

//Computes an elimination ordering according to fillIn heuristic.
template<typename G_t>
void fillIn_ordering(G_t G,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering,
      bool ignore_isolated_vertices=false)
{
    std::vector<bool> visited(boost::num_vertices(G), false);
    unsigned int n_ = 0;

    //Mark isolated vertices as already visited.
    if(ignore_isolated_vertices){
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            if(boost::degree(*vIt, G) == 0){
                unsigned int pos = noboost::get_pos(*vIt, G);
                visited[pos] = true;
            }
            else{ n_++; }
        }
    }
    else{ n_ = boost::num_vertices(G); }

    elim_ordering.resize(n_);

    _fillIn_ordering(G, elim_ordering, visited);
}



template <typename G_t>
int get_width_of_elimination_ordering(G_t &G,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering)
{
    int width = -1;

    for(unsigned int i = 0; i < elimination_ordering.size(); i++){
        unsigned int deg = noboost::eliminate_vertex(elimination_ordering[i], G);

        width = (width > (int)deg)? width : (int)deg;
    }

    return width;
}


template <typename G_t, typename T_t>
void _ordering_to_treedec(G_t &G,
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering,
    T_t &T, unsigned int idx, bool ignore_isolated_vertices)
{
    if(idx == elimination_ordering.size()){
        return;
    }

    if(ignore_isolated_vertices && boost::degree(elimination_ordering[idx], G) == 0){
        treedec::_ordering_to_treedec(G, elimination_ordering, T, idx+1, ignore_isolated_vertices);
        return;
    }

    typename noboost::treedec_traits<T_t>::bag_type bag;
    noboost::fetch_neighbourhood(bag, boost::adjacent_vertices(elimination_ordering[idx], G), G);

    noboost::make_clique(boost::adjacent_vertices(elimination_ordering[idx], G), G);

    boost::clear_vertex(elimination_ordering[idx], G);

    treedec::_ordering_to_treedec(G, elimination_ordering, T, idx+1, ignore_isolated_vertices);

    treedec::glue_bag(bag, elimination_ordering[idx], T);
}

template <typename G_t, typename T_t>
void ordering_to_treedec(G_t &G,
                         std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering,
                         T_t &T, bool ignore_isolated_vertices=false)
{
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    treedec::_ordering_to_treedec(G, elimination_ordering, T, 0, ignore_isolated_vertices);
}


template <typename G_t, typename T_t>
void _treedec_to_ordering(T_t &T,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering)
{
    bool leaf_found = false;

    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    typename boost::graph_traits<T_t>::vertex_descriptor leaf, parent;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        if(boost::out_degree(*tIt, T) <= 1 && !noboost::bag(*tIt, T).empty()){
            leaf = *tIt;
            leaf_found = true;
            break;
        }
    }

    if(leaf_found){
        typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
        boost::tie(nIt, nEnd) = boost::adjacent_vertices(leaf, T);
        parent = *nIt;

        typename noboost::treedec_traits<T_t>::bag_type difference;

        if(boost::out_degree(leaf, T) == 1){
            if(!std::includes(noboost::bag(parent, T).begin(),
                              noboost::bag(parent, T).end(),
                              noboost::bag(leaf, T).begin(),
                              noboost::bag(leaf, T).end()))
            {
                std::set_difference(noboost::bag(leaf, T).begin(),
                                    noboost::bag(leaf, T).end(),
                                    noboost::bag(parent, T).begin(),
                                    noboost::bag(parent, T).end(),
                                    std::inserter(difference, difference.begin()));
            }
            boost::clear_vertex(leaf, T);
        }
        else{
            difference = MOVE(noboost::bag(leaf, T));
        }

        for(typename noboost::treedec_traits<T_t>::bag_type::iterator sIt = difference.begin();
            sIt != difference.end(); sIt++)
        {
            elimination_ordering.push_back(*sIt);
        }

        noboost::bag(leaf, T).clear();

        _treedec_to_ordering<G_t, T_t>(T, elimination_ordering);
    }
}

template <typename G_t, typename T_t>
void treedec_to_ordering(T_t &T,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering)
{
    if(boost::num_vertices(T) == 0){
        return;
    }

    _treedec_to_ordering<G_t, T_t>(T, elimination_ordering);
}

//Make G a filled graph according to the provided elimination_ordering. Stores the cliques in C and the additional
//edges in F.
template <typename G_t>
void make_filled_graph(G_t &G,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering,
      std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > &C,
      std::vector<std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > > &F)
{
    C.resize(elim_ordering.size());
    F.resize(elim_ordering.size());

    std::vector<bool> visited(boost::num_vertices(G), false);

    for(unsigned int i = 0; i < elim_ordering.size(); i++){
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N_i, E_i;
        C[i].insert(elim_ordering[i]);

        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(elim_ordering[i], G); nIt != nEnd; nIt++){
            unsigned int pos = noboost::get_pos(*nIt, G);
            if(!visited[pos]){
                C[i].insert(*nIt);
            }
        }

        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt1 =
            C[i].begin(); sIt1 != C[i].end(); sIt1++)
        {
            typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt2 = sIt1;
            sIt2++;
            for(; sIt2 != C[i].end(); sIt2++){
                if(!boost::edge(*sIt1, *sIt2, G).second){
                    typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> edge(2);
                    edge[0] = *sIt1;
                    edge[1] = *sIt2;
                    F[i].push_back(edge);
                    boost::add_edge(*sIt1, *sIt2, G);
                }
            }
        }

        unsigned int pos = noboost::get_pos(elim_ordering[i], G);
        visited[pos] = true;
    }
}

template <typename G_t>
void LEX_M_fill_in(G_t &G,
     std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > &fill_in_edges)
{
    unsigned int max = boost::num_vertices(G);
    std::vector<bool> visited(max+1);
    std::vector<float> label(max+1);
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> alpha_inv(max+1);
    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > reached_i(max+1);

    //Initializing.
    unsigned int i = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int pos = noboost::get_pos(*vIt, G);
        label[pos] = 1.0;
        alpha_inv[i++] = 0;
        visited[pos] = false;
    }

    unsigned int k = 1;

    for(int i = boost::num_vertices(G)-1; i >= 0; i--){
        typename boost::graph_traits<G_t>::vertex_descriptor v = *vIt;
        unsigned int max = 0;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int pos = noboost::get_pos(*vIt, G);
            if(alpha_inv[pos] == 0){
                if(label[pos] > max){
                    max = (unsigned int) label[pos];
                    v = *vIt;
                }
            }
        }
        unsigned int pos = noboost::get_pos(v, G);
        visited[pos] = true;
        alpha_inv[pos] = i+1;

        for(unsigned int j = 0; j < k; j++){
            reached_i[j].clear();
        }

        for(unsigned int j = 0; j < alpha_inv.size(); j++){
            if(alpha_inv[j] == 0){
                visited[j] = false;
            }
        }

        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; nIt++){
            unsigned int posn = noboost::get_pos(*nIt, G);
            if(alpha_inv[posn] == 0){
                reached_i[(int)label[posn]-1].push_back(*nIt);
                visited[posn] = true;
                label[posn] += 0.5;
            }
        }

        for(unsigned int j = 0; j < k; j++){
            while(reached_i[j].size() != 0){
                typename boost::graph_traits<G_t>::vertex_descriptor w = reached_i[j].back();

                reached_i[j].pop_back();
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(w, G); nIt != nEnd; nIt++){
                    unsigned int posn = noboost::get_pos(*nIt, G);
                    if(visited[posn]){
                        continue;
                    }

                    visited[posn] = true;
                    if((unsigned int)label[posn]-1 > j){
                        reached_i[(int)label[posn]].push_back(*nIt);
                        label[posn] += 0.5;
                        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> edge(2);
                        edge[0] = v;
                        edge[1] = *nIt;
                        fill_in_edges.push_back(edge);
                    }
                    else{
                        reached_i[j].push_back(*nIt);
                    }
                }
            }
        }

        for(unsigned int j = 0; j < label.size(); j++){
            label[j] = (float)roundf(label[j]);
            k = (k > (unsigned int)label[j])? k : (unsigned int)label[j];
        }
    }
}

template <typename G_t>
void LEX_M_minimal_ordering(G_t &G,
     typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &alpha)
{
    unsigned int max = boost::num_vertices(G);
    alpha.resize(boost::num_vertices(G));
    std::vector<bool> visited(max+1);
    std::vector<float> label(max+1);
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> alpha_inv(max+1);
    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > reached_i(max+1);

    unsigned int i = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int pos = noboost::get_pos(*vIt, G);
        label[pos] = 1.0;
        alpha_inv[i++] = 0;
        visited[pos] = false;
    }

    unsigned int k = 1;

    for(int i = boost::num_vertices(G)-1; i >= 0; i--){
        typename boost::graph_traits<G_t>::vertex_descriptor v=*vEnd;
        unsigned int max = 0;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int pos = noboost::get_pos(*vIt, G);
            if(alpha_inv[pos] == 0){
                if((unsigned int)label[pos] > max){
                    max = (unsigned int) label[pos];
                    v = *vIt;
                }
            }
        }
        unsigned int posv = noboost::get_pos(v, G);
        visited[posv] = true;
        alpha[i] = v;
        alpha_inv[posv] = i+1;

        for(unsigned int j = 0; j < k; j++){
            reached_i[j].clear();
        }

        for(unsigned int j = 0; j < alpha_inv.size(); j++){
            if(alpha_inv[j] == 0){
                visited[j] = false;
            }
        }

        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; nIt++){
            unsigned int posn = noboost::get_pos(*nIt, G);
            if(alpha_inv[posn] == 0){
                reached_i[(int)label[posn]-1].push_back(posn);
                visited[posn] = true;
                label[posn] += 0.5;
            }
        }

        for(unsigned int j = 0; j < k; j++){
            while(reached_i[j].size() != 0){
                typename boost::graph_traits<G_t>::vertex_descriptor w = reached_i[j].back();

                reached_i[j].pop_back();
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(w, G); nIt != nEnd; nIt++){
                    unsigned int posn = noboost::get_pos(*nIt, G);
                    if(visited[posn]){
                        continue;
                    }

                    visited[posn] = true;
                    if((unsigned int)label[posn]-1 > j){
                        reached_i[(int)label[posn]].push_back(*nIt);
                        label[posn] += 0.5;
                    }
                    else{
                        reached_i[j].push_back(*nIt);
                    }
                }
            }
        }

        for(unsigned int j = 0; j < label.size(); j++){
            label[j] = (float)roundf(label[j]);
            k = (k > (unsigned int)label[j])? k : (unsigned int)label[j];
        }
    }
}

} //namespace treedec

#endif //TD_ELIMINATION_ORDERING

// vim:ts=8:sw=4:et
