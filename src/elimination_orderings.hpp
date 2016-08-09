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
 * Offers functionality to compute tree decompositions of graphs
 * according to various heuristic, and functions which convert
 * tree decompositions to elimination orderings and vice versa.
 * Also the LEX-M algorithm is included in this header
 *
 * These functions are most likely to be interesting for outside use:
 *
 * - void minDegree_decomp(G_t &G, T_t &T)
 * - void fillIn_decomp(G_t &G, T_t &T)
 * - void minDegree_ordering(G_t &G,
 *         typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering)
 * - void fillIn_ordering(G_t& G,
 *         typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering)
 * - void ordering_to_treedec(G_t &G,
 *         typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering, T_t &T)
 * - void treedec_to_ordering<G_t, T_t>(T_t &T,
 *         typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering)
 * - void LEX_M_minimal_ordering(G_t &G,
 *         typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering)
 *
*/

#ifndef TD_ELIMINATION_ORDERINGS
#define TD_ELIMINATION_ORDERINGS

#include <cmath>
#include <climits>

#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

#include <boost/graph/minimum_degree_ordering.hpp>

#include "trace.hpp"
#include "preprocessing.hpp"
#include "simple_graph_algos.hpp"
#include "misc.hpp"

#include "graph.hpp"
#include "fill.hpp"
#include "platform.hpp"

namespace treedec{

//register a 2-neigborhood to FILL
template<class U, class G_t, class B, class F>
std::pair<typename boost::graph_traits<G_t>::vertex_descriptor,unsigned>
        refill(U, G_t &G, B& range, F& fill)
{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;

    while(!range.empty()){
        BOOST_AUTO(b, *range.begin());
        range.erase(range.begin());

        unsigned current_fill=treedec::count_missing_edges(b, G);

        if(!current_fill){
            return std::make_pair(b,0);
        }else{
            fill.reg(b, current_fill);
        }

    }
    return std::make_pair(vertex_descriptor(),-1);
}

namespace detail{
// update vertices adjacent to neighbors that have just been connected.
// outedgee to center node has ideally been removed already.
template<typename G_t>
struct fill_update_cb : public graph_callback<G_t>{
    typedef typename boost::graph_traits<G_t>::edge_descriptor edge_descriptor;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename fill_chooser<G_t>::type fill_type;

    fill_update_cb(fill_type* d, G_t const& g) :
        _fill(d), G(g){}

    void operator()(vertex_descriptor v)
    {
        _fill->q_eval(v);
    }
    void operator()(edge_descriptor edg)
    {
        assert(boost::source(edg, G) < boost::target(edg, G));
        // e has just been inserted.
        BOOST_AUTO(cni, common_out_edges(boost::source(edg, G), boost::target(edg, G), G));
        BOOST_AUTO(i, cni.first);
        BOOST_AUTO(e, cni.second);
        for(; i!=e; ++i){
            assert(*i != boost::source(edg, G));
            assert(*i != boost::target(edg, G));
//            no. maybe theres only half an edge.
//            assert(boost::edge(boost::source(edg, G), *i, G).second);
//            assert(boost::edge(boost::target(edg, G), *i, G).second);

            // BUG: *i might be within 1-neighborhood.
            _fill->q_decrement(*i);
        }
    }
private:
    fill_type* _fill;
    G_t const& G;
};

}// detail

namespace detail{

template <typename G_t, typename T_t, typename B_t, typename O_t>
void skeletal_to_treedec(G_t &G, T_t &T, B_t &B, O_t &O, unsigned n_)
{
    typedef typename treedec_traits<T_t>::bag_type bag_type;

    std::vector<unsigned int> inv_O(boost::num_vertices(G), n_+1);
    for(unsigned u = 0; u < n_; u++){
        typename treedec_chooser<G_t>::value_type e=O[u];
        unsigned pos = get_pos(e, G);
        inv_O[pos] = u;
    }

    //Bag for the u-th elimination vertex will be stored in T[u].
    for(unsigned u = 0; u < n_; u++){
        boost::add_vertex(T);
    }


    //Since we made the neighbourhood N of the u-th vertex a clique,
    //the bag of the neighbour of this vertex with lowest elimination index
    //will have N as a subset.
    unsigned max = n_-1u;
    for(unsigned u = 0; u < max; u++){
        unsigned min_index = max; //note: if there's an empty bag, we can glue
                                  //it toghether with an arbitrary bag.
        for(typename bag_type::iterator bIt = B[u].begin(); bIt != B[u].end(); bIt++){
           unsigned pos = get_pos(*bIt, G);
           unsigned index = inv_O[pos];
           if(index < min_index){
               min_index = index;
           }
        }
        //(min_index, u) will lead to a connected directed graph, if G_t is
        //directed.
        boost::add_edge(min_index, u, T);
    }

    //Bag for the u-th elimination vertex will be stored in T[u].
    for(unsigned u = 0; u < n_; u++){
        bag(u, T) = MOVE(B[u]);
        bag(u, T).insert(O[u]); //printer variant without this inserting?
    }
}

} //namespace detail


namespace impl{

//Construct a tree decomposition T of G using the elimination ordering
//obtained by the minimum-degree heuristic. Ignore isolated vertices.
#if __cplusplus >= 201103L
template <typename G_t, typename T_t=typename treedec_chooser<G_t>::type, typename O_t>
typename boost::graph_traits<G_t>::vertices_size_type
   minDegree_decomp(G_t &G, T_t *T=NULL, O_t *O=NULL, unsigned ub=UINT_MAX)
#else
template <typename G_t, typename T_t, typename O_t>
typename boost::graph_traits<G_t>::vertices_size_type
   minDegree_decomp(G_t &G, T_t *T, O_t *O, unsigned ub=UINT_MAX)
#endif
{
    typedef typename treedec_chooser<G_t>::value_type my_vd;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::adjacency_iterator adjacency_iterator;
    typedef typename deg_chooser<G_t>::type degs_type;
    typedef typename treedec_traits<T_t>::bag_type bag_type;
    std::vector<bag_type> bags;
    bag_type bag_i;
    bag_type* bags_i=&bag_i;
    std::vector<my_vd> elim_vertices;

    degs_type degs(G);
    const degs_type& cdegs(degs);

    typename boost::graph_traits<G_t>::vertices_size_type num_vert=boost::num_vertices(G);

    if(num_vert == 0){
        if(T){
            boost::add_vertex(*T);
            return 0;
        }
        else{
            return 0; //should be -1 ..
        }
    }

#ifndef NDEBUG
    check(G);
#endif

    unsigned int i = 0;
    unsigned min_ntd = 1; // minimum nontrivial vertex degree
    unsigned upper_bound = 0; // computed, if T

    if(T){
        bags.resize(num_vert);
        elim_vertices.resize(num_vert);

        BOOST_AUTO(it, cdegs[0].begin());
        for(; it!=cdegs[0].end(); ++it){
            elim_vertices[i++] = get_vd(G, *it);
        }
    }

    while(boost::num_edges(G) > 0){
        INTERRUPTION_POINT;
        assert(min_ntd != num_vert);

        // recompute ntd can only increase from here
        vertex_descriptor c;
        if(min_ntd>1){
            --min_ntd;
        }
        boost::tie(c, min_ntd) = degs.pick_min(min_ntd, num_vert);
        assert(min_ntd == boost::degree(c,G));

        //Abort if the width of this decomposition would be larger than 'ub'.
        if(min_ntd > ub){
            T->clear();
            throw exception_unsuccessful();
        }

        if(T){
            bags_i = &bags[i];
            elim_vertices[i] = get_vd(G, c);
        }
        else if(min_ntd > upper_bound){
            upper_bound = min_ntd;
        }
        else{}

        ++i; // number of nodes in tree decomposition tree

        adjacency_iterator I, E;
        for(boost::tie(I, E) = boost::adjacent_vertices(c, G); I!=E; ++I){
            assert(*I!=c); // no self loops...
            vertex_descriptor w=*I;
            degs.unlink(w);
        }

        make_clique_and_detach(c, G, *bags_i);
#ifndef NDEBUG // safety net.
        check(G);
#endif

        redegree(NULL, G, *bags_i, degs);
        if(!T){
            bags_i->clear();
        }

        degs.unlink(c, min_ntd);

        assert(boost::degree(c, G)==0);

        degs.flush();
    }
    assert(boost::num_edges(G)==0);

    BOOST_AUTO(it, cdegs[0].begin());
    for(; it!=cdegs[0].end(); ++it){
        elim_vertices[i++] = get_vd(G, *it);
    }

    //Build a treedecomposition.
    if(T){
        assert(i == num_vert);
        treedec::detail::skeletal_to_treedec(G, *T, bags, elim_vertices, num_vert);

        return 0;
    }
    else{
        return upper_bound;
    }
}

template <typename G_t, typename T_t>
typename boost::graph_traits<G_t>::vertices_size_type
  minDegree_decomp(G_t &G, T_t *T, unsigned ub=UINT_MAX)
{
    return minDegree_decomp(G, T, (typename std::vector<typename treedec_chooser<G_t>::value_type>*)NULL);

}

#if __cplusplus < 201103L
template <typename G_t>
typename boost::graph_traits<G_t>::vertices_size_type
   minDegree_decomp(G_t &G)
{
    return minDegree_decomp(G, (typename treedec_chooser<G_t>::type*)NULL);
}
#endif

} // namespace impl

//Constructs a tree decomposition from the elimination ordering obtained by the
//minimum-degree heuristic.
template <typename G_t, typename T_t>
void minDegree_decomp(G_t &G, T_t &T, unsigned ub=UINT_MAX){
    impl::minDegree_decomp(G, &T);
}

namespace impl{

template <typename G_t>
void boost_minDegree_ordering(G_t &G, std::vector<int> &O){
    typedef typename boost::graph_traits<G_t>::edges_size_type edges_size_type;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;

    vertices_size_type n = boost::num_vertices(G);
    edges_size_type e = boost::num_edges(G);


    O.resize(n);

    unsigned i = 0;
    // boost does not like cliques.
    if((n*(n-1u)) == boost::num_edges(G)){
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            O[i++] = *vIt;
        }
    }

    std::vector<int> inverse_perm(n, 0);
    std::vector<int> supernode_sizes(n, 1);
    typename boost::property_map<G_t, boost::vertex_index_t>::type id = boost::get(boost::vertex_index, G);
    std::vector<int> degree(n, 0);

    boost::minimum_degree_ordering
             (G,
              boost::make_iterator_property_map(&degree[0], id, degree[0]),
              &inverse_perm[0],
              &O[0],
              boost::make_iterator_property_map(&supernode_sizes[0], id, supernode_sizes[0]),
              0,
              id);
}


template <typename G_t, typename T_t>
void endless_minDegree_decomp(G_t &G, T_t &T){
    unsigned best = UINT_MAX;
    G_t Gbak(G);
    while(true){
        try{
            T_t U;
            best = treedec::minDegree_decomp(G, U, best);
            T = MOVE(U);
        }
        catch(...){
        }
        boost::copy_graph(Gbak, G);
    }
}

} //namespace impl

namespace impl{

//Construct a tree decomposition from the elimination ordering obtained by the
//fill-in heuristic. Ignores isolated vertices.
//return the treewidth.
#if __cplusplus >= 201103L
template <typename G_t, typename T_t=typename treedec_chooser<G_t>::type, typename O_t=
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor>
    >
typename boost::graph_traits<G_t>::vertices_size_type
   fillIn_decomp(G_t &G, T_t *T=NULL, O_t *O=NULL, unsigned ub = UINT_MAX)
#else
template <typename G_t, typename T_t, typename O_t>
typename boost::graph_traits<G_t>::vertices_size_type
   fillIn_decomp(G_t &G, T_t *T, O_t *O, unsigned ub = UINT_MAX)
#endif
{
    typedef typename treedec_chooser<G_t>::value_type my_vd;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename fill_chooser<G_t>::type fill_type;
    typedef typename treedec_traits<T_t>::bag_type bag_type;
    std::vector<bag_type> bags;
    bag_type bag_i;
    bag_type* bags_i = &bag_i;
    std::vector<my_vd> elim_vertices;

    std::set<vertex_descriptor> refill_q;

    typename boost::graph_traits<G_t>::vertices_size_type num_vert = boost::num_vertices(G);
    if(T){
        bags.resize(num_vert);
        elim_vertices.resize(num_vert);
    }

    fill_type fill(G);
    detail::fill_update_cb<G_t> cb(&fill, G);

    unsigned int i = 0;
    unsigned int min_fill = -1;
    unsigned int upper_bound = 0; // computed, if T

    vertex_descriptor v;

    while(boost::num_edges(G) > 0){
        INTERRUPTION_POINT;
        //Find a vertex v such that least edges are missing for making the
        //neighbourhood of v a clique.
        //
        fill.check();
        boost::tie(v, min_fill) = fill.pick_min(0, -1, true);
        fill.check();
        assert(is_valid(v,G));

        BOOST_AUTO(deg, boost::degree(v, G));

        //Abort if the width of this decomposition would be larger than 'ub'.
        if(deg > ub){
            T->clear();
            throw exception_unsuccessful();
        }

        if(T){
            assert(i<bags.size());
            bags_i = &bags[i];
            elim_vertices[i] = get_vd(G, v);
        }else{untested();
        }

        fill.mark_neighbors(v, min_fill);

        assert(!bags_i->size());


#ifndef NDEBUG
        size_t newedges = make_clique_and_detach(v, G, *bags_i, &cb);
        if(newedges == min_fill){
        }else{ untested();
            assert(false); // for now.
            // something is terribly wrong.
            // or some extra-heuristics is active
        }
#else
        make_clique_and_detach(v, G, *bags_i, &cb);
#endif

        fill.unmark_neighbours(*bags_i);


        if(!T){
            bags_i->clear();
            if(deg > upper_bound){
                upper_bound = deg;
            }
        }

        assert(boost::degree(v, G)==0);
        ++i; // number of nodes in tree decomposition tree
    }

    if(T){
        treedec::detail::skeletal_to_treedec(G, *T, bags, elim_vertices, i);

        return 0;
    }
    else{
        return upper_bound;
    }
}

template <typename G_t, typename T_t>
void fillIn_decomp(G_t &G, T_t *T, unsigned ub = UINT_MAX)
{
    fillIn_decomp(G, T, (typename std::vector<typename treedec_chooser<G_t>::value_type>*)NULL, ub);
}

#if __cplusplus < 201103L
template <typename G_t>
typename boost::graph_traits<G_t>::vertices_size_type
   fillIn_decomp(G_t &G)
{
    return fillIn_decomp(G, (typename treedec_chooser<G_t>::type*)NULL,
                            (typename std::vector<typename treedec_chooser<G_t>::value_type>*)NULL);
}
#endif

} //namespace impl


//Construct a tree decomposition from the elimination ordering obtained by the
//fill-in heuristic.
template <typename G_t, typename T_t>
void fillIn_decomp(G_t &G, T_t &T, unsigned ub = UINT_MAX){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    std::vector< boost::tuple<typename treedec_traits<T_t>::bag_type::value_type,
                              typename treedec_traits<T_t>::bag_type> > bags;

    treedec::impl::Islet(G, bags);
    impl::fillIn_decomp(G, &T, ub);
    treedec::glue_bags(bags, T);
}

namespace detail{

//Compute an elimination ordering according to the minDegree heuristic
//(version used for postprocessing algorithms).
// TODO use impl::minDegree (how?)
template<typename G_t>
void minDegree_ordering(G_t G,
       std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering,
       std::vector<bool> &visited)
{
    unsigned int i = 0;
    while(true){
        INTERRUPTION_POINT;
        //Search a minimum degree vertex.
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;

        unsigned int min_degree = UINT_MAX;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int pos = get_pos(*vIt, G);
            if(visited[pos]){
                continue;
            }

            unsigned int degree = boost::degree(*vIt, G);
            if(degree < min_degree){
                min_degree  = degree;
                min_vertex = *vIt;
            }
        }

        if(min_degree == UINT_MAX){
            return;
        }

        make_clique(boost::adjacent_vertices(min_vertex, G), G, (graph_callback<G_t>*)NULL);

        elim_ordering[i++] = min_vertex;
        visited[get_pos(min_vertex, G)] = true;

        boost::clear_vertex(min_vertex, G);
    }
}

} //namespace detail

//Computes an elimination ordering according to minDegree heuristic.
template<typename G_t>
void minDegree_ordering(G_t& G,
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
                unsigned int pos = get_pos(*vIt, G);
                visited[pos] = true;
            }
            else{ n_++; }
        }
    }
    else{ n_ = boost::num_vertices(G); }

    elim_ordering.resize(n_);

    detail::minDegree_ordering(G, elim_ordering, visited);
}

namespace detail{

//Compute an elimination ordering according to fillIn heuristic (version used
//for postprocessing algorithms).
template<typename G_t>
void fillIn_ordering(G_t &G,
       std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering,
       std::vector<bool> &visited)
{
    unsigned int i = 0;
    while(true){
        //Search a vertex v such that least edges are missing for making the
        //neighbourhood of v a clique.

        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;

        unsigned int min_fill = UINT_MAX;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int pos = get_pos(*vIt, G);
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

        make_clique(boost::adjacent_vertices(min_vertex, G), G);

        elim_ordering[i++] = min_vertex;
        visited[get_pos(min_vertex, G)] = true;

        boost::clear_vertex(min_vertex, G);
    }
}

} //detail

//Computes an elimination ordering according to fillIn heuristic.
template<typename G_t>
void fillIn_ordering(G_t& G,
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
                unsigned int pos = get_pos(*vIt, G);
                visited[pos] = true;
            }
            else{
                n_++;
            }
        }
    }
    else{ n_ = boost::num_vertices(G); }

    elim_ordering.resize(n_);

    detail::fillIn_ordering(G, elim_ordering, visited);
}



template <typename G_t>
int get_width_of_elimination_ordering(G_t &G,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering)
{
    int width = -1;

    for(unsigned int i = 0; i < elimination_ordering.size(); i++){
        unsigned deg=boost::degree(elimination_ordering[i], G);

        typename outedge_set<G_t>::type xbag;
        treedec::make_clique_and_detach(elimination_ordering[i], G, xbag);
        xbag.clear(); // provide interface with clear included? (not urgent)

        width = (width > (int)deg)? width : (int)deg;
    }

    return width;
}

namespace impl{

template <typename G_t, typename T_t>
void ordering_to_treedec(G_t &G,
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &O,
    T_t &T)
{
    unsigned n = O.size();

    typename std::vector<typename treedec_traits<T_t>::bag_type> bags(n);

    for(unsigned int i = 0; i < O.size(); i++){
        make_clique_and_detach(O[i], G, bags[i]);
    }

    treedec::detail::skeletal_to_treedec(G, T, bags, O, n);
}

} //namespace impl

template <typename G_t, typename T_t>
void ordering_to_treedec(G_t &G,
                         std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &O,
                         T_t &T)
{
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    treedec::impl::ordering_to_treedec(G, O, T);
}

template <typename G_t, typename T_t>
void ordering_to_treedec(G_t &G, std::vector<int> &O, T_t &T)
{
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> O_(O.size());
    for(unsigned int i = 0; i < O.size(); i++){ O_[i] = O[i]; }

    ordering_to_treedec(G, O_, T);
}

namespace draft{

// TODO: tries hard to add edges at the end. does that make sense?
// TODO: what are the preconditions?!
// TODO: can order be an input range?
template <typename G_t, typename O_t, class T_t>
void vec_ordering_to_tree(G_t &G, O_t &O, T_t& t, O_t* io=NULL)
{
    size_t num_vert = boost::num_vertices(G);
    assert(num_vert = O.size());
    O_t iOlocal;

    if(io){
        assert(io->size()==num_vert);
    }else{
        iOlocal.resize(num_vert);
        io=&iOlocal;
    }
    O_t& iO=*io;

    size_t invalid=num_vert;
    std::vector<unsigned> edges(num_vert-1u, invalid);
    assert(edges.size()==num_vert-1);
    typedef boost::adjacency_matrix<boost::directedS> bamd;
    bamd bags(num_vert);

    for(unsigned i = 0; i < num_vert; i++){
        iO[O[i]] = i;
    }

    for(unsigned i = 0; i < num_vert; i++){
        BOOST_AUTO(R, boost::adjacent_vertices(O[i], G));
        for(;R.first!=R.second;++R.first) {
            unsigned n_node = *R.first;
            if((unsigned)iO[n_node] > i){
                boost::add_edge(i, n_node, bags);
            }
        }
    }

    for(unsigned i = 0; i < num_vert; i++){
        std::vector<unsigned> N;
        for(unsigned j = 0; j < num_vert; j++){
            if(boost::edge(i, j, bags).second){
                N.push_back(j);
                unsigned iO_n_node = iO[j];
                if(iO_n_node < edges[i]){
                    edges[i] = iO_n_node;
                }
            }
        }

        for(unsigned j = 0; j < N.size(); j++){
            for(unsigned k = 0; k < N.size(); k++){
                if(iO[N[k]] > iO[N[j]]){
                    boost::add_edge(iO[N[j]], N[k], bags);
                    if((unsigned)iO[N[k]] < edges[iO[N[j]]]){
                        edges[iO[N[j]]] = iO[N[k]];
                    }
                }
            }
        }
    }

    for(unsigned i = 0; i < num_vert; i++){
        boost::add_vertex(t);
        typename treedec_traits<T_t>::bag_type& b=bag(i,t);
        duh::push(b,O[i]);
        for(unsigned j = 0; j < num_vert; j++){
            if(boost::edge(i, j, bags).second){
                duh::push(b,j);
             }
         }
     }

    for(unsigned i = 0; i < num_vert-1u; i++){
        assert(edges[i]>i || edges[i]==invalid);
        if(edges[i]!=invalid){
            // normal edge, as computed above.
            boost::add_edge(i,edges[i],t);
        }else if(i+1!=num_vert){
            // edge to next component
            boost::add_edge(i,i+1,t);
        }else{
            // exiting last connected component.
            // ... dont connect
        }
    }

}

} // draft

namespace impl{

template <typename G_t, typename T_t>
void treedec_to_ordering(T_t &T,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &O)
{
    bool leaf_found = false;

    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    typename boost::graph_traits<T_t>::vertex_descriptor leaf, parent;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        if(boost::degree(*tIt, T) <= 1 && !bag(*tIt, T).empty()){
            leaf = *tIt;
            leaf_found = true;
            break;
        }
    }

    if(leaf_found){
        typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
        boost::tie(nIt, nEnd) = boost::adjacent_vertices(leaf, T);
        parent = *nIt;

        typename treedec_traits<T_t>::bag_type difference;

        if(boost::degree(leaf, T) == 1){
            if(!std::includes(bag(parent, T).begin(),
                              bag(parent, T).end(),
                              bag(leaf, T).begin(),
                              bag(leaf, T).end()))
            {
                std::set_difference(bag(leaf, T).begin(),
                                    bag(leaf, T).end(),
                                    bag(parent, T).begin(),
                                    bag(parent, T).end(),
                                    std::inserter(difference, difference.begin()));
            }
            boost::clear_vertex(leaf, T);
        }
        else{
            difference = MOVE(bag(leaf, T));
        }

        for(typename treedec_traits<T_t>::bag_type::iterator sIt = difference.begin();
            sIt != difference.end(); sIt++)
        {
            O.push_back(*sIt);
        }

        bag(leaf, T).clear();

        treedec_to_ordering<G_t, T_t>(T, O);
    }
}

} //namespace impl

template <typename G_t, typename T_t>
void treedec_to_ordering(T_t &T,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering)
{
    if(boost::num_vertices(T) == 0){
        return;
    }

    treedec::impl::treedec_to_ordering<G_t, T_t>(T, elimination_ordering);
}

//Make G a filled graph according to the provided elimination_ordering. Stores
//the cliques in C and the additional edges in F.
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
            unsigned int pos = get_pos(*nIt, G);
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

        unsigned int pos = get_pos(elim_ordering[i], G);
        visited[pos] = true;
    }
}

template <typename G_t>
void LEX_M_fill_in(G_t &G,
     std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > &fill_in_edges)
{
    unsigned int nv = boost::num_vertices(G);
    std::vector<bool> visited(nv);
    std::vector<float> label(nv);
    std::vector<bool> alpha_inv(nv);
    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > reached_i(nv);

    //Initializing.
    unsigned int i = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int pos = get_pos(*vIt, G);
        label[pos] = 1.0;
        alpha_inv[i++] = false;
        visited[pos] = false;
    }

    unsigned int k = 1;

    for(int i = boost::num_vertices(G)-1; i >= 0; i--){
        typename boost::graph_traits<G_t>::vertex_descriptor v = *vIt;
        unsigned int max = 0;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int pos = get_pos(*vIt, G);
            if(!alpha_inv[pos]){
                if(label[pos] > max){
                    max = (unsigned int) label[pos];
                    v = *vIt;
                }
            }
        }
        unsigned int pos = get_pos(v, G);
        visited[pos] = true;
        alpha_inv[pos] = true;

        for(unsigned int j = 0; j < k; j++){
            reached_i[j].clear();
        }

        for(unsigned int j = 0; j < alpha_inv.size(); j++){
            if(!alpha_inv[j]){
                visited[j] = false;
            }
        }

        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; nIt++){
            unsigned int posn = get_pos(*nIt, G);
            if(!alpha_inv[posn]){
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
                    unsigned int posn = get_pos(*nIt, G);
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
    unsigned int nv = boost::num_vertices(G);
    alpha.resize(boost::num_vertices(G));
    std::vector<bool> visited(nv);
    std::vector<float> label(nv);
    std::vector<bool> alpha_inv(nv);
    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > reached_i(nv);

    unsigned int i = 0;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int pos = get_pos(*vIt, G);
        label[pos] = 1.0;
        alpha_inv[i++] = 0;
        visited[pos] = false;
    }

    unsigned int k = 1;

    for(int i = boost::num_vertices(G)-1; i >= 0; i--){
        typename boost::graph_traits<G_t>::vertex_descriptor v=*vEnd;
        unsigned int max = 0;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int pos = get_pos(*vIt, G);
            if(!alpha_inv[pos]){
                if((unsigned int)label[pos] > max){
                    max = (unsigned int) label[pos];
                    v = *vIt;
                }
            }
        }
        unsigned int posv = get_pos(v, G);
        visited[posv] = true;
        alpha[i] = v;
        alpha_inv[posv] = true;

        for(unsigned int j = 0; j < k; j++){
            reached_i[j].clear();
        }

        for(unsigned int j = 0; j < alpha_inv.size(); j++){
            if(!alpha_inv[j]){
                visited[j] = false;
            }
        }

        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; nIt++){
            unsigned int posn = get_pos(*nIt, G);
            if(!alpha_inv[posn]){
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
                    unsigned int posn = get_pos(*nIt, G);
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

#endif //TD_ELIMINATION_ORDERINGS

// vim:ts=8:sw=4:et
