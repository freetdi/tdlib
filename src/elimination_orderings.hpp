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
 * For a graph class G_t, a suitable tree decomposition class T_t and a
 * sequence container O_t, the following functions are meant for outside use.
 * TODO/CHECK: these should be the only ones exported.
 *
 * - void minDegree_decomp(G_t &G, T_t &T)
 * - void fillIn_decomp(G_t &G, T_t &T)
 * - void minDegree_ordering(G_t &G, O_t& elim_ordering)
 * - unsigned boost_minDegree_ordering(G_t &G, O_t &elim_ordering)
 * - unsigned boost_minDegree_ordering(G_t &G, O_t &elim_ordering, O_t &inv_elim_ordering)
 * - void fillIn_ordering(G_t& G, O_t &elim_ordering)
 * - void ordering_to_treedec(G_t &G, O_t &elim_ordering, T_t &T)
 *   BUG: G_t is a template arg, why?
 * - void treedec_to_ordering<G_t, T_t>(T_t &T, O_t& elim_ordering)
 * - void LEX_M_minimal_ordering(G_t &G, O_t& elim_ordering)
 *
*/

#ifndef TD_ELIMINATION_ORDERINGS
#define TD_ELIMINATION_ORDERINGS

#include <cmath>
#include <climits>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>

#include <boost/graph/copy.hpp>
// #include <boost/graph/minimum_degree_ordering.hpp>

#include "trace.hpp"
#include "preprocessing.hpp"
#include "simple_graph_algos.hpp"
#include "misc.hpp"

#include "graph.hpp"
#include "fill.hpp"
#include "platform.hpp"

#include "preprocessing.hpp"
#include "simple_graph_algos.hpp"
#include "misc.hpp"

#ifndef MINIMUM_DEGREE_ORDERING_HPP
# include "minimum_degree_ordering.hpp"
# define HAVE_MINDEGREE_FORK
#endif


namespace treedec{ //

// FIXME: move to detail header.
// or to greedy_heuristics class?
namespace detail{

template <typename G_t, typename T_t, typename B_t, typename O_t>
void skeleton_to_treedec(G_t &G, T_t &T, B_t &B, O_t &O, unsigned n_)
{
    typedef typename treedec_traits<T_t>::bag_type bag_type;

    std::vector<unsigned int> inv_O(boost::num_vertices(G), n_+1);
    for(unsigned u = 0; u < n_; u++){ itested();
        typename treedec_chooser<G_t>::value_type e=O[u];
        unsigned pos = get_pos(e, G);
        inv_O[pos] = u;
    }

    if(n_==0){
        return;
    }

    //Bag for the u-th elimination vertex will be stored in T[u].
    for(unsigned u = 0; u < n_; u++){ itested();
        boost::add_vertex(T);
    }


    //Since we made the neighbourhood N of the u-th vertex a clique,
    //the bag of the neighbour of this vertex with lowest elimination index
    //will have N as a subset.
    unsigned max = n_-1u;
    for(unsigned u = 0; u < max; u++){ itested();
        unsigned min_index = max; //note: if there's an empty bag, we can glue
                                  //it toghether with an arbitrary bag.
        for(typename bag_type::iterator bIt = B[u].begin(); bIt != B[u].end(); bIt++){ itested();
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
    for(unsigned u = 0; u < n_; u++){ itested();
        bag(u, T) = MOVE(B[u]);
        insert(bag(u, T), O[u]); //printer variant without this inserting?
    }
}

} //namespace detail

} // namespace treedec

#include "impl/greedy_heuristic.hpp"

namespace treedec{ //

//Construct a tree decomposition from the elimination ordering obtained by the
//minimum-degree heuristic.
//
template <typename G_t, typename T_t, typename O_t>
typename boost::graph_traits<G_t>::vertices_size_type
  minDegree_decomp(G_t &G, T_t &T, O_t *O, //BUG: should be optional//,
                      unsigned ub=UINT_MAX /* FIXME: move to backend */,
                      bool ignore_isolated_vertices=false /* FIXME: move to backend */)
{
    // here?
    if(boost::num_vertices(G) == 0){ untested();
        boost::add_vertex(T);
        return 0;
    }else{
    }

    impl::minDegree<G_t, T_t, O_t> MD(G, &T, O, ub, ignore_isolated_vertices);
    MD.do_it();
    MD.tree_decomposition();
    return MD.get_bagsize();
}

template <typename G_t, typename T_t>
typename boost::graph_traits<G_t>::vertices_size_type
  minDegree_decomp(G_t &G, T_t &T)
{

    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return 0;
    }else{
    }

    // hrghs
    typedef typename std::vector<typename treedec_chooser<G_t>::value_type> O_t;

    impl::minDegree<G_t, T_t, O_t> MD(G, &T, (O_t*)NULL, -1u, false);
    MD.do_it();
    MD.tree_decomposition();
    return MD.get_bagsize(); // BUG (still used!)
}


// BUG: duplicate. use impl.
template <typename G_t, typename O_t>
int boost_minDegree_ordering(G_t &G, O_t &O, O_t &iO, unsigned ub = UINT_MAX)
{ untested();
    unsigned n = boost::num_vertices(G);
    unsigned e = boost::num_edges(G);

    O.resize(n);
    unsigned i = 0;
    if(n == 0) { untested();
    }else if(n*(n-1u)==e || e==0){ untested();
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            O[i++] = *vIt;
        }
        if(e==0){ untested();
            return 1;
        }else{ untested();
            return n;
        }
    }else{ untested();
    }

    std::vector<int> inverse_perm(n, 0); // FIXME: use signed_type(vertex_index_t)
    std::vector<int> supernode_sizes(n, 1);
    // typename boost::property_map<G_t, boost::vertex_index_t>::type ...
    auto id = boost::get(boost::vertex_index, G);
    std::vector<int> degree(n, 0);

    /*
     * (Graph& G, 
     *  DegreeMap degree, 
     *  InversePermutationMap inverse_perm, 
     *  PermutationMap perm, 
     *  SuperNodeMap supernode_size, 
     *  int delta, 
     *  VertexIndexMap vertex_index_map)
     */

    int w =
#ifndef HAVE_MINDEGREE_FORK
        0;
    untested();
#endif
    boost::minimum_degree_ordering
             (G,
              boost::make_iterator_property_map(&degree[0], id, degree[0]),
              &iO[0],
              &O[0],
              boost::make_iterator_property_map(&supernode_sizes[0], id, supernode_sizes[0]),
              0,
              id
#ifdef HAVE_MINDEGREE_FORK
              , ub
#endif
              );

#ifdef HAVE_MINDEGREE_FORK
    return w;
#else
    untested();
#endif
}

template <typename G_t, typename O_t>
int boost_minDegree_ordering(G_t &G, O_t &O, unsigned ub=UINT_MAX)
{
    O_t iO(boost::num_vertices(G), 0);
    return boost_minDegree_ordering(G, O, iO, ub);
}


namespace hack_cleanup_later{
template<class ARG>
struct dummy_callback{
    public:
        void operator()(ARG){}
};

}

namespace impl{

template <typename G_t, typename T_t, class CB=hack_cleanup_later::dummy_callback<unsigned>,
        class Graphtype=G_t >
void endless_minDegree_decomp(G_t &G, T_t &T, CB* cb=NULL, Graphtype* =NULL)
{ untested();
    typedef typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> O_t;
    unsigned best = UINT_MAX;
    unsigned n = boost::num_vertices(G);
    unsigned ne = boost::num_edges(G);
    Graphtype Gbak(G);
    T_t U;
    impl::minDegree<G_t, T_t, O_t> MD(G, &U, (O_t*)NULL, best);
    assert(n==boost::num_vertices(G));
    assert(n==boost::num_vertices(Gbak));
    assert(ne==boost::num_edges(G));
    assert(ne==boost::num_edges(Gbak));
    while(true){
        try{
            MD.do_it();
            MD.tree_decomposition(); // FIXME: not here.
            best = get_bagsize(U); // stupid?
            T = U;
            assert(boost::num_vertices(T));
            assert(n==boost::num_vertices(G));
            if(cb){ untested();
                // this is sort of an interruption point.
                (*cb)(best);
            }else{ itested();
            }
        }catch(exception_unsuccessful){
            (*cb)(0);
            // probably not better. retry.
        }
        boost::copy_graph(Gbak, G);
        MD.reset();
        assert(ne==boost::num_edges(G));
        assert(n==boost::num_vertices(G));
    }
}


template <typename G_t>
typename boost::graph_traits<G_t>::vertices_size_type
  boost_minDegree_ordering(G_t &G, std::vector<int> &O)
{
    typedef typename boost::graph_traits<G_t>::edges_size_type edges_size_type;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;

    vertices_size_type n = boost::num_vertices(G);
    edges_size_type e = boost::num_edges(G);

    O.resize(n);

    unsigned i = 0;
    if(n == 0){
        return 0;
    }else if(n*(n-1u) == boost::num_edges(G) || e == 0){
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            O[i++] = *vIt;
        }
        if(e==0){
            return 1;
        }else{
            return n;
        }
    }

    std::vector<int> inverse_perm(n, 0);
    std::vector<int> supernode_sizes(n, 1);
    typename boost::property_map<G_t, boost::vertex_index_t>::type id = boost::get(boost::vertex_index, G);
    std::vector<int> degree(n, 0);

    unsigned w = boost::minimum_degree_ordering
             (G,
              boost::make_iterator_property_map(&degree[0], id, degree[0]),
              &inverse_perm[0],
              &O[0],
              boost::make_iterator_property_map(&supernode_sizes[0], id, supernode_sizes[0]),
              0,
              id);

    return w;
}

} //namespace impl

namespace impl{

template <typename G_t, typename T_t>
void fillIn_decomp(G_t &G, T_t *T, unsigned ub=UINT_MAX)
{
    assert(T);
    typedef typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> O_t;
    fillIn<G_t, T_t, O_t> FI(G, T, (O_t*)NULL, ub);
    FI.do_it();
    FI.tree_decomposition();
}

template <typename G_t, typename T_t>
void fillIn_decomp(G_t &G, T_t &T, unsigned ub=UINT_MAX)
{
    return fillIn_decomp(G, &T, ub);
}

template <typename G_t, class O, class CB=hack_cleanup_later::dummy_callback<unsigned>,
        class Graphtype=G_t >
void endless_fillIn_ordering(G_t &G, O &o, CB* cb=NULL, Graphtype* =NULL)
{ untested();
    typedef typename treedec::graph_traits<G_t>::treedec_type T;
    T t;
    unsigned best = -1u;
    unsigned n = boost::num_vertices(G);
    unsigned ne = boost::num_edges(G);
    Graphtype Gbak(G);
    O U;
    fillIn<G_t, T, O> fill(G, &t, &U, best);
    while(true){ untested();
        try{ untested();
            fill.do_it();
            best = fill.get_bagsize();
            o = U;
            if(cb){ untested();
                // this is sort of an interruption point.
                (*cb)(best);
            }else{ untested();
            }
        }catch(exception_unsuccessful){ untested();
            (*cb)(0);
            // probably not better. retry.
        }
        boost::copy_graph(Gbak, G);
        assert(ne==boost::num_edges(G));
        assert(n==boost::num_vertices(G));
    }
}

} //namespace impl

//Construct a tree decomposition from the elimination ordering obtained by the
//fill-in heuristic.
template <typename G_t, typename T_t>
void fillIn_decomp(G_t &G, T_t &T, unsigned ub=UINT_MAX)
{
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    typename treedec_traits<T_t>::bag_type B;
    BOOST_AUTO(vIt, boost::vertices(G).first);
    BOOST_AUTO(vEnd, boost::vertices(G).second);
    for(; vIt != vEnd; vIt++){
        if(boost::degree(*vIt, G) == 0){
            typename treedec_traits<T_t>::vd_type vd=get_vd(G, *vIt);
            insert(B, vd);
        }
    }

    impl::fillIn_decomp(G, &T, ub);

    auto f=boost::add_vertex(T);
    if(boost::num_vertices(T)!=1){
        boost::add_edge(f, *boost::vertices(T).first, T);
    }

    auto v=B.begin();
    while(true){
        if(v==B.end()){
            break;
        }
        insert(bag(f, T), *v);
        auto b=boost::add_vertex(T);
        boost::add_edge(f, b, T);
        f = b;
        ++v;
    }
}

namespace detail{ //

// Compute an elimination ordering according to minDegree heuristic.
// optionally, treat isolated vertices as deleted.
template<typename G_t>
void minDegree_ordering(G_t& G,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering,
      bool ignore_isolated_vertices=false)
{
    if(ignore_isolated_vertices){ untested();
        // this is not in use... yet?
    }else{
    }
    typedef typename treedec::graph_traits<G_t>::treedec_type T;
    typedef typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> O_t;
    impl::minDegree<G_t, T, O_t> MD(G, elim_ordering, ignore_isolated_vertices);
    MD.do_it();
    MD.elimination_ordering();
}

}

//Compute an elimination ordering according to minDegree heuristic.
template<typename G_t, typename O_t>
void minDegree_ordering(G_t& G, O_t& O)
{
    return detail::minDegree_ordering(G, O, false);
}

namespace detail{ //

//Compute an elimination ordering according to fillIn heuristic (version used
//for postprocessing algorithms).
template<typename G_t, typename O_t>
void fillIn_ordering(G_t &G, O_t &elim_ordering, bool ignore_isolated_vertices=false)
{
    trace3("fillIn_ordering", ignore_isolated_vertices, boost::num_vertices(G), elim_ordering.size());

    typedef typename treedec::graph_traits<G_t>::treedec_type T;
    impl::fillIn<G_t, T, O_t> FI(G, elim_ordering, ignore_isolated_vertices, -1u);
    FI.do_it();
    FI.elimination_ordering();
    assert(elim_ordering.size()==boost::num_vertices(G) || ignore_isolated_vertices);
}

} //detail

//Compute an elimination ordering according to fillIn heuristic.
template<typename G_t>
void fillIn_ordering(G_t& G,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering,
      bool ignore_isolated_vertices=false /* bug, not in frontend! */)
{
    detail::fillIn_ordering(G, elim_ordering, ignore_isolated_vertices);
}

// incomplete: inefficient. see some.h
// (seems to be not in use, remove?)
template <typename G_t, typename O_t>
int get_width_of_elimination_ordering(G_t &G, O_t& elimination_ordering)
{ untested();
    int width = -1;

    for(unsigned int i = 0; i < elimination_ordering.size(); i++){ untested();
        unsigned deg=boost::degree(elimination_ordering[i], G);

        typename graph_traits<G_t>::outedge_set_type xbag;
        treedec::make_clique_and_detach(elimination_ordering[i], G, xbag);
        xbag.clear(); // provide interface with clear included? (not urgent)

        width = (width > (int)deg)? width : (int)deg;
    }

    return width;
}

namespace impl{

template <typename G_t, typename V_t, typename T_t>
void ordering_to_treedec(G_t &G, V_t const& O, T_t &T)
{
    unsigned n = O.size();

    typename std::vector<typename treedec_traits<T_t>::bag_type> bags(n);

    for(unsigned int i = 0; i < O.size(); i++){
        make_clique_and_detach(O[i], G, bags[i]);
    }

    treedec::detail::skeleton_to_treedec(G, T, bags, O, n);
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
{ untested();
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> O_(O.size());
    for(unsigned int i = 0; i < O.size(); i++){ O_[i] = O[i]; }

    ordering_to_treedec(G, O_, T);
}

namespace draft{ //

// TODO: tries hard to add edges at the end. does that make sense?
// TODO: what are the preconditions?!
// TODO: can order be an input range?
// TODO: a test...
template <typename G_t, typename O_t, class T>
void vec_ordering_to_tree(G_t const &G, O_t &O, T& t, O_t* io=NULL,
        boost::adjacency_matrix<boost::directedS> *em=NULL )
{
    size_t num_vert = boost::num_vertices(G);

    if(num_vert == 0){
        boost::add_vertex(t);
        return;
    }

    assert(num_vert = O.size());
    O_t iOlocal;
    typedef boost::adjacency_matrix<boost::directedS> bamd;


    bamd* b;
    if(em){ untested();
        b = em;
    }else{
        // fixme:free
        b = new boost::adjacency_matrix<boost::directedS>(num_vert);
    }
    bamd& bags=*b;

    if(io){
        assert(io->size()==num_vert);
    }else{
        iOlocal.resize(num_vert);
        io=&iOlocal;
    }
    O_t& iO=*io;

    // FIXME: use adjacency matrix
    auto invalid=num_vert;
    std::vector<unsigned> edges(num_vert-1u, invalid);
    assert(edges.size()==num_vert-1);
//    std::vector<std::vector<bool> > bags(num_vert);
//    typedef boost::adjacency_matrix<boost::undirectedS> bamu;

    for(unsigned i = 0; i < num_vert; i++){ itested();
        iO[O[i]] = i;
#ifdef DO_NOT_USE_BAMU
//        bags[i].assign(num_vert, false);
#endif
    }

    for(unsigned i = 0; i < num_vert; i++){ itested();
        auto R=boost::adjacent_vertices(O[i], G);
        for(;R.first!=R.second;++R.first) { itested();
            unsigned n_node = *R.first;
            if((unsigned)iO[n_node] > i){ itested();
                boost::add_edge(i, n_node, bags);
            }
        }
    }

    for(unsigned i = 0; i < num_vert; i++){ itested();
        std::vector<unsigned> N;
        for(unsigned j = 0; j < num_vert; j++){ itested();
            if(boost::edge(i, j, bags).second){ itested();
                N.push_back(j);
                unsigned iO_n_node = iO[j];
                if(iO_n_node < edges[i]){ itested();
                    edges[i] = iO_n_node;
                }
            }
        }

        for(unsigned j = 0; j < N.size(); j++){ itested();
            for(unsigned k = 0; k < N.size(); k++){ itested();
                if(iO[N[k]] > iO[N[j]]){ itested();
                //    bags[iO[N[j]]][N[k]] = true;
                    boost::add_edge(iO[N[j]], N[k], bags);
                    if((unsigned)iO[N[k]] < edges[iO[N[j]]]){ itested();
                        edges[iO[N[j]]] = iO[N[k]];
                    }else{ itested();
                    }
                }
            }
        }
    }

    for(unsigned i = 0; i < num_vert; i++){ itested();
        boost::add_vertex(t);
        auto& b=bag(i,t);
        push(b, O[i]);
        for(unsigned j = 0; j < num_vert; j++){ itested();
            if(boost::edge(i, j, bags).second){ itested();
                push(b, j);
            }
         }
     }

    for(unsigned i = 0; i < num_vert-1u; i++){ itested();
        assert(edges[i]>i || edges[i]==invalid);
        if(edges[i]!=invalid){ itested();
            // normal edge, as computed above.
            boost::add_edge(i, edges[i], t);
        }
        else if(i+1!=num_vert){
            // edge to next component
            boost::add_edge(i, i+1, t);
        } else{ untested();
            // exiting last connected component.
            // ... dont connect
        }
    }

    if(!em){
        delete &bags;
    }
}

} // draft

namespace impl{ //

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

        impl::treedec_to_ordering<G_t, T_t>(T, O);
    }
}

} //namespace impl

template <typename G_t, typename T_t>
void treedec_to_ordering(T_t &T,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &O)
{
    if(boost::num_vertices(T) == 0){ untested();
        return;
    }else if(boost::num_vertices(T) == 1){
        typename boost::graph_traits<T_t>::vertex_descriptor t =
                                                   *(boost::vertices(T).first);
        for(typename treedec_traits<T_t>::bag_type::iterator sIt =
                            bag(t, T).begin(); sIt != bag(t, T).end(); sIt++)
        {
            O.push_back(*sIt);
        }
        return;
    }

    treedec::impl::treedec_to_ordering<G_t, T_t>(T, O);
}

//Make G a filled graph according to the provided elimination_ordering. Stores
//the cliques in C and the additional edges in F.
// BUG: dont use ^make_ if nothing is made.
// BUG: not listed above and not used. remove?
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

// TODO: what does this function do?
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

// TODO: what does this function do?
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

namespace impl{ //

/*
template <typename iO_t, typename M_t, typename G_t>
class elim_predicate{ untested();
public:
    elim_predicate(iO_t &_iO, M_t &_M, G_t &_G) : iO(_iO), M(_M), G(_G){}

    template <typename E_t>
    bool operator()(const E_t &e){ untested();
        unsigned iO_id = iO[boost::source(*e, G)];
        if(M[boost::target(*e, G)] && iO[boost::target(*e, G)] > iO_id){ untested();
            return true;
        }
        return false;
    }
private:
    const M_t &M;
    const iO_t &iO;
    const G_t &G;
};
*/

} //impl

} //namespace treedec

#endif //TD_ELIMINATION_ORDERINGS

// vim:ts=8:sw=4:et
