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

#include "trace.hpp"
#include "simple_graph_algos.hpp"
#include "misc.hpp"

#include "fill.hpp"
#include "platform.hpp"

#include "simple_graph_algos.hpp"
#include "misc.hpp"
#include "skeleton.hpp"
#include "treedec.hpp"

#ifndef MINIMUM_DEGREE_ORDERING_HPP
# include "minimum_degree_ordering.hpp"
# define HAVE_MINDEGREE_FORK
#endif

#include "impl/greedy_heuristic.hpp"

#define get_pos(a,b) ( boost::get(boost::vertex_index, b, a) )

namespace treedec{ //

//Construct a tree decomposition from the elimination ordering obtained by the
//minimum-degree heuristic.
//
template <typename G_t, typename T_t, typename O_t>
typename boost::graph_traits<G_t>::vertices_size_type
  minDegree_decomp(G_t &G, T_t &T, O_t *, //TODO: should be optional//,
                      unsigned ub=UINT_MAX /* TODO: move to backend */,
                      bool ignore_isolated_vertices=false /* TODO: move to backend */)
{
    if(boost::num_vertices(G) == 0){ untested();
        boost::add_vertex(T);
        return 0;
    }

    impl::minDegree<G_t> MD(G, ub, ignore_isolated_vertices);
    MD.do_it();
    MD.get_tree_decomposition(T);
    return MD.get_bagsize()-1;
}

template <typename G_t, typename T_t>
typename boost::graph_traits<G_t>::vertices_size_type
  minDegree_decomp(G_t &G, T_t &T)
{

    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return 0;
    }

    impl::minDegree<G_t> MD(G, -1u, false);
    MD.do_it();
    MD.get_tree_decomposition(T);
    return MD.get_bagsize()-1;
}


// TODO: duplicate. use impl.
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
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){ untested();
            O[i++] = *vIt;
        }
        if(e==0){ untested();
            return 1;
        }
        else{ untested();
            return n;
        }
    }

    std::vector<int> inverse_perm(n, 0); //TODO: use signed_type(vertex_index_t)
    std::vector<int> supernode_sizes(n, 1);
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
{ untested();
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
    }
    else if(n*(n-1u) == boost::num_edges(G) || e == 0){
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
void fillIn_decomp(G_t &G, T_t *T, unsigned ub=UINT_MAX, bool ignore_isolated=false)
{
    assert(T);
    fillIn<G_t> FI(G, ub, ignore_isolated);
    FI.do_it();
    FI.get_tree_decomposition(*T);
}

template <typename G_t, typename T_t>
void fillIn_decomp(G_t &G, T_t &T, unsigned ub=UINT_MAX, bool ignore_isolated=false)
{
    return fillIn_decomp(G, &T, ub, ignore_isolated);
}


} //namespace impl

//Construct a tree decomposition from the elimination ordering obtained by the
//fill-in heuristic.
template <typename G_t, typename T_t>
void fillIn_decomp(G_t &G, T_t &T, unsigned ub=UINT_MAX, bool ignore_isolated=false)
{
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }

    impl::fillIn_decomp(G, &T, ub, ignore_isolated);
}

namespace detail{

// Compute an elimination ordering according to minDegree heuristic.
// optionally, treat isolated vertices as deleted.
template<typename G_t>
typename boost::graph_traits<G_t>::vertices_size_type
  minDegree_ordering(G_t& G,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering,
      bool ignore_isolated_vertices=false)
{
    if(ignore_isolated_vertices){ untested();
        //TODO/CHECK: this is not in use... yet?
    }

    impl::minDegree<G_t> MD(G, ignore_isolated_vertices);
    MD.do_it();
    auto o=MD.get_elimination_ordering();
    elim_ordering = o; // HACK

    return MD.get_bagsize()-1;
}

}

//Compute an elimination ordering according to minDegree heuristic.
template<typename G_t, typename O_t>
typename boost::graph_traits<G_t>::vertices_size_type
 minDegree_ordering(G_t& G, O_t& O)
{
    return detail::minDegree_ordering(G, O, false);
}

namespace detail{

//Compute an elimination ordering according to fillIn heuristic (version used
//for postprocessing algorithms).
template<typename G_t, typename O_t>
typename boost::graph_traits<G_t>::vertices_size_type
  fillIn_ordering(G_t &G, O_t &elim_ordering, bool ignore_isolated_vertices=false)
{
    trace3("fillIn_ordering", ignore_isolated_vertices, boost::num_vertices(G), elim_ordering.size());

    impl::fillIn<G_t> FI(G, ignore_isolated_vertices, -1u);
    FI.do_it();
    auto o=FI.get_elimination_ordering();
    elim_ordering = o; // HACK
    assert(elim_ordering.size()==boost::num_vertices(G) || ignore_isolated_vertices);
    return FI.get_bagsize()-1;
}

} //detail

//Compute an elimination ordering according to fillIn heuristic.
template<typename G_t>
typename boost::graph_traits<G_t>::vertices_size_type
 fillIn_ordering(G_t& G,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elim_ordering,
      bool ignore_isolated_vertices=false /* fixme, not in frontend! */)
{
    return detail::fillIn_ordering(G, elim_ordering, ignore_isolated_vertices);
}

// TODO: inefficient
template <typename G_t, typename O_t>
int get_width_of_elimination_ordering(G_t &G, O_t& elimination_ordering)
{
    int width = -1;
    for(unsigned int i = 0; i < elimination_ordering.size(); i++){
        unsigned deg=boost::out_degree(elimination_ordering[i], G);

        typename graph_traits<G_t>::outedge_set_type xbag;

        treedec::make_clique_and_detach(elimination_ordering[i], G, xbag);
        xbag.clear(); // provide interface with clear included? (not urgent)

        width = (width > (int)deg)? width : (int)deg;
    }

    return width;
}


template <typename G_t, typename O_t>
unsigned get_bagsize_of_elimination_ordering(G_t &G, O_t& elimination_ordering)
{
    return get_width_of_elimination_ordering(G, elimination_ordering)+1;
}


namespace impl{

template <typename G_t, typename V_t, typename T_t>
void ordering_to_treedec(G_t &G, V_t const& O, T_t &T)
{
    unsigned n = O.size();
    typedef unsigned vertex_descriptor;

    typename std::vector<
        std::pair<vertex_descriptor,
        typename treedec_traits<T_t>::bag_type>
            > bags(n);

    // stuff center and friends into "skeleton"
    for(unsigned int i = 0; i < O.size(); i++){
        bags[i].first = O[i];
        make_clique_and_detach(O[i], G, bags[i].second);
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

namespace draft{

//TODO
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
        // TODO: free!
        b = new boost::adjacency_matrix<boost::directedS>(num_vert);
    }
    bamd& bags=*b;

    if(io){
        assert(io->size()==num_vert);
    }
    else{
        iOlocal.resize(num_vert);
        io=&iOlocal;
    }
    O_t& iO=*io;

    //TODO: use adjacency matrix
    auto invalid=num_vert;
    std::vector<unsigned> edges(num_vert-1u, invalid);
    assert(edges.size()==num_vert-1);

    for(unsigned i = 0; i < num_vert; i++){
        iO[O[i]] = i;
    }

    for(unsigned i = 0; i < num_vert; i++){
        auto R=boost::adjacent_vertices(O[i], G);
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
        auto& b=boost::get(treedec::bag_t(), t, i);
        push(b, O[i]);
        for(unsigned j = 0; j < num_vert; j++){
            if(boost::edge(i, j, bags).second){
                push(b, j);
            }
         }
     }

    for(unsigned i = 0; i < num_vert-1u; i++){
        assert(edges[i]>i || edges[i]==invalid);
        if(edges[i]!=invalid){
            // normal edge, as computed above.
            boost::add_edge(i, edges[i], t);
        }
        else if(i+1!=num_vert){
            // edge to next component
            boost::add_edge(i, i+1, t);
        }
        else{ untested();
            // exiting last connected component.
            // ... dont connect
        }
    }

    if(!em){
        delete &bags;
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
        if(boost::out_degree(*tIt, T) <= 1 && !bag(*tIt, T).empty()){
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

        if(boost::out_degree(leaf, T) == 1){
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
    }
    else if(boost::num_vertices(T) == 1){
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
template <typename G_t>
void make_filled_graph(G_t &G,
      std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> const &elim_ordering,
      std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > &C,
      std::vector<std::vector<std::pair<
      typename boost::graph_traits<G_t>::vertex_descriptor,
      typename boost::graph_traits<G_t>::vertex_descriptor> > > &F)
{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    C.resize(elim_ordering.size());
    F.resize(elim_ordering.size());

    std::vector<BOOL> visited(boost::num_vertices(G), false);

    for(unsigned int i = 0; i < elim_ordering.size(); i++){
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        std::set<vertex_descriptor> N_i, E_i;
        C[i].insert(elim_ordering[i]);

        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(elim_ordering[i], G); nIt != nEnd; nIt++){
            auto pos=boost::get(boost::vertex_index, G, *nIt);
            if(!visited[pos]){
                C[i].insert(*nIt);
            }
        }

        for(typename std::set<vertex_descriptor>::iterator sIt1 =
            C[i].begin(); sIt1 != C[i].end(); sIt1++)
        {
            typename std::set<vertex_descriptor>::iterator sIt2 = sIt1;
            sIt2++;
            for(; sIt2 != C[i].end(); sIt2++){
                if(!boost::edge(*sIt1, *sIt2, G).second){
                    typename std::pair<vertex_descriptor, vertex_descriptor> edge;
                    edge.first = *sIt1;
                    edge.second = *sIt2;
                    F[i].push_back(edge);
                    boost::add_edge(*sIt1, *sIt2, G);
                }
            }
        }

        auto pos=boost::get(boost::vertex_index, G, elim_ordering[i]);
        visited[pos] = true;
    }
}

template <typename G_t, typename E_t>
void LEX_M_fill_in(G_t &G, E_t &fill_in_edges)
{
    unsigned int nv = boost::num_vertices(G);
    std::vector<BOOL> visited(nv);
    std::vector<float> label(nv);
    std::vector<BOOL> alpha_inv(nv);
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
                        auto edge = std::make_pair(v, *nIt);
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
void LEX_M_minimal_ordering(const G_t &G,
     typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &alpha)
{
    unsigned int nv = boost::num_vertices(G);
    alpha.resize(boost::num_vertices(G));
    std::vector<BOOL> visited(nv);
    std::vector<float> label(nv);
    std::vector<BOOL> alpha_inv(nv);
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
        unsigned max = 0;
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

/*
    unsigned max = 0;
    for(unsigned int j = 0; j < label.size(); j++){ untested();
        max = (label[j] > max)? label[j] : max;
    }

    return (int)max-1; //width of new ordering
*/
}

} //namespace treedec

#undef get_pos

#endif //TD_ELIMINATION_ORDERINGS

// vim:ts=8:sw=4:et
