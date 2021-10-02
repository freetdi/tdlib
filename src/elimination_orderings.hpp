// Lukas Larisch, 2014 - 2017
// Felix Salfelder, 2016 - 2017, 2021
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

#ifndef TREEDEC_ELIMINATION_ORDERINGS_HPP
#define TREEDEC_ELIMINATION_ORDERINGS_HPP

#include <cmath>
#include <climits>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adjacency_matrix.hpp>
// #include <boost/graph/graph_utility.hpp> // print

#include <boost/graph/copy.hpp>

#include "trace.hpp"
#include "simple_graph_algos.hpp"
#include "misc.hpp"

#include "fill.hpp"
#include "platform.hpp"

#include "misc.hpp"
#include "skeleton.hpp"
#include "treedec.hpp"

#ifndef MINIMUM_DEGREE_ORDERING_HPP
# include "minimum_degree_ordering.hpp"
# define HAVE_MINDEGREE_FORK
#endif

#include "impl/greedy_heuristic.hpp"

#ifdef HAVE_GALA_GRAPH_H
#include <gala/sfinae.h>
#endif

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

namespace draft{

template<class G, class O>
class io_smaller_than{
    typedef typename boost::property_map< G, boost::vertex_index_t >::const_type::value_type vertex_index_type;
public:
    explicit io_smaller_than(size_t k, O const& o, G const& g) : _k(k), _o(o), _g(g) {
    }
    template<class E>
    bool operator()(E const& e) const{
        auto t = boost::target(e, _g);
        auto idm = boost::get(boost::vertex_index, _g);
        if(idm[_o[t]] < _k){
            return true;
        }else{
            return false;
        }
    }
private:
    vertex_index_type _k;
    O const& _o;
    G const& _g;
};

template<class N>
class mapped_order{
public:
    mapped_order(N const& n):_n(n){}

    template<class V, class W>
    bool operator()(V const& a, W const& b) const{
        return _n[a] < _n[b];
    }
private:
    N const& _n;
};

template<class I, class N, class G, class O, class M>
void cleanup_bmdo(I j, N const& numbering, G& g, O const&, M const& my_numbering_order)
{
    io_smaller_than<G, O> P(numbering[j], numbering, g);
    vertex_descriptor_G b(j);
    boost::remove_out_edge_if(b, P, g);

    std::sort(g->vertices()[j].begin(), g->vertices()[j].end(), my_numbering_order);
}

template <typename G, typename O_t, class T>
void inplace_bmdo_tree(G &g, O_t const& O, T& t, size_t bagsize, O_t const& io_)
{
    auto const* io = &io_;
    typedef typename boost::property_map<G, boost::vertex_index_t>::type::value_type vertex_index_type;
    size_t num_vert = boost::num_vertices(g);

    if(num_vert == 0){ untested();
        boost::add_vertex(t);
    }else{

        assert(num_vert == O.size());
        O_t iOlocal;

        if(io){
            trace2("DBG", io->size(), num_vert);
            assert(io->size()==num_vert);
            for(vertex_index_type i = 0; i < num_vert; i++){
                // bug? does O map to vertex_descriptors?
                assert(vertex_index_type((*io)[O[i]]) == i);
            }
        }else{ untested();
            iOlocal.resize(num_vert);
            io=&iOlocal;
            for(unsigned i = 0; i < num_vert; i++){
                iOlocal[O[i]] = i;
            }
        }
        //    O_t const& iO=*io;
        O_t const& numbering=*io;

        auto invalid=num_vert;
        std::vector<unsigned> edges(num_vert-1u, invalid);
        assert(edges.size()==num_vert-1);

        mapped_order<O_t> my_numbering_order(numbering);

        for(unsigned j = 0; j < num_vert-bagsize; j++){
            cleanup_bmdo(O[j], numbering, g, O, my_numbering_order);
        }

        std::vector<vertex_descriptor_G> buf;
        auto nodes_left = O.size();

        for(auto oi_ = O.begin(); ; ++oi_){
            auto oi = *oi_;
            --nodes_left;
            auto i = numbering[oi];
            auto R = boost::adjacent_vertices(oi, g);

            for(;R.first!=R.second;++R.first) {
                auto j = *R.first;
                unsigned iO_n_node = numbering[j];
                if(iO_n_node < edges[i]){
                    edges[i] = iO_n_node;
                }else{
                }
            }

            if(bagsize == nodes_left + 1){
                // the rest is one big bag, no matter what.
                auto& last_adj = g->vertices()[oi];
                last_adj.resize(bagsize-1);
                size_t f = 0;
                while(++oi_ != O.end()){
                    last_adj[f++] = *oi_;
                }
                assert(last_adj.size() + 1 == bagsize);
                break;
            }else{
            }

            auto const& NN = g->vertices()[oi];
            // NN is now "bag O[i]"
            // R = boost::adjacent_vertices(oi, g);

            // rewire bag at O[i]. will become bag i.
            for( auto j_ = NN.begin(); j_ != NN.end(); ++j_ ) {
                auto k_ = j_;
                ++k_;
                auto j = *j_;

                buf.resize(0);

                auto Aj = boost::adjacent_vertices(j, g);

                // could try canonical order... but then NN is wrong.
                //
                if(numbering[j] + bagsize < num_vert){
                    std::set_union(k_, NN.end(), Aj.first, Aj.second, std::back_inserter(buf), my_numbering_order);
                    std::swap(g->vertices()[j], buf);
                }else{
                }

                for( ; k_!=NN.end() ; ++k_ ){
                    auto k = *k_;
                    assert(numbering[j] < numbering[k]);

                    if((unsigned)numbering[k] < edges[numbering[j]]){
                        edges[numbering[j]] = numbering[k];
                    }else{
                    }
                }
            }
        }

        trace2("loop done", nodes_left, num_vert);
        assert(! boost::num_vertices(t));

        for(unsigned i = 0; i < num_vert-nodes_left; ++i){
            boost::add_vertex(t);
            assert(i+1 == boost::num_vertices(t));
            auto& b = boost::get(treedec::bag_t(), t, i);

            auto& NN = g->vertices()[O[i]];
            assign(b, std::move(NN));
            push(b, O[i]);
        }

        assert(boost::num_vertices(t) == num_vert-nodes_left);

        // invert edge direction?
        for(unsigned i = 0; i < num_vert-nodes_left-1u; ++i){
            assert(edges[i]>i || edges[i]==invalid);
            if(edges[i] >= num_vert-nodes_left){
                edges[i] = num_vert-nodes_left-1;
                boost::add_edge(i, edges[i], t);
            }else if(edges[i]!=invalid){
                // normal edge, as computed above.
                boost::add_edge(i, edges[i], t);
            }else if(i+1!=num_vert){ untested();
                // edge to next component
                boost::add_edge(i, i+1, t);
            }else{ untested();
                incomplete(); // ?
                trace2("edging exit", i, edges[i]);
                // exiting last connected component.
                // ... dont connect
            }
        }

        assert(boost::num_vertices(t) == num_vert-nodes_left);
        assert(boost::num_edges(t) +1 == boost::num_vertices(t));

        for(unsigned i = 0; i < num_vert-nodes_left; i++){
            auto& b=boost::get(treedec::bag_t(), t, i);
            assert(b.size());
            treedec::sort(b); // bug in check_tree_decomp, need sorted bags...
        }
    }

} // inplace_bmdo

template <typename G_t, typename O_t, class T, class N>
void vec_ordering_to_tree(G_t const &G, O_t const& O, T& t, N* io=NULL,
        boost::adjacency_matrix<boost::directedS> *em=NULL )
{
    size_t num_vert = boost::num_vertices(G);

    if(num_vert == 0){
        boost::add_vertex(t);
        return;
    }

    assert(num_vert == O.size());
    O_t iOlocal;
    typedef boost::adjacency_matrix<boost::directedS> bamd;


    bamd* b{nullptr};
    if(em){ untested();
        b = em;
    }else{
        // TODO: free!
        b = new boost::adjacency_matrix<boost::directedS>(num_vert);
    }
    bamd& bags=*b;

    if(io){
        assert(io->size()==num_vert);
    }else{
        iOlocal.resize(num_vert);
        io=&iOlocal;
        for(unsigned i = 0; i < num_vert; i++){
            iOlocal[O[i]] = i;
        }
    }
    N& iO=*io;

    //TODO: use adjacency matrix
    auto invalid=num_vert;
    std::vector<unsigned> edges(num_vert-1u, invalid);
    assert(edges.size()==num_vert-1);


    for(unsigned i = 0; i < num_vert; i++){
        auto R = boost::adjacent_vertices(O[i], G);
        for(; R.first!=R.second; ++R.first) {
            unsigned n_node = *R.first;
            if((unsigned)iO[n_node] > i){
                boost::add_edge(i, n_node, bags);
            }else{
            }
        }
    }

    for(unsigned i = 0; i < num_vert; i++){
        std::vector<unsigned> neigh;
        for(unsigned j = 0; j < num_vert; j++){
            if(boost::edge(i, j, bags).second){
                neigh.push_back(j);
                unsigned iO_n_node = iO[j];
                if(iO_n_node < edges[i]){
                    edges[i] = iO_n_node;
                }
            }
        }

        for(unsigned j = 0; j < neigh.size(); j++){
            for(unsigned k = 0; k < neigh.size(); k++){
                if(iO[neigh[k]] > iO[neigh[j]]){
                    boost::add_edge(iO[neigh[j]], neigh[k], bags);
                    if((unsigned)iO[neigh[k]] < edges[iO[neigh[j]]]){
                        edges[iO[neigh[j]]] = iO[neigh[k]];
                    }else{
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
            }else{
            }
         }
     }

    for(unsigned i = 0; i < num_vert-1u; i++){
        assert(edges[i]>i || edges[i]==invalid);
        if(edges[i]!=invalid){
            // normal edge, as computed above.
            boost::add_edge(i, edges[i], t);
        }else if(i+1!=num_vert){
            // edge to next component
            boost::add_edge(i, i+1, t);
        }else{ untested();
            // exiting last connected component.
            // ... dont connect
        }
    }

    if(!em){
        delete b;
    }else{
    }
}

} // draft

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
              &inverse_perm[0], // numbering. node n is at position ip[n]
              &O[0],   // "ordering" as in eliminate O[0] then O[1] ...
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
    treedec::obsolete::fillIn<G_t> FI(G, ub, ignore_isolated);
    FI.do_it();
    FI.get_tree_decomposition(*T);
}

template <typename G_t, typename T_t>
void fillIn_decomp(G_t &G, T_t &T, unsigned ub=UINT_MAX, bool ignore_isolated=false)
{
    return fillIn_decomp(G, &T, ub, ignore_isolated);
}


} //namespace impl

namespace impl{

template<class G, class T, class X=void>
struct bmdo_{
    template<class D, class O, class N>
    static void gtd(D const& g, O o, T& t, size_t, N const& numbering){
        // ordering_to_treedec(_g, *_o, t);
        treedec::draft::vec_ordering_to_tree(g, o, t, &numbering);
    }
};

#ifdef HAVE_GALA_GRAPH_H // TODO: define graph trait
template<class G, class T>
struct bmdo_<G, T, typename gala::sfinae::is_vector<typename G::vertex_container_type>::type>
{
    template<class D, class O, class N>
    static void gtd(D& g, O o, T& t, size_t bs, N const& numbering){
        treedec::draft::inplace_bmdo_tree(g, o, t, bs, numbering);
    }
};
#endif


// hack?
//using boost::target;
//using boost::source;

// returns BAGSIZE
template <typename G_t>
class bmdo{
public:
    typedef typename boost::graph_traits<G_t>::edges_size_type edges_size_type;
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
    typedef treedec::draft::directed_view<G_t> D_t;
public:
    bmdo(G_t &G, std::vector<int> &O)
      : _g(G),
        _o(&O),
        _ub(vertices_size_type(-1u)) { untested();
    }
    bmdo(G_t &G)
      : _g(G),
        _o(new std::vector<int>()),
        _own_o(true),
        _ub(vertices_size_type(-1)) {
    }
    ~bmdo(){
        if(_own_o){
            delete _o;
        }else{
        }
    }

public:
    vertices_size_type bagsize() const{
        return _bs;
    }
    unsigned lower_bound_bagsize() const{
        incomplete();
        return 0;
    }
    template<class T>
    void get_tree_decomposition(T& t) const{
        incomplete();
    }
    template<class T>
    void get_tree_decomposition(T& t){
        bmdo_<G_t, T>::gtd(_g, o(), t, _bs, _inverse_perm);
    }
    template<class O>
    void get_elimination_ordering(O& v) const{ untested();
        if(_o){ untested();
            v = O(_o->begin(), _o->end());
        }else{ untested();
            incomplete();
        }
    }
    void do_it();
private:
    std::vector<int>& o(){assert(_o); return *_o;}
    std::vector<int> const& o() const{assert(_o); return *_o;}
private:
    D_t _g;
    std::vector<int> _inverse_perm; //TODO: use signed_type(vertex_index_t)
    std::vector<int>* _o{nullptr};
    bool _own_o{false};
    vertices_size_type _bs;
    vertices_size_type _ub;
}; // bmdo

template<class G_t>
void bmdo<G_t>::do_it()
{
    vertices_size_type n=boost::num_vertices(_g);
    edges_size_type e=boost::num_edges(_g);

    o().resize(n);
    _inverse_perm.resize(n);

    // check: is this still necessary?!
    unsigned i = 0;
    if(n == 0){ untested();
        _bs = 0;
    }else if(n*(n-1u) == boost::num_edges(_g) || e == 0){
        auto p=boost::vertices(_g);
        for(; p.first!=p.second; ++p.first){
            _inverse_perm[*p.first] = i;
            o()[i++] = *p.first;
        }
        if(e==0){
            _bs = 1;
        }else{
            _bs = n;
        }
    }else{

        std::vector<int> supernode_sizes(n, 1);
        auto id = boost::get(boost::vertex_index, _g);
        std::vector<int> degree(n, 0);

        _bs =
#ifndef HAVE_MINDEGREE_FORK
        0;
    untested();
#endif
        boost::minimum_degree_ordering
                 (_g,
                  boost::make_iterator_property_map(&degree[0], id, degree[0]),
                  &_inverse_perm[0],
                  &o()[0],
                  boost::make_iterator_property_map(&supernode_sizes[0], id, supernode_sizes[0]),
                  0,
                  id
#ifdef HAVE_MINDEGREE_FORK
              , _ub
#endif
                  );
    }
} // bmdo::do_it

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

    treedec::impl::minDegree<G_t> MD(G, ignore_isolated_vertices);
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

    obsolete::fillIn<G_t> FI(G, ignore_isolated_vertices, -1u);
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
        trace2("picked", O[i], bags[i].second.size());
    }

    treedec::detail::skeleton_to_treedec(G, T, bags, O, n);
}

} //namespace impl

template <typename G_t, typename T_t, class O_>
void ordering_to_treedec(G_t &G, O_ const& O, T_t &T)
{
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
    }else{
        treedec::impl::ordering_to_treedec(G, O, T);
    }
}

template <typename G_t, typename T_t>
void ordering_to_treedec(G_t &G, std::vector<int> &O, T_t &T)
{ untested();
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> O_(O.size());
    for(unsigned int i = 0; i < O.size(); i++){ untested();
        O_[i] = O[i];
    }

    ordering_to_treedec(G, O_, T);
}

namespace impl{

template <class T_t, class O_t>
void treedec_to_ordering(T_t &T, O_t& O)
{
    bool leaf_found = false;

    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    typename boost::graph_traits<T_t>::vertex_descriptor leaf, parent;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        auto const& b=boost::get(bag_t(), T, *tIt);
        if(boost::out_degree(*tIt, T) <= 1 && !b.empty()){
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
            auto const& pb=boost::get(bag_t(), T, parent);
            auto const& lb=boost::get(bag_t(), T, leaf);
            if(!std::includes(pb.begin(), pb.end(), lb.begin(), lb.end())) {
                std::set_difference(lb.begin(), lb.end(), pb.begin(), pb.end(),
                                    std::inserter(difference, difference.begin()));
            }
            boost::clear_vertex(leaf, T);
        }
        else{
            difference = MOVE(boost::get(bag_t(), T, leaf));
        }

        for(typename treedec_traits<T_t>::bag_type::iterator sIt = difference.begin();
            sIt != difference.end(); sIt++)
        {
            O.push_back(*sIt);
        }

        auto& b=boost::get(bag_t(), T, leaf);
        b.clear();

        impl::treedec_to_ordering<T_t, O_t>(T, O);
    }else{
    }
}

} //namespace impl

template <typename T_t, typename O_t>
void treedec_to_ordering(T_t &T, O_t& O)
{
    if(boost::num_vertices(T) == 0){ untested();
        return;
    }
    else if(boost::num_vertices(T) == 1){
        typename boost::graph_traits<T_t>::vertex_descriptor t =
                                                   *(boost::vertices(T).first);
        auto& b=boost::get(bag_t(), T, t);
        for(auto sIt=b.begin(); sIt != b.end(); ++sIt) {
            O.push_back(*sIt);
        }
        return;
    }

    treedec::impl::treedec_to_ordering<T_t, O_t>(T, O);
}

//Make G a filled graph according to the provided elimination_ordering. Stores
//the cliques in C and the additional edges in F.
template <typename G_t, class O_>
void make_filled_graph(G_t &G, O_ const& elim_ordering,
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

template <typename G_t, class O_>
void LEX_M_minimal_ordering(const G_t &G, O_& alpha)
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

// start heuristics.hpp?
namespace he {
    // yes, no camel case, as all the others.
    template<class x, template<class G_, class ...> class C=treedec::algo::default_config>
    using fill_in=treedec::impl::fillIn<x, C>;
}

template <typename G_t, typename O_t>
int boost_minDegree_ordering(G_t &G, O_t &O, unsigned ub=UINT_MAX)
{ untested();
    impl::bmdo<G_t> b(G);
    b.set_ub(ub);
    b.do_it();
    b.get_elimination_ordering(O);
    return b.get_bagsize()-1;
}


} //namespace treedec

#undef get_pos

#endif //TREEDEC_ELIMINATION_ORDERINGS_HPP

// vim:ts=8:sw=4:et
