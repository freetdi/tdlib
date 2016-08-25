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
 *
 * Offers functionality to compute a minimal seperator of two vertex sets
 *
 * These functions are most likely to be interesting for outside use:
 *
 * - void seperate_vertices(G_t &G, std::set<unsigned int> &X,
 *                          std::set<unsigned int> &Y, std::set<unsigned int> &S)
 *
 * Computes a seperator S up to size k, aborts and returns false if S would be
 * greater than k:
 *
 * - bool seperate_vertices(G_t &G, std::set<unsigned int> &X,
 *                          std::set<unsigned int> &Y, std::set<unsigned int> &S,
 *                          unsigned int k)
 *
 */

/*
This algorithm is based on Menger's theorem (1927):

  Let G = (V, E) be a graph and A, B subsets of V. Then the minimum number
  of vertices seperating A from B in G is equal to the maximum number of
  disjoint A-B paths in G.

Let P be a family of pairwise disjoint paths in G from X to Y. A
P-alternating walk is a sequence Q = w_1 . . . w_m of vertices of G such
that {w_i, w_(i+1)} is in E for all i in {1, .., m-1} and: (i)  No edge
occurs twice on Q; that is, {w_i, w_(i+1)}  != {w_j , w_(j+1)} for all
distinct i, j in {1, .., m-1}.  (ii) If w_i occurs on a path P = v_1 .. v_l
in P, say w_i = v_j, then w_(i+1) = v_(j−1) or w_(i-1) = v_(j+1).

The algorithm below computes the maximum number of disjoint A-B paths in G
successivly extending a family of disjoint paths P by computing a
P-alternating walk W, if such one exists. If such a walk exists, P can be
extended to P' such that P' containes |P| + 1 disjoint paths.

For some proofs of Menger's theorem, including a contructive one, see

  Reinhard Diestel: Graph Theory, 4th Edition. Graduate texts in mathematics
                    173, Springer 2012, ISBN 978-3-642-14278-9

For a proof of correctness of the algorithm below, see e.g.

  J. Flum and M. Grohe. 2006. Parameterized Complexity Theory (Texts in
                         Theoretical Computer Science. an EATCS Series).
     Springer-Verlag New York, Inc., Secaucus, NJ, USA.
*/

#ifndef TD_NETWORK_FLOW
#define TD_NETWORK_FLOW

#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>

#include "graph.hpp"

#if 0 // huh, defined in graph.hpp
struct Vertex_NF{
    bool visited;
    int predecessor;
};

struct Edge_NF{
    bool path; //true if a path uses the edge
};
#endif

typedef boost::adjacency_list<boost::vecS, boost::vecS,
                          boost::bidirectionalS, treedec::Vertex_NF, treedec::Edge_NF> digraph_t;

namespace treedec{

template<class D>
void check_dis(D dis, size_t num)
{
#ifndef NDEBUG
    typename D::const_iterator i=dis.begin();
    unsigned n=0;
    for(;i!=dis.end();++i) if(*i) ++n;
    assert(n==num);
#endif
}



//Copies the induced subgraph of G formed by disabled into diG. The graph diG
//is a digraph, that is for each edge of the induced subgraph we have an edge
//and a reverse edge. diG has vertex and edge properties, in which the
//iterativly computed edge set of the disjoint paths is stored. A new vertex in
//diG called 'source' is connected to each vertex in X (no reverse edges) and
//each vertex in Y is connected to another new vertex called 'sink' (no
//reverse edges). 'idxMap' is used for the conversion from vertex descriptors
//of diG to vertex descriptors of G, needed to translate a seperator of diG to
//a seperator of G.  Complexity: O(|V| + |E|)
//
// (If G had the properties 'visited' and 'predecessor' on vertices and 'path'
// on edges, the copy step would not be necessary)
//

template <typename G_t>
std::pair<typename boost::graph_traits<typename graph_traits<G_t>::directed_overlay>::vertex_descriptor,
          typename boost::graph_traits<typename graph_traits<G_t>::directed_overlay>::vertex_descriptor>
    make_digraph_with_source_and_sink(G_t const &G, std::vector<bool> const &disabled,
                 unsigned num_dis,
                 typename graph_traits<G_t>::directed_overlay& diG,
                 /* FIXME: use property... */
                 std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &idxMap,
                 typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> const &X,
                 typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> const &Y)
{
    typedef typename graph_traits<G_t>::directed_overlay digraph_t;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;

    check_dis(disabled, num_dis);

    //G may be a graph with ids not in range [0, |V(G)|). The maximum id of a
    //vertex in G is disabled.size().
    std::vector<typename digraph_t::vertex_descriptor> internal_idxMap(disabled.size()+3);
    //needed for linear copy of the edge set of G

    assert(boost::num_vertices(G)>=num_dis);
    unsigned num_dig_verts = boost::num_vertices(G)+2-num_dis;
//    std::cerr << "digaph verts " << num_dig_verts << "\n";

    diG.clear();
    // no resize?! roll out own resize... later.
    // boost::resize(num_dig_verts, diG);

    while(boost::num_vertices(diG)<num_dig_verts){
        boost::add_vertex(diG);
    }
    while(boost::num_vertices(diG)>num_dig_verts){
        // uuh. necessary?
        boost::remove_vertex(*boost::vertices(diG).first, diG);
    }

    unsigned int j = 0;
    BOOST_AUTO(dv, boost::vertices(diG).first);
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int pos = get_pos(*vIt, G);
        if(!disabled[pos]){
            internal_idxMap[pos] = *dv;
            ++dv;
            boost::get(&Vertex_NF::visited, diG, j) = false;
            boost::get(&Vertex_NF::predecessor, diG, j) = j;
            ++j;
            idxMap.push_back(*vIt);
        }
    }

    typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
    for(boost::tie(eIt, eEnd) = boost::edges(G); eIt!=eEnd; ++eIt){
        unsigned int sid = get_pos(boost::source(*eIt, G), G);
        unsigned int tid = get_pos(boost::target(*eIt, G), G);
        if(!disabled[sid] && !disabled[tid]){

            typename digraph_t::edge_descriptor e1 =
                boost::add_edge(internal_idxMap[sid], internal_idxMap[tid], diG).first;
            boost::get(&Edge_NF::path, diG, e1) = false;
            typename digraph_t::edge_descriptor e2 =
                boost::add_edge(internal_idxMap[tid], internal_idxMap[sid], diG).first;
            boost::get(&Edge_NF::path, diG, e2) = false;
        }
    }

    typename digraph_t::vertex_descriptor source = *dv;
    ++dv;
    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt =
            X.begin(); sIt != X.end(); sIt++){
        unsigned int pos = get_pos(*sIt, G);
        typename digraph_t::edge_descriptor e =
            boost::add_edge(source, internal_idxMap[pos], diG).first;
        boost::get(&Edge_NF::path, diG, e) = false;
    }

    typename digraph_t::vertex_descriptor sink = *dv;
    ++dv;
    assert(dv == boost::vertices(diG).second);
    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt =
            Y.begin(); sIt != Y.end(); sIt++){
        unsigned int pos = get_pos(*sIt, G);
        typename digraph_t::edge_descriptor e =
            boost::add_edge(internal_idxMap[pos], sink, diG).first;
        boost::get(&Edge_NF::path, diG, e) = false;
    }

    diG[j].visited = false;
    diG[j].predecessor = j;

    j++;
    diG[j].visited = false;
    diG[j].predecessor = j;

    return std::make_pair(source, sink);
}

namespace impl{
    template <typename G_t>
    class disjoint_ways{
    public:
        // hmm does not work.
        // disjoint_ways(G_t const& g) : dg(g) {}
    public:
        typedef typename graph_traits<G_t>::directed_overlay digraph_t;
        typename graph_traits<G_t>::directed_overlay dg;
        typedef typename boost::graph_traits<digraph_t>::vertex_descriptor divd;
        std::set<divd> dangerous;
        std::vector<std::vector<unsigned int> > P;
    };
} // impl

namespace detail{

//Build the computed disjoint paths by following the edge set of the disjoint
//paths stored in the edge properties of diG, starting in 'source'.
//Complexity: O(|V| + k*|E|), where k is the parameter of the function
//'seperate_vertices' or tw(G), if k is not given.
template<class digraph_t>
static void make_paths(
        digraph_t const &diG, unsigned int source, unsigned int sink,
        std::vector<std::vector<unsigned int> > &P)
{
    typedef typename boost::graph_traits<digraph_t>::out_edge_iterator out_edge_iterator;
    out_edge_iterator nIt1, nEnd1, nIt2, nEnd2;

    unsigned int i = 0;
    if(i<P.size()){
        P[i].clear();
    }else{
    }
    for(boost::tie(nIt1, nEnd1) = boost::out_edges(source, diG); nIt1!=nEnd1; ++nIt1){
        if(boost::get(&Edge_NF::path, diG, *nIt1)){
            typename digraph_t::vertex_descriptor v = boost::target(*nIt1, diG);
            while(true){
                for(boost::tie(nIt2, nEnd2) = boost::out_edges(v, diG); nIt2!=nEnd2; ++nIt2){
                    if(boost::get(&Edge_NF::path, diG, *nIt2)){
                        P[i].push_back(v);
                        v = boost::target(*nIt2, diG);
                        if(v == sink){
                            i++;
                            if(i<P.size()){
                                P[i].clear();
                            }else{
                            }
                            goto NEXT_ITER;
                        }
                        break;
                    }else{
                    }
                }
            }
            // unreachable
        }else{
            // edge source -> v (for some v) but no path....?
        }
        NEXT_ITER:
        ;
    }
}

//An extended depth-first-search: Do a depth-first-search-step if a vertex is
//visited, that is not on a path in P. If a path is entered by visiting a
//vertex, say v, the vertex that has to be visited next is the predecessor of v
//on the path, in which v is contained. Once a P-alternating walk W has been
//computed, the extended family of of disjoint paths P' with respect to W is
//formed by the symetric difference of the edge set of paths in P and the edge
//set of W. The symmetric difference will be immediatly by manipulating the
//edge property of edges, that have been used.  Complexity: O(|V| + |E|) for a
//complete search.
template <typename G_t, typename digraph_t>
static bool t_search_disjoint_ways(
        digraph_t &diG,
        typename boost::graph_traits<digraph_t>::vertex_descriptor v,
        typename boost::graph_traits<digraph_t>::vertex_descriptor sink,
        bool edge_used,
        typename boost::graph_traits<digraph_t>::vertex_descriptor source,
        std::set<typename boost::graph_traits<digraph_t>::vertex_descriptor> &dangerous,
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> const &idxMap,
        G_t const &G)
{
    boost::get(&Vertex_NF::visited, diG, v) = true;
    bool on_a_path = boost::get(&Vertex_NF::predecessor, diG, v) != v;

    //The walk has reached the sink. We can extend the set of disjoint paths by another path.
    if(v == sink){
        return true;
    }

    //Case, that v is on a path in P and the last visited vertex was not the
    //predecessor of v on this path in P.  This vertex could be possibly
    //reached by the predecessor of v on the path at a later time.
    if(on_a_path && !edge_used){
        boost::get(&Vertex_NF::visited, diG, v) = false;
        dangerous.insert(v);

        BOOST_AUTO(&pre, boost::get(&Vertex_NF::predecessor, diG, v));
        BOOST_AUTO(vis, boost::get(&Vertex_NF::visited, diG, pre));
        if(vis){
            return false;
        }else{
            //If a P-alternating walk can be computed by taking this 'reverse edge', P' will not
            //contain {v, w}, where w is the predecessor of v on the path in P that contains v.
            if(t_search_disjoint_ways(diG, pre, sink, true, source, dangerous, idxMap, G)){
                BOOST_AUTO(e, boost::edge(pre, v, diG).first);
                boost::get(&Edge_NF::path, diG, e) = false;
                pre = v;

                return true;
            }
            return false;
        }
    }

    //Do a 'normal' depth-first-search, and ensure that no edge, that is
    //contained on some path in P will be used.
    typename boost::graph_traits<digraph_t>::out_edge_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd)=boost::out_edges(v, diG); nIt!=nEnd; ++nIt){
        BOOST_AUTO(vis, boost::get(&Vertex_NF::visited, diG, boost::target(*nIt, diG)));
        BOOST_AUTO(path,boost::get(&Edge_NF::path, diG, *nIt));
        if(!vis && !path){
            BOOST_AUTO(&pre, boost::get(&Vertex_NF::predecessor, diG, v));
            BOOST_AUTO(T, boost::target(*nIt, diG));
            assert(v==boost::source(*nIt, diG));
            assert(boost::edge(v, T, diG).second);
            bool edge_used_ = pre == (int)T;

            //Recursivly build the walk
            if(t_search_disjoint_ways(diG, T, sink, edge_used_, source, dangerous, idxMap, G)){
                BOOST_AUTO(redge, boost::edge(T, v, diG).first);
                BOOST_AUTO(e, *nIt);

                if(v==source || T==sink){
                    boost::get(&Edge_NF::path, diG, e) = true;
                    boost::get(&Vertex_NF::predecessor, diG, T) = v;
                }else{
                    assert(boost::edge(T, v, diG).second);
                    bool rpath = boost::get(&Edge_NF::path, diG, redge);
                    if(rpath){
                        boost::get(&Edge_NF::path, diG, redge) = false;
                        pre = v;
                    }else{
                        boost::get(&Edge_NF::path, diG, e) = true;
                        boost::get(&Vertex_NF::predecessor, diG, T) = v;
                    }
                }

                return true;
            }
        }else{
            // no path?
        }
    }
    return false;
}

// compute disjoint ways in G \ { v | disabled[v] }
template <typename G_t>
bool disjoint_ways(G_t const &G, std::vector<bool> const &disabled,
        unsigned num_dis,
        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> const &X,
        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> const &Y,
        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &S,
        unsigned int k,
        typename impl::disjoint_ways<G_t>* self)
{
    assert(self);
    typedef typename graph_traits<G_t>::directed_overlay digraph_t;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<digraph_t>::vertex_descriptor divd;
    typedef typename boost::graph_traits<digraph_t>::vertex_iterator divi;

    std::vector<vertex_descriptor> idxMap;

    digraph_t& diG = self->dg;
    std::vector<std::vector<unsigned int> > &P = self->P;
    std::set<divd>& dangerous = self->dangerous;
    divd source, sink;
    boost::tie(source, sink) =
        make_digraph_with_source_and_sink(G, disabled, num_dis, diG, idxMap, X, Y);

    //Main loop of algorithm. min{k+1, |X|+1} iterations are sufficient (one
    //for the unavailing try).
    unsigned int iter = 0;
    for(; iter < X.size()+1; iter++){
        if(S.size()+iter > k){
            return false;
        }

        dangerous.clear();
        //start extended DFS at source
        if(!t_search_disjoint_ways(diG, source, sink, false, source, dangerous, idxMap, G)){
            for(typename std::set<divd>::iterator sIt =
                  dangerous.begin(); sIt != dangerous.end(); sIt++){
                boost::get(&Vertex_NF::visited, diG, *sIt) = true;
            }

            break;
        }

        //undo visited
        divi vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(diG); vIt != vEnd; vIt++){
            boost::get(&Vertex_NF::visited, diG, *vIt) = false;
        }
    }

    P.resize(iter);

    make_paths(diG, source, sink, P);

    //Compute a separator by taking the last vertex on each path, that could be
    //reached by a P-alternating walk. If no such separator exists, the first
    //vertex on a path is taken.
    BOOST_AUTO(Pi, P.begin());
    for(; Pi!=P.end(); ++Pi){
        BOOST_AUTO(Pij, Pi->rbegin());
        for(; Pij!=Pi->rend(); ++Pij){
            if(boost::get(&Vertex_NF::visited, diG, *Pij)){
                S.insert(idxMap[*Pij]);
                goto there;
            }else{
            }
        }
        if(Pij==Pi->rend()){
            S.insert(idxMap[*Pi->begin()]);
        }else{
        }
there:
        ;
    }

    return true;
}

} //namespace detail

//This version immediatly aborts after at most k+1 iterations. The return value
//indicates, whether a seperator of size at most k
//exists or not.
template <typename G_t, typename S_t>
bool seperate_vertices(
        G_t const &G, std::vector<bool> &disabled, unsigned &num_dis,
        S_t const &X,
        S_t const &Y,
        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &S,
        unsigned int k,
        impl::disjoint_ways<G_t>* dw=NULL)
{
    assert(dw);
    //Common neighbours must be contained in a seperator.
    std::set_intersection(X.begin(), X.end(), Y.begin(), Y.end(), std::inserter(S, S.begin()));

    typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> X_, Y_;
    std::set_difference(X.begin(), X.end(), S.begin(), S.end(), std::inserter(X_, X_.begin()));
    std::set_difference(Y.begin(), Y.end(), S.begin(), S.end(), std::inserter(Y_, Y_.begin()));

    if(S.size() > k){
        return false;
    }

    if(X_.size() == 0 || Y_.size() == 0){
        return true;
    }

    //disables/deletes the vertices in the intersection of X and Y.
    for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt =
          S.begin(); sIt != S.end(); sIt++){
        unsigned int pos = get_pos(*sIt, G);
        assert(!disabled[pos]);
        ++num_dis;
        disabled[pos] = true;
    }

    return detail::disjoint_ways(G, disabled, num_dis, X_, Y_, S, k, dw);
}

//Version that computes a X-Y-seperator S without aborting after k iterations (S really will be a seperator).
template <typename G_t, typename S_t, typename Sep_t>
void seperate_vertices(G_t &G, std::vector<bool> &disabled, unsigned num_dis,
        S_t const &X,
        S_t const &Y,
        Sep_t &S,
        impl::disjoint_ways<G_t>* dw=NULL)
{
    typedef impl::disjoint_ways<G_t> dw_t;
    bool own=false;
    if(!dw){
        own=true;
        dw=new dw_t;
    }
    assert(dw);
    seperate_vertices(G, disabled, num_dis, X, Y, S, UINT_MAX, dw);
    if(own){
        delete dw;
    }
}

} //namespace treedec

#endif //TD_NETWORK_FLOW

// vim:ts=8:sw=4:et
