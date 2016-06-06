// Lukas Larisch, 2014 - 2015
// Felix Salfelder 2016
//
// (c) 2014-2015 Goethe-Universität Frankfurt
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
// A tree decomposition is a graph that has a set of vertex indices as bundled property, e.g.:
//
// struct tree_dec_node
// {
//  std::set<unsigned int> bag;
// };
// typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, tree_dec_node> tree_dec_t;
//
// typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> TD_graph_t;
//
//
// These function is most likely to be interesting for outside use:
//
// void separator_algorithm(G_t &G, T_t &T)
//
//
// For more information, see:
//
//   J. Flum and M. Grohe. 2006. Parameterized Complexity Theory (Texts in Theoretical Computer Science. an EATCS Series).
//   Springer-Verlag New York, Inc., Secaucus, NJ, USA.
//
//   Bruce A. Reed. 1992. Finding approximate separators and computing tree width quickly. In Proceedings of the twenty-fourth
//   annual ACM symposium on Theory of computing (STOC '92). ACM, New York, NY, USA, 221-228.
//
//

#ifndef TD_SEPARATOR_ALGORITHM
#define TD_SEPARATOR_ALGORITHM

#include <set>
#include <vector>

#include "network_flow.hpp"
#include "simple_graph_algos.hpp"
#include "misc.hpp"
#include "graph.hpp"
#include "iter.hpp"

namespace treedec{

//copies all pairs of disjoint subsets of 'X' of size 'min_card' up to
//'max_card' and stores it in 'subsX' and 'subsY'.
//
// outcome
// subsX: all the subsets s of X with min_card <= size(s) <= max_card
// for some i:
// subsX[i] a set
// subsY[i]: all the subsets of X\subsX[i], sized respectively.
// for some j:
// subsY[i][j]: a set equal to subsX[n] for some n (recheck).
//
// could subsY use elements in subsX?
// TODO: use deques?
template <typename X_t>
void disjoint_subsets(X_t const &X, unsigned int min_card, unsigned int max_card,
                      std::vector<typename X_t::value_type> &sub, // internal use?
                      std::vector<X_t> &subsX,
                      std::vector<std::vector<X_t> > &subsY)
{
    assert(!sub.size());
    for(unsigned int i = min_card; i <=max_card; i++){
        subsets(X, X.size(), i, 0, subsX, &sub);
    }
    BOOST_AUTO(I, subsX.begin());
    for(; I!=subsX.end(); ++I){
        X_t difference;

        // compute the complement of subsX[i] wrt X
        std::set_difference(X.begin(), X.end(),
                            I->begin(), I->end(),
                       std::inserter(difference, difference.begin()));

        unsigned int maximum =
             (difference.size() > max_card)? max_card : difference.size();

        std::vector<X_t> subsXY;
        for(unsigned int t = 1; t <= maximum; t++){
            assert(!sub.size());
            subsets(difference, difference.size(), t, 0, sub, subsXY);
        }
        subsY.push_back(subsXY);
    }
    assert(!sub.size());
}

//Collects some vertices of 'V' in 'X' until |X| = size.
template <typename T>
void superset(T &X, T const &V, unsigned int size){
    assert(V.size()>=size); // might not be possible otherwise
    assert(X.size()<=size); // will never terminate...
    typename T::iterator sIt = V.begin();
    while(X.size() != size){
        assert(sIt!=V.end());
        X.insert(*sIt);
        ++sIt;
    }
}

// (proof of concept, generalize later)
// return true, if every element of W is
// an element of S xor an element of X
// aka intersect(W, S) == W \setminus X
template<class W_t, class S_t, class X_t>
inline bool this_intersection_thing(W_t const& W, S_t const& S, X_t const& X)
{
    typename S_t::const_iterator s=S.begin();
    typename X_t::const_iterator x=X.begin();

    for(typename W_t::const_iterator w=W.begin(); w!=W.end(); ++w){

        while(*s<*w && s!=S.end()){
            ++s;
        }
        bool inS=s!=S.end() && *s==*w;

        while(*x<*w && x!=X.end()){
            ++x;
        }
        bool inX=x!=X.end() && *x==*w;

        if(inX==inS){ untested();
            return false;
        }
    }

    return true;
}

//Find a nearly balanced seperator S of W by doing an extended
//deepth-first-search which finds a minimal X-Y-seperator.
//The sets X and Y are all possible disjoint subsets of W of size 1 to 2k.
template <typename G_t, typename W_t, typename S_t>
bool nearly_balanced_seperator(G_t const &G, W_t const &W, S_t &S,
    std::vector<bool> const &disabled, unsigned int k)
{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename std::set<vertex_descriptor> vertex_set;

    typedef std::vector<vertex_descriptor> diff_container_t;
    diff_container_t difference;

#if __cplusplus >= 201103L
    typename subsets_iter<typename W_t::const_iterator>::scratch_type scratch1;
    std::vector<typename std::vector<vertex_descriptor>::iterator > scratch2;
#endif

    // vector?
    vertex_set sX, sY, X_Y;

    for(unsigned s=1; s<=2*k; ++s){
        BOOST_AUTO(P, make_subsets_iter(W.begin(), W.end(), s, s
#if __cplusplus >= 201103L
                , &scratch1
#endif
                ));
        BOOST_AUTO(I, P.first);

        // TODO: don't visit a combination twice.
        for(; I!=W.end(); ++I) {
            difference.resize(0);

        // N = I \setunion neigh(I)
            BOOST_AUTO(N, make_neighbourhood01_iter((*I).first, (*I).second, G, s));

            assert(std::includes( N.first, N.second, N.first, N.second ));
            assert(std::includes( N.first, N.second, (*I).first, (*I).second ));

            // compute W \ (I v neigh(I))
            std::set_difference(W.begin(), W.end(), N.first, N.second,
                           std::inserter(difference, difference.begin()));
            assert(std::includes(W.begin(), W.end(), difference.begin(), difference.end()));

            for(unsigned Js=1; Js<=2*k; ++Js){
                BOOST_AUTO(PP, make_subsets_iter(
                        difference.begin(), difference.end(), Js, Js
#if __cplusplus >= 201103L
                        , &scratch2
#endif
                        ));
                BOOST_AUTO(J, PP.first);
                BOOST_AUTO(e, difference.end());

                for(; J!=e; ++J){
                    S.clear();
                    std::vector<bool> disabled_(disabled);
                    sX.clear();
                    sY.clear();
                    X_Y.clear();

                    // TODO. don't instanciate X_Y. just iterate.
                    // (do we need ordered iterator?)
                    std::set_union((*I).first, (*I).second,
                                   (*J).first, (*J).second,
                                   std::inserter(X_Y, X_Y.begin()));

                    typename vertex_set::const_iterator sIt=X_Y.begin();
                    for(; sIt!=X_Y.end(); ++sIt){
                        unsigned int pos = get_pos(*sIt, G);
                        disabled_[pos] = true;
                    }

                    //Do the extended deepth-first-search on the neighbours of vertices in X and Y
                    get_neighbourhood(G, disabled_, (*I).first, (*I).second, sX);
                    get_neighbourhood(G, disabled_, (*J).first, (*J).second, sY);

                    std::set_difference(W.begin(), W.end(), X_Y.begin(), X_Y.end(),
                                        std::inserter(sX, sX.begin()));
                    std::set_difference(W.begin(), W.end(), X_Y.begin(), X_Y.end(),
                                        std::inserter(sY, sY.begin()));

                    //status1 = nf1::seperate_vertices(G, disabled_, sX, sY, S_, k+1);
                    // network_flow here.
                    if(!treedec::seperate_vertices(G, disabled_, sX, sY, S, k+1)){
                        continue;
                    }

                    //S now is a sX-sY-seperator. Check if S holds the remaining
                    //property of a nearly balanced seperator.

                    // if  ( S \intersect W == W \ X_Y )
                    if(this_intersection_thing(W, S, X_Y)){
                        return true;
                    }else{ untested();
                    }
                }
            } // inner loop
        }
    } // outer loop
    return false;
}

//Glues the 'bag' with 'glueBag' in the current tree decomposition 'T'.
template <typename T_t>
void sep_glue_bag(typename treedec_traits<T_t>::bag_type &b,
                  typename treedec_traits<T_t>::bag_type &glueBag, T_t &T){
    if(boost::num_vertices(T) == 0){
        boost::add_vertex(T);
    }

    typename boost::graph_traits<T_t>::vertex_iterator vertexIt, vertexEnd;
    for(boost::tie(vertexIt, vertexEnd) = boost::vertices(T); vertexIt != vertexEnd; vertexIt++){
        if(bag(*vertexIt, T) == glueBag){
            typename boost::graph_traits<T_t>::vertex_descriptor t_dec_node = boost::add_vertex(T);
            bag(t_dec_node, T) = b;
            boost::add_edge(t_dec_node, *vertexIt, T);
            return;
        }
    }
}

//The main procedure of the seperator algorithm.
// return true if finished (note to self: what does it mean?)
template <typename G_t, typename T_t, class W_t, class P_t, class V_t>
bool sep_decomp(G_t const &G, T_t &T,
        W_t &W, // a vertex set
        P_t const &parent, // a vertex set
        V_t &vertices, // a vertex set
        std::vector<bool> &disabled, unsigned int k)
{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename std::set<vertex_descriptor> vertex_set;
#if 0 // not yet
    typedef typename treedec_traits<T_t>::vd_type vd_type; // vertex identifier. possibly shorter than
                                                           // vertex_descriptor
    typedef typename std::set<vd_type> vd_set;             // just treedec bag_type?
#endif

    //tw(G) > k - one could replace this with a better lower bound (see lower_bounds.hpp).
    if(boost::num_edges(G) > k*boost::num_vertices(G)){
        return false;
    }

    //Passing if V(G) is a subset of W.
    if(std::includes(W.begin(), W.end(), vertices.begin(), vertices.end())){
        return true;
    }

    typename treedec_traits<T_t>::bag_type B1, B2;
    treedec::map_descriptors_to_bags<G_t>(parent, B2);

    //Trivial decomposition
    if(vertices.size() < 4*k + 2){
        treedec::map_descriptors_to_bags<G_t>(vertices, B1);
        treedec::sep_glue_bag(B1, B2, T);
        return true;
    }else if(k==0){ untested();
    }else if(k==1){
    }else{
    }

    //Turn W into a superset of W of size 3k + 1.
    treedec::superset(W, vertices, 3*k + 1);

    vertex_set S;

    //If a nearly balanced seperator S of W' exists, proceed with the graphs
    //induced by the resulting components and the seperator recursively, add a
    //bag containing the union of W and S to the decomposition,
    //connected with the bag, created in the 'parent-call' of the procedure.
    if(nearly_balanced_seperator(G, W, S, disabled, k)){
        std::vector<vertex_set> C;

        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt
                = S.begin(); sIt != S.end(); sIt++) {
            unsigned int pos = get_pos(*sIt, G);
            disabled[pos] = true;
        }

        treedec::get_components_provided_map(G, C, disabled);

        vertex_set union_W_S;
        std::set_union(W.begin(), W.end(), S.begin(), S.end(), std::inserter(union_W_S, union_W_S.begin()));

        //Create a bag (W' v S) and connect it with the bag containing parent.
        treedec::map_descriptors_to_bags<G_t>(union_W_S, B1);
        treedec::sep_glue_bag(B1, B2, T);

        for(unsigned int i = 0; i < C.size(); i++){
            vertex_set union_C_i_S, is_C_i_W, newW;
            std::set_union(C[i].begin(), C[i].end(), S.begin(), S.end(),
                           std::inserter(union_C_i_S, union_C_i_S.begin()));

            std::set_intersection(C[i].begin(), C[i].end(), W.begin(), W.end(),
                                  std::inserter(is_C_i_W, is_C_i_W.begin()));

            std::set_union(is_C_i_W.begin(), is_C_i_W.end(), S.begin(), S.end(),
                           std::inserter(newW, newW.begin()));

            std::vector<bool> disabled_(boost::num_vertices(G), true);
            for(typename vertex_set::iterator sIt
                    = union_C_i_S.begin(); sIt != union_C_i_S.end(); sIt++) {
                unsigned int pos = get_pos(*sIt, G);
                disabled_[pos] = false;
            }

            //Reject if no seperator can be found in one of the ongoing recursive calls.
            if(!sep_decomp(G, T, newW, union_W_S, union_C_i_S, disabled_, k)){
                return false;
            }
        }
        return true;
    }
    return false;
}

//Starts the seperator algorithm, and tries k = 0,1,2,.. until the whole
//graph could be decomposed.
template <typename G_t, typename T_t>
void separator_algorithm(G_t const &G, T_t &T)
{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename std::set<vertex_descriptor> vertex_set;

    unsigned int k = 0;
    bool finished = false;

    typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> vertices;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        vertices.insert(vertices.end(), *vIt);
    }

    while(!finished){
        std::vector<bool> disabled(boost::num_vertices(G), false);
        vertex_set emptySet, parent;
        finished = sep_decomp(G, T, emptySet, parent, vertices, disabled, k);
        k++;

        if(!finished){
            T.clear();
        }
    }
}

} //namespace treedec

#endif //TD_SEPERATOR_ALGORITHM

// vim:ts=8:sw=4:et:
