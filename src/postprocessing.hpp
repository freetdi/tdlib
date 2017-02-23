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
 * Offers functionality to possibly reduce the width of a tree decomposition of a given graph.
 *
 * These functions are most likely to be interesting for outside use:
 *
 * - void MSVS(G_t &G, T_t &T)
 * - void minimalChordal(G_t G, std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &old_elimination_ordering,
 *                              std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &new_elimination_ordering)
 *
 */

#ifndef TD_POSTPROCESSING
#define TD_POSTPROCESSING

#include <boost/graph/adjacency_list.hpp>
#include "elimination_orderings.hpp"
#include "network_flow.hpp"
#include "simple_graph_algos.hpp"
#include "misc.hpp"
#include "graph.hpp"
#include "overlay.hpp"

namespace treedec{

//Check a modified induced subgraph of the bag 'bag(t_desc, T)' for
//possible improvement.
template <typename B_t, typename S_t, typename vd_t>
bool is_improvement_bag(B_t const &H,
                        std::vector<bool> &disabled,
                        S_t &X,
                        S_t &Y,
                        vd_t a,
                        vd_t b)
{

    typedef typename boost::graph_traits<B_t> graph_traits;
    typedef typename graph_traits::adjacency_iterator adjacency_iterator_H;

    disabled.assign(boost::num_vertices(H), false);

    //Test for completeness.
    typename graph_traits::vertices_size_type nv=boost::num_vertices(H);
    typename graph_traits::edges_size_type ne=boost::num_edges(H);
    if(nv*(nv-1u) == 2*ne){
        //There is no remaining non-edge.
        X.clear();
        Y.clear();
        return false;
    }

    assert(a!=b);
    assert(!boost::edge(a, b, H).second);

    //Collect the neighbours of x and y, resulting in the sets X and Y.
    adjacency_iterator_H nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(a, H); nIt!=nEnd; ++nIt){
        X.push_back(*nIt);
    }

    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(b, H); nIt!=nEnd; ++nIt){
        Y.push_back(*nIt);
    }

    unsigned int pos1 = get_pos(a, H);
    unsigned int pos2 = get_pos(b, H);
    disabled[pos1] = true;
    disabled[pos2] = true;

    return true;
}

/* MinimalSeparatingVertexSet(MSVS)-algorithm
 *
 * Tries to find a minimal separator S in a graph H, that
 *    (1) contains the induced subgraph of a maximum-sized bag B(t) of T,
 *    (2) has an edge {x,y}, if {x,y} is a subset of B(t') for some neighbour t' of t in T.
 * If no separator can be found for none of the maximum-sized bags, the algorithm stops. Otherwise,
 * the tree decomposition T is refined according to S.
 *
 */
template <typename G_t, typename T_t>
void MSVS(G_t const &G, T_t &T)
{
    assert(is_valid_treedecomposition(G, T));
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor G_vertex_descriptor;
    typedef typename boost::graph_traits<T_t>::vertex_descriptor bag_descriptor;
    typedef typename boost::graph_traits<T_t>::vertex_iterator bag_iterator;
    typedef typename graph_traits<G_t>::immutable_type immutable_type;
    typedef typename boost::graph_traits<immutable_type>::vertex_descriptor imm_vertex_descriptor;

    std::vector<bool> disabled;
    std::vector<bool> disabled_;
    unsigned width = treedec::get_width(T);
    std::set<vertex_descriptor> S;
    std::vector<imm_vertex_descriptor> X, Y;

    immutable_type H; // malloc/free, where?

    std::set<typename immutable_type::vertex_descriptor> S_;
    std::vector<G_vertex_descriptor> vdMap_, vdMap;
    std::set<typename boost::graph_traits<G_t>::vertex_descriptor> component;

    // FIXME: (terribly) inefficient:
    std::vector<std::set<vertex_descriptor> > union_S_W_i;

    std::vector<bag_descriptor> newN;
    std::vector<bag_descriptor> oldN;
    typename treedec_traits<T_t>::bag_type intersection;

    while(true){
        width = treedec::get_width(T);

        PROPAGATION_POINT;
        INTERRUPTION_POINT;

        //Check all maximum sized bags, whether they can be improved or not. Take the first improvable.
        H.clear();
        X.clear();
        Y.clear();
        disabled.resize(0);
        disabled_.resize(0);
        vdMap.resize(0);

        bag_iterator tIt, tEnd;
        bag_descriptor refinement_vertex;
        immutable_type const* HI=NULL;
        bool status=false;

        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt!=tEnd; ++tIt){
            INTERRUPTION_POINT;
            if(bag(*tIt, T).size() == width+1){
                disabled_.resize(0);
                vdMap_.resize(0);

                is_in_neighbour_bd<vertex_descriptor, T_t> cb(T, *tIt);
                BOOST_AUTO(mybag, bag(*tIt, T));
                HI = &immutable_clone(G, H, mybag.begin(), mybag.end(), mybag.size(), &vdMap_, &cb);
                status = is_improvement_bag
                  <immutable_type, 
                   std::vector<imm_vertex_descriptor>,
                   long unsigned /* H positions... */>
                       (*HI, disabled_, X, Y, cb.a, cb.b);

                if(status){
                    refinement_vertex = *tIt;
                    disabled = MOVE(disabled_);
                    vdMap = MOVE(vdMap_);
                    assert(vdMap.size() == boost::num_vertices(*HI));
                    break;
                }
            }
        }

        if(!status){
            return;
        }
        assert(HI);

#ifndef NDEBUG
        std::vector<bool>::const_iterator x=disabled.begin();
        unsigned num_dis=0;
        for(; x!=disabled.end(); ++x){
            if(*x) ++num_dis;
        }
        assert(num_dis==2);
#endif

        //Compute a seperating set S.
        S_.clear();
        seperate_vertices(H, disabled, 2, X, Y, S_);

        //S consists of vertex descriptors of H. Use vd_map to map these to descriptors of G.
        S.clear();
        map_descriptors(S_, S, *HI, vdMap);

        //Mark the vertices of the seperator as visited (used for computing connected components).
        std::vector<bool> visited(boost::num_vertices(H), false);
        BOOST_AUTO(sIt, S_.begin());
        for(; sIt!=S_.end(); ++sIt){
            unsigned int pos = get_pos(*sIt, *HI);
            visited[pos] = true;
        }

        //Convert the descriptors stored in S to a bag and replace
        //the bag of 'refinement_vertex' with this bag.
        typename treedec_traits<T_t>::bag_type B;
        treedec::map_descriptors_to_bags<G_t>(S, B);
        typename treedec_traits<T_t>::bag_type old_bag = bag(refinement_vertex, T);
        bag(refinement_vertex, T) = MOVE(B);

        //Store the connected components of H[V(H)\S] in 'components'.
        typedef typename boost::graph_traits<immutable_type>::vertex_descriptor HI_vertex_descriptor;
        std::vector<std::set<HI_vertex_descriptor> > components;
        treedec::get_components_provided_map(*HI, components, visited);

        //Store the (neighbours of 'refinement_vertex' in T) in 'oldN'.
        typename boost::graph_traits<T_t>::adjacency_iterator t_nIt, t_nEnd;
        oldN.resize(boost::degree(refinement_vertex, T));
        unsigned int c = 0;
        for(boost::tie(t_nIt, t_nEnd) = boost::adjacent_vertices(refinement_vertex, T);
              t_nIt != t_nEnd; t_nIt++){
            oldN[c++] = *t_nIt;
        }

        boost::clear_vertex(refinement_vertex, T);

        //'refinement_vertex' gets |connected_components|-many neighbours, that are new vertices in T.
        //The bag of neighbours i will be (connected_components[i] v S).
        union_S_W_i.resize(components.size());

        // FIXME: proper container!
        BOOST_AUTO(ii, union_S_W_i.begin());
        for(; ii != union_S_W_i.end(); ++ii ){
            ii->clear();
        }
        newN.resize(components.size());
        for(unsigned int i = 0; i < components.size(); i++){
            component.clear();
            treedec::map_descriptors(components[i], component, *HI, vdMap);

            std::set_union(S.begin(), S.end(),
                           component.begin(), component.end(),
                           std::inserter(union_S_W_i[i], union_S_W_i[i].begin()));

            newN[i] = boost::add_vertex(T);
            typename treedec_traits<T_t>::bag_type uB;
            treedec::map_descriptors_to_bags<G_t>(union_S_W_i[i], uB);
            bag(newN[i], T) = MOVE(uB);

            boost::add_edge(refinement_vertex, newN[i], T);
        }

        //Let intersection_i be the intersection of the old bag of 'refinement_vertex' with
        //the bag of the old neighbour i. Add an edge between i and a new neighbour j of
        //'refinement_vertex', if the bag of j includes intersection_i. This can only be done
        //with exactly one new neighbour of 'refinement_vertex'.
        for(unsigned int i = 0; i <  oldN.size(); i++){
            intersection.clear();
            std::set_intersection(old_bag.begin(), old_bag.end(),
                                  bag(oldN[i], T).begin(),
                                  bag(oldN[i], T).end(),
                                  std::inserter(intersection, intersection.begin()));

            for(unsigned int j = 0; j < newN.size(); j++){
                if(std::includes(bag(newN[j], T).begin(), bag(newN[j], T).end(),
                                 intersection.begin(), intersection.end()))
                {
                    boost::add_edge(newN[j], oldN[i], T);
                    break;
                }
            }
        }
    }
}

template <typename G_t>
bool is_candidate_edge(std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &edge, unsigned int i,
                       std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering, G_t &M)
{
    //Position i in 'elimination_ordering_' will store the 'elimination date' of vertex i
    std::vector<unsigned int> elimination_ordering_(boost::num_vertices(M));
    for(unsigned int t = 0; t < elimination_ordering.size(); t++){
        unsigned int pos = get_pos(elimination_ordering[t], M);
        elimination_ordering_[pos] = t;
    }

    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(edge[0], M); nIt != nEnd; nIt++){
        unsigned int pos = get_pos(*nIt, M);
        if(elimination_ordering_[pos] > i && boost::edge(edge[1], *nIt, M).second
       && !boost::edge(*nIt, elimination_ordering[i], M).second)
        {
            return false;
        }
    }

    return true;
}

template <typename G_t>
inline void delete_edges(G_t &G, std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > &edges){
    for(unsigned int i = 0; i < edges.size(); i++){
        boost::remove_edge(edges[i][0], edges[i][1], G);
    }
}

/* minimalChordal-algorithm
 *
 * Computes possibly redundant fill-in-edges and runs LEX-M to check,
 * if the graph after removal of a fill-in-edge is chordal.
 * Finally, the algorithm computes a new perfect elimination ordering, that
 * possibly causes lower width than 'old_elimination_ordering'.
 */
template <typename G_t>
void minimalChordal(G_t &G,
                typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &old_elimination_ordering,
                typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &new_elimination_ordering)
{
    //Make 'G' a filled-in graph according to 'old_elimination_ordering'. This operation stores
    //all new edges in F.
    std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > C;
    std::vector<std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > > F;
    treedec::make_filled_graph(G, old_elimination_ordering, C, F);

    for(int i = old_elimination_ordering.size()-1; i >= 0; i--){
        //Checks if F[i][j] is an candidate edge. If this is the case, F[i][j] will be stored in
        //'candidate'. The endpoints of F[i][j] will be stored in 'incident'.
        std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > candidate;
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> incident;
        for(unsigned int j = 0; j < F[i].size(); j++){
            if(treedec::is_candidate_edge(F[i][j], i, old_elimination_ordering, G)){
                candidate.push_back(F[i][j]);
                incident.insert(F[i][j][0]);
                incident.insert(F[i][j][1]);
            }
        }
        if(candidate.size() != 0){
            //If there is some candidate edge, create W_i := (incident, G[incident] \ candidate)
            //and run the LEX_M algorithm. The algorithm will possibly return some edges not in W_I that
            //have to be added to W_i to make the graph chordal. 
            G_t W_i;
            typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> vdMap;
            treedec::induced_subgraph_omit_edges(W_i, G, incident, candidate, vdMap);

            std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > keep_fill_;
            treedec::LEX_M_fill_in(W_i, keep_fill_);

            //Translate descriptors of W_i to descriptors of G.
            std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > keep_fill(keep_fill_.size());
            for(unsigned int j = 0; j < keep_fill_.size(); j++){
                unsigned int pos1 = get_pos(keep_fill_[j][0], W_i);
                unsigned int pos2 = get_pos(keep_fill_[j][1], W_i);
                keep_fill[j].push_back(vdMap[pos1]);
                keep_fill[j].push_back(vdMap[pos2]);
            }

            //Delete all candidate edges that can be deleted in G according to LEX_M_fill_in.
            for(unsigned int j = 0; j < candidate.size(); j++){
                for(unsigned int k = 0; k < keep_fill.size(); k++){
                    if((candidate[j][0] == keep_fill[k][0] && candidate[j][1] == keep_fill[k][1])
                     ||(candidate[j][0] == keep_fill[k][1] && candidate[j][1] == keep_fill[k][0]))
                    {
                        candidate.erase(candidate.begin()+j);
                        break;
                    }
                }
            }

            treedec::delete_edges(G, candidate);
        }
    }
    treedec::LEX_M_minimal_ordering(G, new_elimination_ordering);
}

} //namespace treedec

#endif //ifdef TD_POSTPROCESSING

// vim:ts=8:sw=4:et:
