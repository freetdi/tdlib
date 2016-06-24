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

namespace treedec{

//Create a modified induced subgraph of the bag 'bag(t_desc, T)'.
template <typename G_t, typename T_t, typename M_t, typename S_t>
bool is_improvement_bag(G_t &H,
                        std::vector<bool> &disabled,
                        S_t &X,
                        S_t &Y,
                        typename boost::graph_traits<T_t>::vertex_descriptor t_desc,
                        M_t & vdMap,
                        G_t const &G, T_t &T)
{

    typedef typename boost::graph_traits<G_t> graph_traits;
    typedef typename graph_traits::vertex_iterator vertex_iterator;
    typedef typename graph_traits::adjacency_iterator adjacency_iterator;
    treedec::induced_subgraph(H, G, bag(t_desc, T), vdMap);
    disabled.assign(boost::num_vertices(H), false);

    //Add an additional edge, if a non-edge 'occurs' in a bag of an adjacent
    //vertex t_desc' of t_desc in T.
    // iterate non_edges(H)?
    vertex_iterator vIt1, vIt2, vEnd;
    for(boost::tie(vIt1, vEnd) = boost::vertices(H); vIt1 != vEnd; vIt1++){
        vIt2 = vIt1;
        vIt2++;
        for(; vIt2 != vEnd; vIt2++){
            if(!boost::edge(*vIt1, *vIt2, H).second){
                typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(t_desc, T); nIt != nEnd; nIt++){
                    unsigned int pos1 = get_pos(*vIt1, H);
                    unsigned int pos2 = get_pos(*vIt2, H);
                    typename M_t::value_type vd1 = vdMap[pos1];
                    typename M_t::value_type vd2 = vdMap[pos2];

                    if(bag(*nIt, T).find(vd1) != bag(*nIt, T).end()
                    && bag(*nIt, T).find(vd2) != bag(*nIt, T).end())
                    {
                        boost::add_edge(*vIt1, *vIt2, H);
                        break;
                    }
                }
            }
        }
    }

    //Find a non-edge {x,y} and collect the neighbours of x and y, resulting in the sets X and Y.
    for(boost::tie(vIt1, vEnd) = boost::vertices(H); vIt1 != vEnd; vIt1++){
        vIt2 = vIt1;
        vIt2++;
        for(; vIt2 != vEnd; vIt2++){
            if(!boost::edge(*vIt1, *vIt2, H).second){
                adjacency_iterator nIt, nEnd;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt1, H); nIt != nEnd; nIt++){
                    X.insert(*nIt);
                }

                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt2, H); nIt != nEnd; nIt++){
                    Y.insert(*nIt);
                }

                unsigned int pos1 = get_pos(*vIt1, H);
                unsigned int pos2 = get_pos(*vIt2, H);
                disabled[pos1] = true;
                disabled[pos2] = true;

                goto BREAK_LOOP;
            }
        }
    }

    BREAK_LOOP:

    //Test for completeness.
    typename graph_traits::vertices_size_type nv=boost::num_vertices(H);
    typename graph_traits::edges_size_type ne=boost::num_edges(H);
    if(nv*(nv-1u) == 2*ne){
        H.clear();
        X.clear();
        Y.clear();
        return false;
    }

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
 * Warning: This function is not tested with directed treedecompositions
 *           (and probably computes an invalid treedecomposition. It should
 *            be possible to fix this by re-rooting the resulting treedecomposition).
 */
template <typename G_t, typename T_t>
void MSVS(G_t const &G, T_t &T)
{
    assert(!is_valid_treedecomposition(G, T));
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<T_t>::vertex_descriptor bag_descriptor;
    typedef typename boost::graph_traits<T_t>::vertex_iterator bag_iterator;

    typename std::vector<vertex_descriptor> vdMap_;
    typename std::vector<vertex_descriptor> vdMap;
    std::vector<bool> disabled;
    std::vector<bool> disabled_;
    unsigned width = treedec::get_width(T);
    std::set<vertex_descriptor> S_;
    std::set<vertex_descriptor> S;
    std::set<vertex_descriptor> X, Y;

    while(true){
        INTERRUPTION_POINT;
        // TODO: propagate results to caller
        width = treedec::get_width(T);
//        std::cerr << "MSVS " << width << "\n";

        //Check all maximum sized bags, whether they can be improved or not. Take the first improvable.
        G_t H;
        X.clear();
        Y.clear();
        disabled.resize(0);
        disabled_.resize(0);
        vdMap.resize(0);

        bag_iterator tIt, tEnd;
        bag_descriptor refinement_vertex;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt!=tEnd; ++tIt){
            INTERRUPTION_POINT;
            if(bag(*tIt, T).size() == width+1){
                disabled_.resize(0);
                vdMap_.resize(0);
                if(is_improvement_bag(H, disabled_, X, Y, *tIt, vdMap_, G, T)){
                    refinement_vertex = *tIt;
                    disabled = MOVE(disabled_);
                    vdMap = MOVE(vdMap_);
                    break;
                }
            }
        }

        if(boost::num_vertices(H) == 0){
            //No improvement possible.
            return;
        }

        //Compute a seperating set S.
        S_.clear();
        treedec::seperate_vertices(H, disabled, X, Y, S_);

        //S consists of vertex descriptors of H. Use vd_map to map these to descriptors of G.
        S.clear();
        treedec::map_descriptors(S_, S, H, vdMap);

        //Mark the vertices of the seperator as visited (used for computing connected components).
        std::vector<bool> visited(boost::num_vertices(H), false);
        for(typename std::set<vertex_descriptor>::iterator sIt = S_.begin(); sIt!=S_.end(); ++sIt){
            unsigned int pos = get_pos(*sIt, H);
            visited[pos] = true;
        }

        //Convert the descriptors stored in S to a bag and replace
        //the bag of 'refinement_vertex' with this bag.
        typename treedec_traits<T_t>::bag_type B;
        treedec::map_descriptors_to_bags<G_t>(S, B);
        typename treedec_traits<T_t>::bag_type old_bag = bag(refinement_vertex, T);
        bag(refinement_vertex, T) = MOVE(B);

        //Store the connected components of H[V(H)\S] in 'components'.
        std::vector<std::set<vertex_descriptor> > components;
        treedec::get_components_provided_map(H, components, visited);

        //Store the (neighbours of 'refinement_vertex' in T) in 'oldN'.
        typename boost::graph_traits<T_t>::adjacency_iterator t_nIt, t_nEnd;
        std::vector<bag_descriptor> oldN(boost::degree(refinement_vertex, T));
        unsigned int c = 0;
        for(boost::tie(t_nIt, t_nEnd) = boost::adjacent_vertices(refinement_vertex, T);
              t_nIt != t_nEnd; t_nIt++){
            oldN[c++] = *t_nIt;
        }

        boost::clear_vertex(refinement_vertex, T);

        //'refinement_vertex' gets |connected_components|-many neighbours, that are new vertices in T.
        //The bag of neighbours i will be (connected_components[i] v S).
        std::vector<std::set<vertex_descriptor> > union_S_W_i(components.size());
        std::vector<bag_descriptor> newN(components.size());
        for(unsigned int i = 0; i < components.size(); i++){
            std::set<typename boost::graph_traits<G_t>::vertex_descriptor> component;
            treedec::map_descriptors(components[i], component, H, vdMap);

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
            typename treedec_traits<T_t>::bag_type intersection;
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
