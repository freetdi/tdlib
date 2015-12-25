// Lukas Larisch, 2014 - 2015
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
// Offers functionality to possibly reduce the width of tree decompositions of graphs
//
// A tree decomposition is a graph that has a set of vertex indices as bundled property, e.g.:
//
// struct tree_dec_node
// {
//  std::set<unsigned int> bag;
// };
// typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, tree_dec_node> tree_dec_t;
//
// Vertices of the input graph have to provide the attribute 'id', e.g.:
//
// struct Vertex
// {
//  unsigned int id;
// };
// typedef boost::adjacency_list<boost::setS, boost::listS, boost::undirectedS, Vertex> TD_graph_t;
//
//
// These functions are most likely to be interesting for outside use:
//
// void MSVS(G_t &G, T_t &T)
// void minimalChordal(G_t G, std::vector<unsigned int> &old_elimination_ordering, std::vector<unsigned int> &new_elimination_ordering)
//

#ifndef TD_POSTPROCESSING
#define TD_POSTPROCESSING

#include <boost/graph/adjacency_list.hpp>
#include "TD_elimination_orderings.hpp"
#include "TD_NetworkFlow.hpp"
#include "TD_simple_graph_algos.hpp"
#include "TD_misc.hpp"

namespace treedec{

//Creates a modified induced subgraph of the bag 'T[t_desc].bag'.
template <typename G_t, typename T_t>
bool is_improvement_bag(G_t &H, std::vector<bool> &disabled, std::set<unsigned int> &X, std::set<unsigned int> &Y, typename boost::graph_traits<T_t>::vertex_descriptor t_desc, G_t &G, T_t &T){
    induced_subgraph(H, G, T[t_desc].bag);

    //Add an additional edge, if a non-edge 'occures' in a bag of an adjacent vertex t_desc' of t_desc in T.
    typename boost::graph_traits<G_t>::vertex_iterator vIt1, vIt2, vEnd;
    for(boost::tie(vIt1, vEnd) = boost::vertices(H); vIt1 != vEnd; vIt1++){
        vIt2 = vIt1;
        vIt2++;
        for(; vIt2 != vEnd; vIt2++){
            if(!boost::edge(*vIt1, *vIt2, H).second){
                typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(t_desc, T); nIt != nEnd; nIt++){
                    if(T[*nIt].bag.find(H[*vIt1].id) != T[*nIt].bag.end() && T[*nIt].bag.find(H[*vIt2].id) != T[*nIt].bag.end()){
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
                typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt1, H); nIt != nEnd; nIt++)
                    X.insert(H[*nIt].id);

                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt2, H); nIt != nEnd; nIt++)
                    Y.insert(H[*nIt].id);

                disabled[H[*vIt1].id] = true;
                disabled[H[*vIt2].id] = true;

                goto BREAK_LOOP;
            }
        }
    }

    BREAK_LOOP:

    //Test for completeness.
    if(boost::num_vertices(H)*(boost::num_vertices(H)-1) == 2*boost::num_edges(H)){
        H.clear();
        X.clear();
        Y.clear();
        return false;
    }

    return true;
}

/* MinimalSeperatingVertexSet(MSVS)-algorithm
 *
 * Tries to find a minimal seperator S in the graph H, that
 *    (1) containes the induced subgraph of a maximum-sized bag B(t) of T,
 *    (2) has an edge {x,y}, if {x,y} is a subset of B(t') for some neighbour t' of t in T.
 * If no seperator can be found for none of the maximum-sized bags, the algorithm stops. Otherwise,
 * the tree decomposition T is refined according to S.
 */
template <typename G_t, typename T_t>
void MSVS(G_t &G, T_t &T){
    while(true){
        unsigned int width = treedec::get_width(T);

        //Check all maximum sized bags, whether they can be improved or not. Take the first improvable.
        G_t H;
        std::set<unsigned int> X, Y;
        std::vector<bool> disabled(boost::num_vertices(G), false);

        typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
        typename boost::graph_traits<T_t>::vertex_descriptor refinement_bag;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            if(T[*tIt].bag.size() == width+1){
                std::vector<bool> disabled_(disabled);
                if(is_improvement_bag(H, disabled_, X, Y, *tIt, G, T)){
                    refinement_bag = *tIt;
                    disabled = disabled_;
                    break;
                }
            }
        }

        //No improvement possible.
        if(boost::num_vertices(H) == 0){
            return;
        }

        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> X_, Y_, S_;
        treedec::id_bag_to_descriptor_bag(H, X_, X);
        treedec::id_bag_to_descriptor_bag(H, Y_, Y);

        //Compute a seperating set S.
        std::set<unsigned int> S;
        treedec::seperate_vertices(H, disabled, X_, Y_, S_);

        treedec::descriptor_bag_to_id_bag(H, S, S_);

        //Do the refinement.
        std::vector<bool> visited(boost::num_vertices(G), true);
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd)= boost::vertices(H); vIt != vEnd; vIt++)
            visited[H[*vIt].id] = false;

        for(std::set<unsigned int>::iterator sIt = S.begin(); sIt != S.end(); sIt++)
            visited[*sIt] = true;

        std::vector<std::set<unsigned int> > components;
        get_components_provided_map(H, components, visited);

        typename boost::graph_traits<T_t>::adjacency_iterator t_nIt, t_nEnd;
        std::vector<typename boost::graph_traits<T_t>::vertex_descriptor> oldN;
        std::vector<typename boost::graph_traits<T_t>::vertex_descriptor> newN(components.size());

        for(boost::tie(t_nIt, t_nEnd) = boost::adjacent_vertices(refinement_bag, T); t_nIt != t_nEnd; t_nIt++)
            oldN.push_back(*t_nIt);

        boost::clear_vertex(refinement_bag, T);

        std::set<unsigned int> old_bag = T[refinement_bag].bag;
        T[refinement_bag].bag = S;

        std::vector<std::set<unsigned int> > union_S_W_i(components.size());

        for(unsigned int i = 0; i < components.size(); i++){
            std::set_union(S.begin(), S.end(), components[i].begin(), components[i].end(), std::inserter(union_S_W_i[i], union_S_W_i[i].begin()));
            newN[i] = boost::add_vertex(T);
            T[newN[i]].bag = union_S_W_i[i];
            boost::add_edge(refinement_bag, newN[i], T);
        }

        for(unsigned int i = 0; i <  oldN.size(); i++){
            std::set<unsigned int> intersection;
            std::set_intersection(old_bag.begin(), old_bag.end(), T[oldN[i]].bag.begin(), T[oldN[i]].bag.end(), std::inserter(intersection, intersection.begin()));
            for(unsigned int j = 0; j < union_S_W_i.size(); j++){
                if(std::includes(union_S_W_i[j].begin(), union_S_W_i[j].end(), intersection.begin(), intersection.end())){
                    boost::add_edge(newN[j], oldN[i], T);
                    break;
                }
            }
        }
    }
}

template <typename G_t>
bool is_candidate_edge(std::vector<unsigned int> &edge, unsigned int i, std::vector<unsigned int> &elimination_ordering, G_t &M, std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &idxMap){
    //Position i in 'elimination_ordering_' will store the 'elimination date' of vertex i
    std::vector<unsigned int> elimination_ordering_(elimination_ordering.size());
    for(unsigned int t = 0; t < elimination_ordering.size(); t++)
        elimination_ordering_[elimination_ordering[t]] = t;

    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(idxMap[edge[0]], M); nIt != nEnd; nIt++){
        if(elimination_ordering_[M[*nIt].id] > i && boost::edge(idxMap[edge[1]], *nIt, M).second && !boost::edge(*nIt, idxMap[elimination_ordering[i]], M).second)
            return false;
    }
    return true;
}



/* minimalChordal-algorithm
 *
 * Computes possibly redundant fill-in-edges and runs LEX-M to check,
 * if the graph after removal of a fill-in-edge is chordal.
 * Finally, the algorithm computes a new perfect elimination ordering, that
 * possibly causes lower width than 'old_elimination_ordering'.
 */
template <typename G_t>
void minimalChordal(G_t &G, std::vector<unsigned int> &old_elimination_ordering, typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &new_elimination_ordering){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    make_index_map(G, idxMap);

    std::vector<std::set<unsigned int> > C;
    std::vector<std::vector<std::vector<unsigned int> > > F;
    make_filled_graph(G, old_elimination_ordering, C, F);

    for(int i = old_elimination_ordering.size()-1; i >= 0; i--){
        std::vector<std::vector<unsigned int> > candidate;
        std::set<unsigned int> incident;
        for(unsigned int j = 0; j < F[i].size(); j++){
            if(is_candidate_edge(F[i][j], i, old_elimination_ordering, G, idxMap)){
                candidate.push_back(F[i][j]);
                incident.insert(F[i][j][0]);
                incident.insert(F[i][j][1]);
            }
        }
        if(candidate.size() != 0){
            G_t W_i;
            induced_subgraph(W_i, G, incident);
            delete_edges(W_i, candidate);

            std::vector<std::vector<unsigned int> > keep_fill;
            LEX_M_fill_in(W_i, keep_fill);

            for(unsigned int j = 0; j < candidate.size(); j++){
                for(unsigned int k = 0; k < keep_fill.size(); k++){
                    if((candidate[j][0] == keep_fill[k][0] && candidate[j][1] == keep_fill[k][1])
                     ||(candidate[j][0] == keep_fill[k][1] && candidate[j][1] == keep_fill[k][0])){
                         candidate.erase(candidate.begin()+j);
                         break;
                     }
                }
            }

            delete_edges(G, candidate);
        }
    }
    LEX_M_minimal_ordering(G, new_elimination_ordering);
}

}

#endif
