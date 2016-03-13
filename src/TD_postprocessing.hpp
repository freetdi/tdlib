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
// typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS, Vertex> TD_graph_t;
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
#include "TD_noboost.hpp"

namespace treedec{

//Creates a modified induced subgraph of the bag 'noboost::bag(T, t_desc)'.
template <typename G_t, typename T_t>
bool is_improvement_bag(G_t &H,
                        std::vector<bool> &disabled,
                        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &X,
                        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &Y,
                        typename boost::graph_traits<T_t>::vertex_descriptor t_desc,
                        typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &vdMap,
                        G_t &G, T_t &T)
{
        std::cout << "bag:" << std::endl;
        for(typename noboost::treedec_traits<T_t>::bag_type::iterator sIt = noboost::bag(T, t_desc).begin(); sIt != noboost::bag(T, t_desc).end(); sIt++){
            std::cout << *sIt << " ";
        } std::cout << std::endl;

    treedec::induced_subgraph(H, G, noboost::bag(T, t_desc), vdMap);
    disabled.assign(boost::num_vertices(H), false);

    std::cout << "H_size: " << boost::num_vertices(H) << std::endl;

    std::cout << "edges: " << boost::num_edges(H) << " " << boost::num_edges(G) << std::endl;

    //Add an additional edge, if a non-edge 'occures' in a bag of an adjacent vertex t_desc' of t_desc in T.
    typename boost::graph_traits<G_t>::vertex_iterator vIt1, vIt2, vEnd;
    for(boost::tie(vIt1, vEnd) = boost::vertices(H); vIt1 != vEnd; vIt1++){
        vIt2 = vIt1;
        vIt2++;
        for(; vIt2 != vEnd; vIt2++){
            if(!boost::edge(*vIt1, *vIt2, H).second){
                typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(t_desc, T); nIt != nEnd; nIt++){
                    unsigned int pos1 = noboost::get_pos(*vIt1, H);
                    unsigned int pos2 = noboost::get_pos(*vIt2, H);
                    typename noboost::treedec_traits<T_t>::bag_type::value_type vd1 = vdMap[pos1];
                    typename noboost::treedec_traits<T_t>::bag_type::value_type vd2 = vdMap[pos2];
                    if(noboost::bag(T, *nIt).find(vd1) != noboost::bag(T, *nIt).end()
                    && noboost::bag(T, *nIt).find(vd2) != noboost::bag(T, *nIt).end())
                    {
                        std::cout << "add an edge" << std::endl;
                        boost::add_edge(*vIt1, *vIt2, H);
                        break;
                    }
                }
            }
        }
    }

    //Find a non-edge {x,y} and collect the neighbours of x and y, resulting in the sets X and Y.
    for(boost::tie(vIt1, vEnd) = boost::vertices(H); vIt1 != vEnd; vIt1++){
        std::cout << "!" << std::endl;
        vIt2 = vIt1;
        vIt2++;
        for(; vIt2 != vEnd; vIt2++){
            std::cout << "." << std::endl;
            if(!boost::edge(*vIt1, *vIt2, H).second){
                typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt1, H); nIt != nEnd; nIt++){
                    X.insert(*nIt);
                }

                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt2, H); nIt != nEnd; nIt++){
                    Y.insert(*nIt);
                }

                unsigned int pos1 = noboost::get_pos(*vIt1, H);
                unsigned int pos2 = noboost::get_pos(*vIt2, H);
                disabled[pos1] = true;
                disabled[pos2] = true;

                std::cout << "disable: " << pos1 << " " << pos2 << std::endl;

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
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> X, Y;
        std::vector<bool> disabled;
        typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> vdMap;

        typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
        typename boost::graph_traits<T_t>::vertex_descriptor refinement_bag;
        for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
            if(noboost::bag(T, *tIt).size() == width+1){
                std::vector<bool> disabled_;
                typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> vdMap_;
                if(is_improvement_bag(H, disabled_, X, Y, *tIt, vdMap_, G, T)){
                    refinement_bag = *tIt;
                    disabled = MOVE(disabled_);
                    vdMap = MOVE(vdMap_);
                    break;
                }
            }
        }

        //No improvement possible.
        if(boost::num_vertices(H) == 0){
            return;
        }

        std::cout << "X:" << std::endl;
        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = X.begin(); sIt != X.end(); sIt++){
            std::cout << *sIt << " ";
        } std::cout << std::endl;

        std::cout << "Y:" << std::endl;
        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = Y.begin(); sIt != Y.end(); sIt++){
            std::cout << *sIt << " ";
        } std::cout << std::endl;

        std::cout << "enabled: " << std::endl;
        for(unsigned int i = 0; i < disabled.size(); i++){
            if(!disabled[i]){
                std::cout << i << " ";
            }
        } std::cout << std::endl;

        //Compute a seperating set S.
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> S_;
        treedec::seperate_vertices(H, disabled, X, Y, S_);

        //Do the refinement.
        std::vector<bool> visited(boost::num_vertices(H), false);

        //S consists of vertex descriptors of H. Use vd_map to map these to descriptors of G.
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> S;
        treedec::map_descriptors(S_, S, H, vdMap);

        if(S.size() == 0){ std::cout << "error: no seperator found" << std::endl; }

        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = S.begin(); sIt != S.end(); sIt++){
            unsigned int pos = noboost::get_pos(*sIt, H);
            visited[pos] = true;
        }

        typename noboost::treedec_traits<T_t>::bag_type old_bag = noboost::bag(T, refinement_bag);
        boost::clear_vertex(refinement_bag, T);

        //Convert the descriptors of G to a bag.
        typename noboost::treedec_traits<T_t>::bag_type B;
        treedec::map_descriptors_to_bags<G_t>(S, B);
        noboost::bag(T, refinement_bag) = MOVE(B);

        std::cout << "seperator:" << std::endl;
        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = S.begin(); sIt != S.end(); sIt++){
            std::cout << H[*sIt].id << std::endl;
        }

        std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > components;
        treedec::get_components_provided_map(H, components, visited);

        typename boost::graph_traits<T_t>::adjacency_iterator t_nIt, t_nEnd;
        std::vector<typename boost::graph_traits<T_t>::vertex_descriptor> oldN(boost::degree(refinement_bag, T));;
        std::vector<typename boost::graph_traits<T_t>::vertex_descriptor> newN(components.size());

        unsigned int c = 0;
        for(boost::tie(t_nIt, t_nEnd) = boost::adjacent_vertices(refinement_bag, T); t_nIt != t_nEnd; t_nIt++){
            oldN[c++] = *t_nIt;
        }

        std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > union_S_W_i(components.size());

        for(unsigned int i = 0; i < components.size(); i++){
            std::set<typename boost::graph_traits<G_t>::vertex_descriptor> component;
            treedec::map_descriptors(components[i], component, H, vdMap);

            std::cout << "comp_" << i << ": " << std::endl;
            for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt 
                 = component.begin(); sIt != component.end(); sIt++)
            { std::cout << G[*sIt].id << " "; } std::cout << std::endl;


            std::set_union(S.begin(), S.end(),
                           component.begin(), component.end(),
                           std::inserter(union_S_W_i[i], union_S_W_i[i].begin()));
            newN[i] = boost::add_vertex(T);

            typename noboost::treedec_traits<T_t>::bag_type uB;
            treedec::map_descriptors_to_bags<G_t>(union_S_W_i[i], uB);

            noboost::bag(T, newN[i]) = MOVE(uB);
            boost::add_edge(refinement_bag, newN[i], T);
        }

        for(unsigned int i = 0; i <  oldN.size(); i++){
            std::set<typename boost::graph_traits<G_t>::vertex_descriptor> intersection;
            std::set_intersection(old_bag.begin(), old_bag.end(),
                                  noboost::bag(T, oldN[i]).begin(),
                                  noboost::bag(T, oldN[i]).end(),
                                  std::inserter(intersection, intersection.begin()));

            for(unsigned int j = 0; j < union_S_W_i.size(); j++){
                if(std::includes(noboost::bag(T, newN[j]).begin(), noboost::bag(T, newN[j]).end(),
                                 intersection.begin(), intersection.end())){
                    boost::add_edge(newN[j], oldN[i], T);
                    break;
                }
            }
        }

        G_t L;
        boost::copy_graph(G, L);
        if (treedec::is_valid_treedecomposition(L, T) < 0){ std::cout << "invalid decomposition" << std::endl;
            for(unsigned int l = 0; l < boost::num_vertices(T); l++){
                for(typename noboost::treedec_traits<T_t>::bag_type::iterator sIt = T[l].bag.begin(); sIt != T[l].bag.end(); sIt++){
                    std::cout << G[*sIt].id << " ";
                } std::cout << std::endl;
            }
            write_dot_td("tmptd.dot", T);
            exit(-1); 
        }

    }
}

template <typename G_t>
bool is_candidate_edge(std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &edge, unsigned int i,
                       std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering, G_t &M)
{
    //Position i in 'elimination_ordering_' will store the 'elimination date' of vertex i
    std::vector<unsigned int> elimination_ordering_(elimination_ordering.size());
    for(unsigned int t = 0; t < elimination_ordering.size(); t++){
        unsigned int pos = noboost::get_pos(elimination_ordering[t], M);
        elimination_ordering_[pos] = t;
    }

    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(edge[0], M); nIt != nEnd; nIt++){
        unsigned int pos = noboost::get_pos(*nIt, M);
        if(elimination_ordering_[pos] > i && boost::edge(edge[1], *nIt, M).second
       && !boost::edge(*nIt, elimination_ordering[i], M).second)
        {
            return false;
        }
    }
    return true;
}

template <typename G_t>
void delete_edges(G_t &G, std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > &edges){
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
    std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > C;
    std::vector<std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > > F;

    treedec::make_filled_graph(G, old_elimination_ordering, C, F);

    for(int i = old_elimination_ordering.size()-1; i >= 0; i--){
        std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > candidate;
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> incident;
        for(unsigned int j = 0; j < F[i].size(); j++){
            if(is_candidate_edge(F[i][j], i, old_elimination_ordering, G)){
                candidate.push_back(F[i][j]);
                incident.insert(F[i][j][0]);
                incident.insert(F[i][j][1]);
            }
        }
        if(candidate.size() != 0){
            G_t W_i;
            typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> vdMap;
            treedec::induced_subgraph(W_i, G, incident, vdMap);
            treedec::delete_edges(W_i, candidate);

            std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > keep_fill_;
            std::cout << "1" << std::endl;
            treedec::LEX_M_fill_in(W_i, keep_fill_);
            std::cout << "2" << std::endl;

            //Translate descriptors of W_i to descriptors of G.
            std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > keep_fill(keep_fill_.size());
            for(unsigned int j = 0; j < keep_fill_.size(); j++){
                unsigned int pos1 = noboost::get_pos(keep_fill_[j][0], W_i);
                unsigned int pos2 = noboost::get_pos(keep_fill_[j][1], W_i);
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
            std::cout << "3" << std::endl;
    treedec::LEX_M_minimal_ordering(G, new_elimination_ordering);
            std::cout << "4" << std::endl;
}

} //namespace treedec

#endif //ifdef TD_POSTPROCESSING

// vim:ts=8:sw=4:et:
