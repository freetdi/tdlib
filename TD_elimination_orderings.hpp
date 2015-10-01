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
// Offers functionality to compute tree decompositions of graphs
// according to various heuristic and functions, that convert
// tree decompositions to elimination orderings and vice versa.
// Also the LEX-M algorithm is included in this header
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
// void minDegree_decomp(G_t G&, T_t &T)
// void fillIn_decomp(G_t G&, T_t &T)
// void minDegree_ordering(G_t G, std::vector<unsigned int> &elim_ordering)
// void fillIn_ordering(G_t G, std::vector<unsigned int> &elim_ordering)
// void make_filled_graph(G_t &G, std::vector<unsigned int> &elim_ordering, std::vector<std::set<unsigned int> > &C, std::vector<std::vector<std::vector<unsigned int> > > &F)
// void ordering_to_treedec(G_t G, std::vector<unsigned int> &elimination_ordering, T_t &T)
// void treedec_to_ordering(T_t T, std::vector<unsigned int> &elimination_ordering)
// void LEX_M_fill_in(G_t &G, std::vector<std::vector<unsigned int> > &fill_in_edges)
// void LEX_M_minimal_ordering(G_t &G, std::vector<unsigned int> &alpha)
//

#ifndef ELIMINATION_ORDERING
#define ELIMINATION_ORDERING

#include <cmath>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include "simple_graph_algos.hpp"
#include "TD_misc.hpp"

namespace treedec{

#ifndef GLUE_BAG
#define GLUE_BAG

//glues a single bag with the current tree decomposition
template <typename T_t>
void glue_bag(std::set<unsigned int> &bag, unsigned int elim_vertex, T_t &T){
    typename boost::graph_traits<T_t>::vertex_descriptor t_dec_node;
    if(boost::num_vertices(T) == 0){
        bag.insert(elim_vertex);
        t_dec_node = boost::add_vertex(T);
        T[t_dec_node].bag = bag;
        bag.clear();
        return;
    }
        
    typename boost::graph_traits<T_t>::vertex_iterator vIt, vEnd;
    
    for(boost::tie(vIt, vEnd) = boost::vertices(T); vIt != vEnd; vIt++){
        std::set<unsigned int>::iterator sIt = bag.begin();
        if(std::includes(T[*vIt].bag.begin(), T[*vIt].bag.end(), sIt, bag.end())){
            bag.insert(elim_vertex);
            t_dec_node = boost::add_vertex(T);
            T[t_dec_node].bag = bag;
            boost::add_edge(*vIt, t_dec_node, T);
            bag.clear();
            return;
        }
    }
    t_dec_node = boost::add_vertex(T);
    bag.insert(elim_vertex);
    T[t_dec_node].bag = bag;
    boost::tie(vIt, vEnd) = boost::vertices(T);
    boost::add_edge(*vIt, t_dec_node, T);
    bag.clear();
}

#endif

#ifndef REMOVE_ISOLATED_VERTICES
#define REMOVE_ISOLATED_VERTICES

template <typename G_t>
void remove_isolated_vertices(G_t &G){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    bool exists_isolated = true;
    while(exists_isolated){
        exists_isolated = false;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){     
            if(boost::out_degree(*vIt, G) == 0){
                boost::remove_vertex(*vIt, G);
                exists_isolated = true;
                break;
            }
        }
    }
}

#endif

//constructs a tree decomposition from the elimination ordering obtained by the
//minimum-degree heuristic 
template <typename G_t, typename T_t>
void minDegree_decomp(G_t &G, T_t &T){
    std::vector<std::set<unsigned int> > bags;
    std::vector<unsigned int> elim_vertices;

    while(boost::num_edges(G) > 0){
        //search a minimum degree vertex
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;
        unsigned int min_degree = boost::num_vertices(G);
    
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int degree = boost::out_degree(*vIt, G);
            if(degree != 0 && degree < min_degree){
                min_degree = degree;
                min_vertex = *vIt;
            }
        }
    
        //collect the neighbours
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> neighbours;
        std::set<unsigned int> bag;
        
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;    

        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, G); nIt != nEnd; nIt++){
            bag.insert(G[*nIt].id);
            neighbours.push_back(*nIt);
        }

	//make the neighbours a clique
        for(unsigned int i = 0; i < neighbours.size(); i++){
            for(unsigned int j = i+1; j < neighbours.size(); j++){
                boost::add_edge(neighbours[i], neighbours[j], G);	
            }
        }

        bags.push_back(bag);
        elim_vertices.push_back(G[min_vertex].id);

        boost::clear_vertex(min_vertex, G);
    }
    
    for(unsigned int i = bags.size(); i > 0; i--)
        glue_bag(bags[i-1], elim_vertices[i-1], T);
}

//constructs a tree decomposition from the elimination ordering obtained by the
//fill-in heuristic 
template <typename G_t, typename T_t>
void fillIn_decomp(G_t &G, T_t &T){
    if(boost::num_edges(G) == 0)
        return;
    
    //search a vertex with least non-adjacent neighbours
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;

    boost::tie(vIt, vEnd) = boost::vertices(G);
    typename boost::graph_traits<G_t>::vertex_descriptor min_vertex = *vIt;

    unsigned int min_fill = boost::num_vertices(G);
    for(; vIt != vEnd; vIt++){
        if(boost::out_degree(*vIt, G) == 0)
            continue;

        unsigned int current_fill = 0;
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++)
            N.insert(*nIt);

        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt1, sIt2;
        for(sIt1 = N.begin(); sIt1 != N.end(); sIt1++){
            sIt2 = sIt1;
            sIt2++;
            for(; sIt2 != N.end(); sIt2++){
                std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(*sIt1, *sIt2, G);
                if(!existsEdge.second)
                    current_fill++;
            }
        }
        if(current_fill < min_fill){
            min_fill = current_fill;
            min_vertex = *vIt;
            if(current_fill == 0)
                break;
        }
    }
    
    //collect the neighbours
    std::set<unsigned int> bag;
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> neighbours;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, G); nIt != nEnd; nIt++){
        bag.insert(G[*nIt].id);
        neighbours.push_back(*nIt);
    }
    
    //make the neighbours a clique
    for(unsigned int i = 0; i < neighbours.size(); i++){
        for(unsigned int j = i+1; j < neighbours.size(); j++)
            boost::add_edge(neighbours[i], neighbours[j], G);
    }
        
    unsigned int elim_vertex = G[min_vertex].id;
    
    boost::clear_vertex(min_vertex, G);
    
    fillIn_decomp(G, T);
    glue_bag(bag, elim_vertex, T);
}

//computes an elimination ordering according to minDegree heuristic (version used for postprocessing algorithms)
template<typename G_t>
void minDegree_ordering(G_t G, std::vector<unsigned int> &elim_ordering){
    while(boost::num_vertices(G) != 0){
        //search a minimum degree vertex
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex = NULL;
        unsigned int min_degree = boost::num_vertices(G);
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int degree = boost::out_degree(*vIt, G);
            if(degree < min_degree){
                min_degree  = degree;
                min_vertex = *vIt;
            }
        }
        
        //collect the neighbours
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> neighbours;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, G); nIt != nEnd; nIt++)
            neighbours.push_back(*nIt);
        
        //make the neighbours a clique
        std::vector<std::vector<unsigned int> > F_i;
        for(unsigned int i = 0; i < neighbours.size(); i++){
            for(unsigned int j = i+1; j < neighbours.size(); j++){
                std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(neighbours[i], neighbours[j], G);
                if(!existsEdge.second)
                    boost::add_edge(neighbours[i], neighbours[j], G);
            }
        }
            
        elim_ordering.push_back(G[min_vertex].id);
    
        boost::clear_vertex(min_vertex, G);
        boost::remove_vertex(min_vertex, G);
    }
}

//computes an elimination ordering according to minDegree heuristic (version used for postprocessing algorithms)
template<typename G_t>
void fillIn_ordering(G_t G, std::vector<unsigned int> &elim_ordering){
    while(boost::num_vertices(G) != 0){
        //search a vertex with least non-adjacent neighbours
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;
        unsigned int min_fill = boost::num_vertices(G);
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int current_fill = 0;
            std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++)
                N.insert(*nIt);
            
            typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt1, sIt2;
            for(sIt1 = N.begin(); sIt1 != N.end(); sIt1++){
                sIt2 = sIt1;
                sIt2++;
                for(; sIt2 != N.end(); sIt2++){
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(*sIt1, *sIt2, G);
                    if(!existsEdge.second)
                        current_fill++;
                }
            }
            if(current_fill < min_fill){
                min_fill = current_fill;
                min_vertex = *vIt;
                if(current_fill == 0)
                    break;
            }
        }
        
        //collect the neighbours
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> neighbours;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, G); nIt != nEnd; nIt++)
            neighbours.push_back(*nIt);
        
    
        //make the neighbours a clique
        for(unsigned int i = 0; i < neighbours.size(); i++){
            for(unsigned int j = i+1; j < neighbours.size(); j++){
                std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(neighbours[i], neighbours[j], G);
                if(!existsEdge.second)
                    boost::add_edge(neighbours[i], neighbours[j], G);
            }
        }
            
        elim_ordering.push_back(G[min_vertex].id);
    
        boost::clear_vertex(min_vertex, G);
        boost::remove_vertex(min_vertex, G);
    }
}

//makes G a filled graph according to the additional edges F
template <typename G_t>
void make_filled_graph(G_t &G, std::vector<unsigned int> &elim_ordering, std::vector<std::set<unsigned int> > &C, std::vector<std::vector<std::vector<unsigned int> > > &F){
    //evaluate maximum id
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    unsigned int max = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        max = (G[*vIt].id > max)? G[*vIt].id : max;
    
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap(max+1); 
    
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        idxMap[G[*vIt].id] = *vIt; 

    for(unsigned int i = 0; i < elim_ordering.size(); i++){
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        std::set<unsigned int> N_i, E_i, intersection;
        N_i.insert(elim_ordering[i]);
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(idxMap[elim_ordering[i]], G); nIt != nEnd; nIt++)
            N_i.insert(G[*nIt].id);
        for(unsigned int j = i; j < elim_ordering.size(); j++)
            E_i.insert(elim_ordering[j]);
        
        std::vector<std::vector<unsigned int> > F_i;
        std::set_intersection(N_i.begin(), N_i.end(), E_i.begin(), E_i.end(), std::inserter(intersection, intersection.begin()));
        C.push_back(intersection);
        
        for(std::set<unsigned int>::iterator sIt1 = intersection.begin(); sIt1 != intersection.end(); sIt1++){
            std::set<unsigned int>::iterator sIt2 = sIt1;
            sIt2++;
            for(; sIt2 != intersection.end(); sIt2++){
                std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(idxMap[*sIt1], idxMap[*sIt2], G);
                if(!existsEdge.second){
                    std::vector<unsigned int> edge;
                    edge.push_back(*sIt1);
                    edge.push_back(*sIt2);
                    F_i.push_back(edge);
                    boost::add_edge(idxMap[*sIt1], idxMap[*sIt2], G);
                }
            }
        }
        
        F.push_back(F_i);
    }
}

template <typename G_t, typename T_t>
void _ordering_to_treedec(G_t &G, std::vector<unsigned int> &elimination_ordering, T_t &T, unsigned int idx){
    if(idx == elimination_ordering.size())
        return;
    
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::vertex_descriptor elim_vertex;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        if(G[*vIt].id == elimination_ordering[idx]){
            elim_vertex = *vIt;
            break;
        }
    }
        
    //collect the neighbours
    std::set<unsigned int> bag;
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> neighbours;
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(elim_vertex, G); nIt != nEnd; nIt++){
        bag.insert(G[*nIt].id);
        neighbours.push_back(*nIt);
    }
    
    //make the neighbours a clique
    for(unsigned int i = 0; i < neighbours.size(); i++){
        for(unsigned int j = i+1; j < neighbours.size(); j++)
            boost::add_edge(neighbours[i], neighbours[j], G);
    }
    
    unsigned int elim_vertex_id = G[elim_vertex].id;
        
    boost::clear_vertex(elim_vertex, G);
    
    _ordering_to_treedec(G, elimination_ordering, T, idx+1);
    glue_bag(bag, elim_vertex_id, T);
}

template <typename G_t, typename T_t>
void ordering_to_treedec(G_t G, std::vector<unsigned int> &elimination_ordering, T_t &T){
    _ordering_to_treedec(G, elimination_ordering, T, 0);
}

template <typename T_t>
void _treedec_to_ordering(T_t &T, std::vector<unsigned int> &elimination_ordering){
    if(boost::num_vertices(T) == 0)
        return;
    
    typename boost::graph_traits<T_t>::vertex_iterator tIt, tEnd;
    typename boost::graph_traits<T_t>::vertex_descriptor leaf, parent;
    for(boost::tie(tIt, tEnd) = boost::vertices(T); tIt != tEnd; tIt++){
        if(boost::out_degree(*tIt, T) <= 1){
            leaf = *tIt;
            break;
        }
    }
    
    if(T[leaf].bag.size() == 0){
        boost::clear_vertex(leaf, T);
        boost::remove_vertex(leaf, T);
    }
    else{
        typename boost::graph_traits<T_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(leaf, T); nIt != nEnd; nIt++)
            parent = *nIt;
        
        std::set<unsigned int> difference;
        
        if(boost::out_degree(leaf, T) != 0){
            if(std::includes(T[parent].bag.begin(), T[parent].bag.end(), T[leaf].bag.begin(), T[leaf].bag.end())){
                boost::clear_vertex(leaf, T);
                boost::remove_vertex(leaf, T);
            }
            else
                std::set_difference(T[leaf].bag.begin(), T[leaf].bag.end(), T[parent].bag.begin(), T[parent].bag.end(), std::inserter(difference, difference.begin()));
        }
        else
            difference = T[leaf].bag;
            
        for(std::set<unsigned int>::iterator sIt = difference.begin(); sIt != difference.end(); sIt++){
            elimination_ordering.push_back(*sIt);
            T[leaf].bag.erase(*sIt);
        }
    }
    
    _treedec_to_ordering(T, elimination_ordering);
}

template <typename T_t>
void treedec_to_ordering(T_t T, std::vector<unsigned int> &elimination_ordering){
    _treedec_to_ordering(T, elimination_ordering);
}

template <typename G_t>
void LEX_M_fill_in(G_t &G, std::vector<std::vector<unsigned int> > &fill_in_edges){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    unsigned int max = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        max = (G[*vIt].id > max)? G[*vIt].id : max;
    
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    make_index_map(G, idxMap);
    
    std::vector<bool> visited(max+1);
    std::vector<float> label(max+1);
    std::vector<unsigned int> alpha_inv(max+1);
    std::vector<std::vector<unsigned int> > reached_i(max+1);
    
    unsigned int i = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        label[G[*vIt].id] = 1.0;
        alpha_inv[i++] = 0;
        visited[G[*vIt].id] = false;
    }
    
    unsigned int k = 1;
    
    for(int i = boost::num_vertices(G)-1; i >= 0; i--){
        typename boost::graph_traits<G_t>::vertex_descriptor v;
        unsigned int max = 0;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            if(alpha_inv[G[*vIt].id] == 0){
                if(label[G[*vIt].id] > max){
                    max = (unsigned int) label[G[*vIt].id];
                    v = *vIt;
                }
            }
        }
        visited[G[v].id] = true; 
        alpha_inv[G[v].id] = i+1;
        
        for(unsigned int j = 0; j < k; j++)
            reached_i[j].clear();
        
        for(unsigned int j = 0; j < alpha_inv.size(); j++){
            if(alpha_inv[j] == 0)
                visited[j] = false;
        }
        
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; nIt++){
            if(alpha_inv[G[*nIt].id] == 0){
                reached_i[(int)label[G[*nIt].id]-1].push_back(G[*nIt].id);
                visited[G[*nIt].id] = true;
                label[G[*nIt].id] += 0.5;
            }
        }

        for(unsigned int j = 0; j < k; j++){
            while(reached_i[j].size() != 0){
                unsigned int w = reached_i[j].back();

                reached_i[j].pop_back();
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(idxMap[w], G); nIt != nEnd; nIt++){
                    if(visited[G[*nIt].id])
                        continue;

                    visited[G[*nIt].id] = true;
                    if((unsigned int)label[G[*nIt].id]-1 > j){
                        reached_i[(int)label[G[*nIt].id]].push_back(G[*nIt].id);
                        label[G[*nIt].id] += 0.5;
                        std::vector<unsigned int> edge;
                        edge.push_back(G[v].id);
                        edge.push_back(G[*nIt].id);
                        fill_in_edges.push_back(edge);
                    }
                    else{
                        reached_i[j].push_back(G[*nIt].id);
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
void LEX_M_minimal_ordering(G_t &G, std::vector<unsigned int> &alpha){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    unsigned int max = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        max = (G[*vIt].id > max)? G[*vIt].id : max;
    
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> idxMap;
    treedec::make_index_map(G, idxMap);
    
    alpha.resize(boost::num_vertices(G));
    
    std::vector<bool> visited(max+1);
    std::vector<float> label(max+1);
    std::vector<unsigned int> alpha_inv(max+1);
    std::vector<std::vector<unsigned int> > reached_i(max+1);
    
    unsigned int i = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        label[G[*vIt].id] = 1.0;
        alpha_inv[i++] = 0;
        visited[G[*vIt].id] = false;
    }
    
    unsigned int k = 1;
    
    for(int i = boost::num_vertices(G)-1; i >= 0; i--){
        typename boost::graph_traits<G_t>::vertex_descriptor v;
        unsigned int max = 0;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            if(alpha_inv[G[*vIt].id] == 0){
                if((unsigned int)label[G[*vIt].id] > max){
                    max = (unsigned int) label[G[*vIt].id];
                    v = *vIt;
                }
            }
        }
        visited[G[v].id] = true; 
        alpha[i] = G[v].id;
        alpha_inv[G[v].id] = i+1;
        
        for(unsigned int j = 0; j < k; j++)
            reached_i[j].clear();
        
        for(unsigned int j = 0; j < alpha_inv.size(); j++){
            if(alpha_inv[j] == 0)
                visited[j] = false;
        }
        
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, G); nIt != nEnd; nIt++){
            if(alpha_inv[G[*nIt].id] == 0){
                reached_i[(int)label[G[*nIt].id]-1].push_back(G[*nIt].id);
                visited[G[*nIt].id] = true;
                label[G[*nIt].id] += 0.5;
            }
        }

        for(unsigned int j = 0; j < k; j++){
            while(reached_i[j].size() != 0){
                unsigned int w = reached_i[j].back();

                reached_i[j].pop_back();
                for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(idxMap[w], G); nIt != nEnd; nIt++){
                    if(visited[G[*nIt].id])
                        continue;

                    visited[G[*nIt].id] = true;
                    if((unsigned int)label[G[*nIt].id]-1 > j){
                        reached_i[(int)label[G[*nIt].id]].push_back(G[*nIt].id);
                        label[G[*nIt].id] += 0.5;
                    }
                    else{
                        reached_i[j].push_back(G[*nIt].id);
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

}

#endif
