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
// Offers functionality to compute lower bounds on tree width
//

/*
 * provides following functions (namespace treedec::lb):
 * 
 * - int delta(G_t &G)
 * - int delta2(G_t &G)
 * - int gamma(G_t &G)
 * - int deltaD(G_t G)
 * - int delta2D(G_t &G)
 * - int gammaD_left(G_t &G)
 * - int gammaD_right(G_t &G)
 * - int gammaD_min_e(G_t &G)
 * - int deltaC_min_d(G_t &G)
 * - int deltaC_max_d(G_t &G)
 * - int deltaC_least_c(G_t &G)
 * 
 * - void k_neighbour_improved_graph(G_t &G, unsigned int k)
 * - int LBN_deltaD(G_t &G)
 * - int LBN_deltaC(G_t &G)
 * - int LBNC_deltaD(G_t &G)
 * - int LBNC_deltaC(G_t &G)
 * - void k_path_improved_graph(G_t &G, unsigned int k)
 * - int LBP_deltaD(G_t &G)
 * - int LBP_deltaC(G_t &G)
 * - int LBPC_deltaD(G_t &G)
 * - int LBPC_deltaC(G_t &G)
 * 
 * - void MCS_random(G_t &G, int &lb)
 * - std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> MCS_all_start_vertices(G_t &G, int &lb)
 * - void MCSC_min_deg(G_t G, int &lb)
 * - void MCSC_last_mcs(G_t G, int &lb)
 * - void MCSC_max_mcs(G_t G, int &lb)
 * 
 * - int relation_edges_vertices(G_t &G)
 */

#ifndef TD_LOWER_BOUNDS
#define TD_LOWER_BOUNDS

#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include "TD_NetworkFlow.hpp"
#include "simple_graph_algos.hpp"

namespace treedec{
    
namespace lb{



/* DEGREE BASED */

//smallest vertex-degree in G
template <typename G_t>
int delta(G_t &G){
    if(boost::num_vertices(G) == 0)
        return -1;

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    unsigned int min = boost::num_vertices(G);
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::out_degree(*vIt, G);
        min = (degree < min)? degree : min;
    }
    return (int)min;
}

//second smallest vertex-degree in G
template <typename G_t>
int delta2(G_t &G){
    if(boost::num_vertices(G) == 0)
        return -1;

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::vertex_descriptor v_min;
    unsigned int min = boost::num_vertices(G);
    unsigned int snd = boost::num_vertices(G);
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::out_degree(*vIt, G);
        if(degree <= min){
            snd = min;
            min = degree;
        }
        if(degree > min && degree < snd)
            snd = degree;
    }
    return (int)snd;
}

//maximum degree of minimal degree of not with an edge incident vertices
template <typename G_t>
int gamma(G_t &G){
    if(boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return boost::num_vertices(G)-1;
    }
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    
    //sort the vertices of G according to rising degree -> degree_sequence
    unsigned int max_degree = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::out_degree(*vIt, G);
        max_degree = (degree>max_degree)? degree : max_degree;
    }
    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > buckets(max_degree+1);
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> degree_sequence;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::out_degree(*vIt, G);
        buckets[degree].push_back(*vIt);
    }
    for(unsigned int i = 0; i <= max_degree; i++){
        for(unsigned int j = 0; j < buckets[i].size(); j++){
            degree_sequence.push_back(buckets[i][j]);
        }
    }
    //take the degree of the right vertex in the first not-edge
    for(unsigned int i = 0; i < boost::num_vertices(G); i++){
        for(unsigned int j = 0; j < i; j++){
            std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(degree_sequence[i], degree_sequence[j], G);
            if(!existsEdge.second){
                unsigned int deg = boost::out_degree(degree_sequence[i], G);
                return (int)deg;
            }
        }
    }
}

//successivly remove a smallest-degree-vertex in G and return the smallest degree seen
template <typename G_t>
int deltaD(G_t &G){
    if(boost::num_vertices(G) == 0)
        return -1;

    G_t H;
    TD_copy_graph(G, H);

    unsigned int maxmin = 0;
    typename boost::graph_traits<G_t>::vertex_iterator hIt, hEnd;
    typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;
    
    while(boost::num_vertices(H) > 0){
        unsigned int min_degree = boost::num_vertices(H);
        for(boost::tie(hIt, hEnd) = vertices(H); hIt != hEnd; hIt++){
            unsigned int degree = boost::out_degree(*hIt, H);
            if(degree < min_degree){
                min_degree = degree;
                min_vertex = *hIt;
            }
        }

        maxmin = (maxmin>min_degree)? maxmin : min_degree;
        
        boost::clear_vertex(min_vertex, H);
        boost::remove_vertex(min_vertex, H);
    }
    return (int)maxmin;
}

//assume each vertex as the minimal one and do deltaD
template <typename G_t>
int delta2D(G_t &G){
    if(boost::num_vertices(G) == 0)
        return -1;

    G_t H;
    TD_copy_graph(G, H);
    
    typename std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> assumed_minimal;
    
    typename boost::graph_traits<G_t>::vertex_iterator hIt, hEnd;
    for(boost::tie(hIt, hEnd) = vertices(H); hIt != hEnd; hIt++)
        assumed_minimal.push_back(*hIt);
    
    unsigned int min_degree, maxmin = 0;
    typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;
    
    for(unsigned int i = 0; i < assumed_minimal.size(); i++){
        while(boost::num_vertices(H) >= 2){
            min_degree = boost::num_vertices(H)+1;

            for(boost::tie(hIt, hEnd) = vertices(H); hIt != hEnd; hIt++){
                if(*hIt == assumed_minimal[i])
                    continue;
                unsigned int degree = boost::out_degree(*hIt, H);
                if(degree < min_degree){
                    min_degree = degree;
                    min_vertex = *hIt;
                }
            }
            maxmin = (maxmin>min_degree)? maxmin : min_degree;
            boost::clear_vertex(min_vertex, H);
            boost::remove_vertex(min_vertex, H);    
        }
        H.clear();
        TD_copy_graph(G, H);
    }
    return (int)maxmin;
}

template <typename G_t>
void _gammaD_left(G_t &G, unsigned int &lb){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    
    //sort the vertices of G according to rising degree -> degree_sequence
    unsigned int max_degree = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::out_degree(*vIt, G);
        max_degree = (degree>max_degree)? degree : max_degree;
    }
    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > buckets(max_degree+1);
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> degree_sequence;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::out_degree(*vIt, G);
        buckets[degree].push_back(*vIt);
    }
    for(unsigned int i = 0; i <= max_degree; i++){
        for(unsigned int j = 0; j < buckets[i].size(); j++){
            degree_sequence.push_back(buckets[i][j]);
        }
    }

    for(unsigned int i = 0; i < degree_sequence.size(); i++){
        for(unsigned int j = 0; j < i; j++){
            std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(degree_sequence[i], degree_sequence[j], G);
            if(existsEdge.second)
                continue;
            
            //gammaD-left heuristic
            unsigned int degree = boost::out_degree(degree_sequence[i], G);
            for(unsigned int k = 0; k < i; k++){
                boost::clear_vertex(degree_sequence[k], G);
                boost::remove_vertex(degree_sequence[k], G);
            }
            lb = (degree > lb)? degree : lb;
            
            _gammaD_left(G, lb);
            return;
        }
    }
}

template <typename G_t>
int gammaD_left(G_t &G){
    if(boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return -1;
    }
    G_t H;
    TD_copy_graph(G, H);

    unsigned int lb = 0;
    _gammaD_left(H, lb);
    return lb;
}

template <typename G_t>
void _gammaD_right(G_t &G, unsigned int &lb){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    
    //sort the vertices of G according to rising degree -> degree_sequence
    unsigned int max_degree = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::out_degree(*vIt, G);
        max_degree = (degree>max_degree)? degree : max_degree;
    }
    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > buckets(max_degree+1);
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> degree_sequence;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::out_degree(*vIt, G);
        buckets[degree].push_back(*vIt);
    }
    for(unsigned int i = 0; i <= max_degree; i++){
        for(unsigned int j = 0; j < buckets[i].size(); j++){
            degree_sequence.push_back(buckets[i][j]);
        }
    }

    for(unsigned int i = 0; i < degree_sequence.size(); i++){
        for(unsigned int j = 0; j < i; j++){
            std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(degree_sequence[i], degree_sequence[j], G);
            if(existsEdge.second)
                continue;
            
            //gammaD-right heuristic
            unsigned int degree = boost::out_degree(degree_sequence[i], G);
            boost::clear_vertex(degree_sequence[i], G);
            boost::remove_vertex(degree_sequence[i], G);
            
            lb = (degree > lb)? degree : lb;
            
            _gammaD_right(G, lb);
            return;
        }
    }
}

template <typename G_t>
int gammaD_right(G_t &G){
    if(boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return -1;
    }
    G_t H;
    TD_copy_graph(G, H);

    unsigned int lb = 0;
    _gammaD_right(H, lb);
    return lb;
}

template <typename G_t>
void _gammaD_min_e(G_t &G, unsigned int &lb){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    
    //sort the vertices of G according to rising degree -> degree_sequence
    unsigned int max_degree = 0;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::out_degree(*vIt, G);
        max_degree = (degree>max_degree)? degree : max_degree;
    }
    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > buckets(max_degree+1);
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> degree_sequence;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        unsigned int degree = boost::out_degree(*vIt, G);
        buckets[degree].push_back(*vIt);
    }
    for(unsigned int i = 0; i <= max_degree; i++){
        for(unsigned int j = 0; j < buckets[i].size(); j++){
            degree_sequence.push_back(buckets[i][j]);
        }
    }

    for(unsigned int i = 0; i < degree_sequence.size(); i++){
        for(unsigned int j = 0; j < i; j++){
            std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(degree_sequence[i], degree_sequence[j], G);
            if(existsEdge.second)
                continue;
            
            //gammaD-min-e heuristic
            unsigned int degree_right = boost::out_degree(degree_sequence[i], G);
            unsigned int degree_left = 0;
            for(unsigned int k = 0; k < i; k++)
                degree_left += boost::out_degree(degree_sequence[k], G);
            
            if(degree_left < degree_right){
                for(unsigned int t = 0; t < i; t++){
                    boost::clear_vertex(degree_sequence[t], G);
                    boost::remove_vertex(degree_sequence[t], G);
                }
            }
            else{
                boost::clear_vertex(degree_sequence[i], G);
                boost::remove_vertex(degree_sequence[i], G);
            }
            
            lb = (degree_right > lb)? degree_right : lb;
            
            _gammaD_min_e(G, lb);
            return;
        }
    }
}

template <typename G_t>
int gammaD_min_e(G_t &G){
    if(boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return -1;
    }
    G_t H;
    TD_copy_graph(G, H);

    unsigned int lb = 0;
    _gammaD_min_e(H, lb);
    return lb;
}

template <typename G_t>
int _deltaC_min_d(G_t &G){
    if(boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return -1;
    }

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;

    unsigned int lb = 0;

    while(boost::num_edges(G) > 0){
        //search a minimum-degree-vertex
        unsigned int min_degree_v = boost::num_vertices(G);

        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int degree = boost::out_degree(*vIt, G);
            if(degree <= min_degree_v){
                min_degree_v = degree;
                min_vertex = *vIt;
            }
        }
        if(min_degree_v == 0){
            boost::remove_vertex(min_vertex, G);
            continue;
        }
        
        lb = (lb>min_degree_v)? lb : min_degree_v;

        //min_d heuristic: search the neighbour of min_vertex with minimal degree
        typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nEnd1, nIt2, nEnd2;
        unsigned int min_degree_w = boost::num_vertices(G);
        typename boost::graph_traits<G_t>::vertex_descriptor w;
        for(boost::tie(nIt1, nEnd1) = boost::adjacent_vertices(min_vertex, G); nIt1 != nEnd1; nIt1++){
            unsigned int degree = boost::out_degree(*nIt1, G);
            if(degree <= min_degree_w){
                min_degree_w = degree;
                w = *nIt1;
            }
        }
        
        //contract the edge between min_vertex and w
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N;
        for(boost::tie(nIt1, nEnd1) = boost::adjacent_vertices(min_vertex, G); nIt1 != nEnd1; nIt1++)
            N.insert(*nIt1);
        for(boost::tie(nIt2, nEnd2) = boost::adjacent_vertices(w, G); nIt2 != nEnd2; nIt2++)
            N.insert(*nIt2);
        
        typename boost::graph_traits<G_t>::vertex_descriptor new_v = boost::add_vertex(G);
        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = N.begin(); sIt != N.end(); sIt++)
            boost::add_edge(new_v, *sIt, G);
        
        boost::clear_vertex(min_vertex, G);
        boost::clear_vertex(w, G);
        boost::remove_vertex(min_vertex, G);
        boost::remove_vertex(w, G);
    }
    return (int)lb;
}

template <typename G_t>
int deltaC_min_d(G_t &G){
    G_t H;
    TD_copy_graph(G, H);
    return _deltaC_min_d(H);
}

template <typename G_t>
int _deltaC_max_d(G_t &G){
    if(boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return -1;
    }

    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;

    unsigned int lb = 0;

    while(boost::num_edges(G) > 0){
        //search a minimum-degree-vertex
        unsigned int min_degree = boost::num_vertices(G);

        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int degree = boost::out_degree(*vIt, G);
            if(degree < min_degree){
                min_degree = degree;
                min_vertex = *vIt;
            }
        }
        
        lb = (lb>min_degree)? lb : min_degree;

        if(min_degree == 0){
            boost::remove_vertex(min_vertex, G);
            continue;
        }

        //max_d heuristic: search the neighbour of min_vertex with maximal degree
        typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nEnd1, nIt2, nEnd2;
        unsigned int max_degree = 0;
        typename boost::graph_traits<G_t>::vertex_descriptor w;
        for(boost::tie(nIt1, nEnd1) = boost::adjacent_vertices(min_vertex, G); nIt1 != nEnd1; nIt1++){
            unsigned int degree = boost::out_degree(*nIt1, G);
            if(degree > max_degree){
                max_degree = degree;
                w = *nIt1;
            }
        }

        //contract the edge between min_vertex and w
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N;
        for(boost::tie(nIt1, nEnd1) = boost::adjacent_vertices(min_vertex, G); nIt1 != nEnd1; nIt1++)
            N.insert(*nIt1);
        for(boost::tie(nIt2, nEnd2) = boost::adjacent_vertices(w, G); nIt2 != nEnd2; nIt2++)
            N.insert(*nIt2);
        
        typename boost::graph_traits<G_t>::vertex_descriptor new_v = boost::add_vertex(G);
        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = N.begin(); sIt != N.end(); sIt++)
            boost::add_edge(new_v, *sIt, G);
        
        boost::clear_vertex(min_vertex, G);
        boost::clear_vertex(w, G);
        boost::remove_vertex(min_vertex, G);
        boost::remove_vertex(w, G);
    }
    return (int)lb;
}

template <typename G_t>
int deltaC_max_d(G_t &G){
    G_t H;
    TD_copy_graph(G, H);
    return _deltaC_max_d(H);
}

template <typename G_t>
int _deltaC_least_c(G_t &G){
    //checking whether G is a complete graph
    if(boost::num_edges(G) == boost::num_vertices(G)*(boost::num_vertices(G)-1)){
        return boost::num_vertices(G)-1;
    }
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;

    unsigned int lb = 0;

    while(boost::num_edges(G) > 0){
        //search a minimum-degree-vertex
        unsigned int min_degree = boost::num_vertices(G);

        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            unsigned int degree = boost::out_degree(*vIt, G);
            if(degree > 0 && degree < min_degree){
                min_degree = degree;
                min_vertex = *vIt;
            }
        }
        
        lb = (lb>min_degree)? lb : min_degree;

        if(min_degree == boost::num_vertices(G))
            return (int)lb;

        //least-c heuristic: search the neighbour of min_vertex such that contracting {min_vertex, w} removes least edges
        typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
        typename boost::graph_traits<G_t>::vertex_descriptor w;
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N;
        
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, G); nIt != nEnd; nIt++)
            N.insert(*nIt);
        
        unsigned int cnt_common, min_common = N.size();
        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt1, sIt2;
        
        for(sIt1 = N.begin(); sIt1 != N.end(); sIt1++){
            cnt_common = 0;
            sIt2 = sIt1;
            sIt2++;
            for(; sIt2 != N.end(); sIt2++){
                std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(*sIt1, *sIt2, G);
                if(existsEdge.second)
                    cnt_common += 1;
            }
            if(cnt_common < min_common){
                w = *sIt1;
                min_common = cnt_common;
            }
        }

        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(w, G); nIt != nEnd; nIt++)
            N.insert(*nIt);

        //contract the edge between min_vertex and w
        typename boost::graph_traits<G_t>::vertex_descriptor new_v = boost::add_vertex(G);

        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = N.begin(); sIt != N.end(); sIt++)
            boost::add_edge(new_v, *sIt, G);

        
        boost::clear_vertex(min_vertex, G);
        boost::clear_vertex(w, G);
    }
    return (int)lb;
}

template <typename G_t>
int deltaC_least_c(G_t &G){
    G_t H;
    TD_copy_graph(G, H);
    return _deltaC_least_c(H);
}


/* IMPROVED GRAPHS */

template <typename G_t>
void k_neighbour_improved_graph(G_t &G, unsigned int k){
    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > not_edges;
    typename boost::graph_traits<G_t>::vertex_iterator vIt1, vIt2, vEnd;
    
    for(boost::tie(vIt1, vEnd) = boost::vertices(G); vIt1 != vEnd; vIt1++){
        vIt2 = vIt1;
        vIt2++;
        for(; vIt2 != vEnd; vIt2++){
            std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(*vIt1, *vIt2, G);
            if(!existsEdge.second){
                std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> not_edge;
                not_edge.push_back(*vIt1);
                not_edge.push_back(*vIt2);
                not_edges.push_back(not_edge);
            }
        }
    }
    
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    
    for(unsigned int i = 0; i < not_edges.size(); i++){
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N1, N2;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(not_edges[i][0], G); nIt != nEnd; nIt++)
            N1.insert(*nIt);
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(not_edges[i][1], G); nIt != nEnd; nIt++)
            N2.insert(*nIt);
        
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> intersection;

        std::set_intersection(N1.begin(), N1.end(), N2.begin(), N2.end(), std::inserter(intersection, intersection.begin()));
        
        if(intersection.size() >= k)
            boost::add_edge(not_edges[i][0], not_edges[i][1], G);
    }
}

template <typename G_t>
int LBN_deltaD(G_t &G){
    int lb = deltaD(G);
    
    while(true){
        G_t H;
        TD_copy_graph(G, H);
        k_neighbour_improved_graph(H, lb+1);

        if(deltaD(H) > lb)
            lb++;
        else
            break;   
    }
    return lb;
}

template <typename G_t>
int LBN_deltaC(G_t &G){
    int lb = deltaC_least_c(G);
    
    while(true){
        G_t H;
        TD_copy_graph(G, H);
        k_neighbour_improved_graph(H, lb+1);
        
        if(deltaC_least_c(H) > lb)
            lb++;
        else
            break;   
    }
    return lb;
}


template <typename G_t>
int LBNC_deltaD(G_t &G){
    int lb = deltaD(G);

    while(true){
        G_t H;
        TD_copy_graph(G, H);
        k_neighbour_improved_graph(H, lb+1);

        int new_lb;

        while(boost::num_edges(H) > 0){
            new_lb = deltaD(H);
            if(new_lb > lb)
                break;
            
            //compute a min-degree-vertex
            typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
            typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;
            unsigned int min_degree = boost::num_vertices(H);
            for(boost::tie(vIt, vEnd) = boost::vertices(H); vIt != vEnd; vIt++){
                unsigned int degree = boost::out_degree(*vIt, H);
                if(degree < min_degree){
                    min_degree = degree;
                    min_vertex = *vIt;
                }
            }
            if(min_degree == 0){
                 boost::remove_vertex(min_vertex, H);
                 continue;
            }
            
            //compute a neighbour of min_vertex, such that they share a minimal count of vertices
            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            typename boost::graph_traits<G_t>::vertex_descriptor w;
            std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N;
                
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, H); nIt != nEnd; nIt++)
                N.insert(*nIt);
                
            unsigned int cnt_common, min_common = N.size();
            typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt1, sIt2;
                
            for(sIt1 = N.begin(); sIt1 != N.end(); sIt1++){
                cnt_common = 0;
                sIt2 = sIt1;
                sIt2++;
                for(; sIt2 != N.end(); sIt2++){
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(*sIt1, *sIt2, H);
                    if(existsEdge.second)
                        cnt_common += 1;
                }
                if(cnt_common < min_common){
                    w = *sIt1;
                    min_common = cnt_common;
                }
            }
                
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(w, H); nIt != nEnd; nIt++)
                N.insert(*nIt);

            //contract the edge between min_vertex and w
            typename boost::graph_traits<G_t>::vertex_descriptor new_v = boost::add_vertex(H);
            H[new_v].id = H[min_vertex].id;
            
            for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = N.begin(); sIt != N.end(); sIt++)
                boost::add_edge(new_v, *sIt, H);
                
            boost::clear_vertex(min_vertex, H);
            boost::clear_vertex(w, H);
            boost::remove_vertex(min_vertex, H);
            boost::remove_vertex(w, H);

            k_neighbour_improved_graph(H, lb+1);
        }
        if(new_lb > lb)
            lb++;
        else
            break;
    }
    return lb;
}

template <typename G_t>
int LBNC_deltaC(G_t &G){
    int lb = deltaC_least_c(G);

    while(true){
        G_t H;
        TD_copy_graph(G, H);
        k_neighbour_improved_graph(H, lb+1);

        int new_lb;

        while(boost::num_edges(H) > 0){
            new_lb = deltaC_least_c(H);
            if(new_lb > lb)
                break;
            
            //compute a min-degree-vertex
            typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
            typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;
            unsigned int min_degree = boost::num_vertices(H);
            for(boost::tie(vIt, vEnd) = boost::vertices(H); vIt != vEnd; vIt++){
                unsigned int degree = boost::out_degree(*vIt, H);
                if(degree < min_degree){
                    min_degree = degree;
                    min_vertex = *vIt;
                }
            }
            if(min_degree == 0){
                 boost::remove_vertex(min_vertex, H);
                 continue;
            }
            
            //compute a neighbour of min_vertex, such that they share a minimal count of vertices
            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            typename boost::graph_traits<G_t>::vertex_descriptor w;
            std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N;
                
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, H); nIt != nEnd; nIt++)
                N.insert(*nIt);
                
            unsigned int cnt_common, min_common = N.size();
            typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt1, sIt2;
                
            for(sIt1 = N.begin(); sIt1 != N.end(); sIt1++){
                cnt_common = 0;
                sIt2 = sIt1;
                sIt2++;
                for(; sIt2 != N.end(); sIt2++){
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(*sIt1, *sIt2, H);
                    if(existsEdge.second)
                        cnt_common += 1;
                }
                if(cnt_common < min_common){
                    w = *sIt1;
                    min_common = cnt_common;
                }
            }
                
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(w, H); nIt != nEnd; nIt++)
                N.insert(*nIt);

            //contract the edge between min_vertex and w
            typename boost::graph_traits<G_t>::vertex_descriptor new_v = boost::add_vertex(H);
            H[new_v].id = H[min_vertex].id;
            
            for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = N.begin(); sIt != N.end(); sIt++)
                boost::add_edge(new_v, *sIt, H);

            boost::clear_vertex(min_vertex, H);
            boost::clear_vertex(w, H);
            boost::remove_vertex(min_vertex, H);
            boost::remove_vertex(w, H);

            k_neighbour_improved_graph(H, lb+1);
        }
        if(new_lb > lb)
            lb++;
        else
            break;
    }
    return lb;
}


template <typename G_t>
void k_path_improved_graph(G_t &G, unsigned int k){
    std::vector<std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> > not_edges;
    
    typename boost::graph_traits<G_t>::vertex_iterator vIt1, vIt2, vEnd;
    for(boost::tie(vIt1, vEnd) = boost::vertices(G); vIt1 != vEnd; vIt1++){
        vIt2 = vIt1;
        vIt2++;
        for(; vIt2 != vEnd; vIt2++){
            std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(*vIt1, *vIt2, G);
            if(!existsEdge.second){
                std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> not_edge;
                not_edge.push_back(*vIt1);
                not_edge.push_back(*vIt2);
                not_edges.push_back(not_edge);
            }
        }
    }
    
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    G_t I;
    TD_copy_graph(G, I);

    for(unsigned int i = 0; i < not_edges.size(); i++){
        std::set<unsigned int> X, Y, S, D;

        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(not_edges[i][0], I); nIt != nEnd; nIt++)
            X.insert(G[*nIt].id);

        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(not_edges[i][1], I); nIt != nEnd; nIt++)
            Y.insert(G[*nIt].id);

        D.insert(I[not_edges[i][0]].id);
        D.insert(I[not_edges[i][1]].id);

        G_t H = graph_after_deletion(I, D);

        seperate_vertices(H, X, Y, S);
        
        
        if(S.size() >= k)
            boost::add_edge(not_edges[i][0], not_edges[i][1], G);
    }
}

template <typename G_t>
int LBP_deltaD(G_t &G){
    int lb = deltaD(G);

    while(true){
        G_t H;
        TD_copy_graph(G, H);
        k_path_improved_graph(H, lb+1);
        
        if (deltaD(H) > lb)
            lb++;
        else
            break;   
    }
    return lb;
}

template <typename G_t>
int LBP_deltaC(G_t &G){
    int lb = deltaC_least_c(G);
    
    while(true){
        G_t H;
        TD_copy_graph(G, H);
        k_path_improved_graph(H, lb+1);

        if (deltaC_least_c(H) > lb)
            lb++;
        else
            break;   
    }
    return lb;
}

template <typename G_t>
int LBPC_deltaD(G_t &G){
    int lb = deltaD(G);

    while(true){
        G_t H;
        TD_copy_graph(G, H);
        k_path_improved_graph(H, lb+1);

        int new_lb;
        
        while(boost::num_edges(H) > 0){
            new_lb = deltaD(H);
            if(new_lb > lb)
                break;
            
            //compute a min-degree-vertex
            typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
            typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;
            unsigned int min_degree = boost::num_vertices(H);
            for(boost::tie(vIt, vEnd) = boost::vertices(H); vIt != vEnd; vIt++){
                unsigned int degree = boost::out_degree(*vIt, H);
                if(degree < min_degree){
                    min_degree = degree;
                    min_vertex = *vIt;
                }
            }
            if(min_degree == 0){
                 boost::remove_vertex(min_vertex, H);
                 continue;
            }
            
            //compute a neighbour of min_vertex, such that they share a minimal count of vertices
            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            typename boost::graph_traits<G_t>::vertex_descriptor w;
            std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N;
                
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, H); nIt != nEnd; nIt++)
                N.insert(*nIt);
                
            unsigned int cnt_common, min_common = N.size();
            typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt1, sIt2;
                
            for(sIt1 = N.begin(); sIt1 != N.end(); sIt1++){
                cnt_common = 0;
                sIt2 = sIt1;
                sIt2++;
                for(; sIt2 != N.end(); sIt2++){
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(*sIt1, *sIt2, H);
                    if(existsEdge.second)
                        cnt_common += 1;
                }
                if(cnt_common < min_common){
                    w = *sIt1;
                    min_common = cnt_common;
                }
            }
                
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(w, H); nIt != nEnd; nIt++)
                N.insert(*nIt);

            //contract the edge between min_vertex and w
            typename boost::graph_traits<G_t>::vertex_descriptor new_v = boost::add_vertex(H);
            H[new_v].id = H[min_vertex].id;
            
            for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = N.begin(); sIt != N.end(); sIt++)
                boost::add_edge(new_v, *sIt, H);
                
            boost::clear_vertex(min_vertex, H);
            boost::clear_vertex(w, H);
            boost::remove_vertex(min_vertex, H);
            boost::remove_vertex(w, H);

            k_path_improved_graph(H, lb+1);

        }
        if(new_lb > lb)
            lb++;
        else
            break;
    }
    return lb;
}

template <typename G_t>
int LBPC_deltaC(G_t &G){
    int lb = deltaC_least_c(G);

    while(true){
        G_t H;
        TD_copy_graph(G, H);
        k_path_improved_graph(H, lb+1);

        int new_lb;

        while(boost::num_edges(H) > 0){
            new_lb = deltaC_least_c(H);

            if(new_lb > lb)
                break;
            
            //compute a min-degree-vertex
            typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
            typename boost::graph_traits<G_t>::vertex_descriptor min_vertex;
            unsigned int min_degree = boost::num_vertices(H);
            for(boost::tie(vIt, vEnd) = boost::vertices(H); vIt != vEnd; vIt++){
                unsigned int degree = boost::out_degree(*vIt, H);
                if(degree < min_degree){
                    min_degree = degree;
                    min_vertex = *vIt;
                }
            }
            if(min_degree == 0){
                 boost::remove_vertex(min_vertex, H);
                 continue;
            }
            
            //compute a neighbour of min_vertex, such that they share a minimal count of vertices
            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            typename boost::graph_traits<G_t>::vertex_descriptor w;
            std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N;
                
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(min_vertex, H); nIt != nEnd; nIt++)
                N.insert(*nIt);
                
            unsigned int cnt_common, min_common = N.size();
            typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt1, sIt2;
                
            for(sIt1 = N.begin(); sIt1 != N.end(); sIt1++){
                cnt_common = 0;
                sIt2 = sIt1;
                sIt2++;
                for(; sIt2 != N.end(); sIt2++){
                    std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(*sIt1, *sIt2, H);
                    if(existsEdge.second)
                        cnt_common += 1;
                }
                if(cnt_common < min_common){
                    w = *sIt1;
                    min_common = cnt_common;
                }
            }
                
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(w, H); nIt != nEnd; nIt++)
                N.insert(*nIt);

            //contract the edge between min_vertex and w
            typename boost::graph_traits<G_t>::vertex_descriptor new_v = boost::add_vertex(H);
            H[new_v].id = H[min_vertex].id;
            
            for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = N.begin(); sIt != N.end(); sIt++)
                boost::add_edge(new_v, *sIt, H);
                
            boost::clear_vertex(min_vertex, H);
            boost::clear_vertex(w, H);
            boost::remove_vertex(min_vertex, H);
            boost::remove_vertex(w, H);

            k_path_improved_graph(H, lb+1);
        }
        if(new_lb > lb)
            lb++;
        else
            break;
    }
    return lb;
}




/* Maximum Cardinality Search-based */

//does a maximum cardinality search and returns the vertex descriptors of the maximum visited degree-vertex and the last visited vertex
template <typename G_t>
std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> MCS_single(G_t &G, int &lb, typename boost::graph_traits<G_t>::vertex_descriptor start_vertex){
    std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N, visited;
    visited.insert(start_vertex);
    int current_visited_degree, max_visited_degree = 0;
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    typename boost::graph_traits<G_t>::vertex_descriptor current_visited_vertex, max_visited_vertex;
    current_visited_vertex = start_vertex;
 
    for(unsigned int i = 1; i < boost::num_vertices(G); i++){
        current_visited_degree = 0;

        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(current_visited_vertex, G); nIt != nEnd; nIt++){
            if(visited.find(*nIt) == visited.end())
                N.insert(*nIt);
        }
        
        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = N.begin(); sIt != N.end(); sIt++){
            std::set<typename boost::graph_traits<G_t>::vertex_descriptor> n;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*sIt, G); nIt != nEnd; nIt++){
                n.insert(*nIt);
            }
            
            std::set<typename boost::graph_traits<G_t>::vertex_descriptor> intersection;
            std::set_intersection(visited.begin(), visited.end(), n.begin(), n.end(), std::inserter(intersection, intersection.begin()));
            
            if((int)intersection.size() > current_visited_degree){
                current_visited_degree = (int)intersection.size();
                current_visited_vertex = *sIt;
            }
        }
        if(current_visited_degree > max_visited_degree){
            max_visited_degree = current_visited_degree;
            max_visited_vertex = current_visited_vertex;
        }    

        visited.insert(current_visited_vertex);
        N.erase(current_visited_vertex);
    }
    lb = (lb > max_visited_degree) ? lb : max_visited_degree;
    
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> rtn;
    rtn.push_back(max_visited_vertex);
    rtn.push_back(current_visited_vertex);
    return rtn;
}

template <typename G_t>
std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> MCS_all_start_vertices(G_t &G, int &lb){
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    int last_lb = lb;
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> data, rtn;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
        data = MCS_single(G, lb, *vIt);
        if(last_lb < lb)
            rtn = data;
        last_lb = lb;
    }
    return rtn;
}

template <typename G_t>
void MCS_random(G_t &G, int &lb){
    std::srand(time(NULL));
    
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    boost::tie(vIt, vEnd) = boost::vertices(G);
    std::advance(vIt, rand() % boost::num_vertices(G)); 
    MCS_single(G, lb, *vIt);
}

template <typename G_t>
void MCSC_min_deg(G_t H, int &lb){
    typename boost::graph_traits<G_t>::vertex_descriptor v, w;
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> data;
    data = MCS_all_start_vertices(H, lb);
    v = data[0];
    
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    while(boost::num_vertices(H) > 2){        
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, H); nIt != nEnd; nIt++)
            N.insert(*nIt);
        
        int cnt_common, min_common = (int)N.size();
        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt1, sIt2;
        
        for(sIt1 = N.begin(); sIt1 != N.end(); sIt1++){
            cnt_common = 0;
            sIt2 = sIt1;
            sIt2++;
            for(; sIt2 != N.end(); sIt2++){
                std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(*sIt1, *sIt2, H);
                if(existsEdge.second)
                    cnt_common += 1;
            }
            if(cnt_common < min_common){
                w = *sIt1;
                min_common = cnt_common;
            }
        }
        
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(w, H); nIt != nEnd; nIt++)
            N.insert(*nIt);
        
        typename boost::graph_traits<G_t>::vertex_descriptor new_v = boost::add_vertex(H);
        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = N.begin(); sIt != N.end(); sIt++)
            boost::add_edge(new_v, *sIt, H);
        
        boost::clear_vertex(v, H);
        boost::clear_vertex(w, H);
        boost::remove_vertex(v, H);
        boost::remove_vertex(w, H);
        
        int deg, min_deg = (int)boost::num_vertices(H);
        for(boost::tie(vIt, vEnd) = boost::vertices(H); vIt != vEnd; vIt++){
            deg = boost::out_degree(*vIt, H);
            if(deg < min_deg){
                min_deg = deg;
                v = *vIt;
            }
        }
    }
}

template <typename G_t>
void MCSC_last_mcs(G_t H, int &lb){
    typename boost::graph_traits<G_t>::vertex_descriptor v, w;
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> data;
    data = MCS_all_start_vertices(H, lb);
    v = data[0];
    
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    while(boost::num_vertices(H) > 2){
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, H); nIt != nEnd; nIt++)
            N.insert(*nIt);
        
        int cnt_common, min_common = (int)N.size();
        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt1, sIt2;
        
        for(sIt1 = N.begin(); sIt1 != N.end(); sIt1++){
            cnt_common = 0;
            sIt2 = sIt1;
            sIt2++;
            for(; sIt2 != N.end(); sIt2++){
                std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(*sIt1, *sIt2, H);
                if(existsEdge.second)
                    cnt_common += 1;
            }
            if(cnt_common < min_common){
                w = *sIt1;
                min_common = cnt_common;
            }
        }
        
        for(boost::tie(nIt, nEnd) = adjacent_vertices(w, H); nIt != nEnd; nIt++)
            N.insert(*nIt);
        
        typename boost::graph_traits<G_t>::vertex_descriptor new_v = boost::add_vertex(H);
        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = N.begin(); sIt != N.end(); sIt++)
            boost::add_edge(new_v, *sIt, H);
        
        boost::clear_vertex(v, H);
        boost::clear_vertex(w, H);
        boost::remove_vertex(v, H);
        boost::remove_vertex(w, H);
        
        data = MCS_single(H, lb, new_v);
        v = data[1];
    }
}

template <typename G_t>
void MCSC_max_mcs(G_t H, int &lb){
    typename boost::graph_traits<G_t>::vertex_descriptor v, w;
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> data;
    data = MCS_all_start_vertices(H, lb);
    v = data[0];
    
    typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    while(boost::num_vertices(H) > 2){
        std::set<typename boost::graph_traits<G_t>::vertex_descriptor> N;
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, H); nIt != nEnd; nIt++)
            N.insert(*nIt);
        
        int cnt_common, min_common = (int)N.size();
        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt1, sIt2;
        
        for(sIt1 = N.begin(); sIt1 != N.end(); sIt1++){
            cnt_common = 0;
            sIt2 = sIt1;
            sIt2++;
            for(; sIt2 != N.end(); sIt2++){
                std::pair<typename G_t::edge_descriptor, bool> existsEdge = boost::edge(*sIt1, *sIt2, H);
                if(existsEdge.second)
                    cnt_common += 1;
            }
            if(cnt_common < min_common){
                w = *sIt1;
                min_common = cnt_common;
            }
        }
        
        for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(w, H); nIt != nEnd; nIt++)
            N.insert(*nIt);
        
        typename boost::graph_traits<G_t>::vertex_descriptor new_v = boost::add_vertex(H);
        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = N.begin(); sIt != N.end(); sIt++)
            boost::add_edge(new_v, *sIt, H);
        
        boost::clear_vertex(v, H);
        boost::clear_vertex(w, H);
        boost::remove_vertex(v, H);
        boost::remove_vertex(w, H);
        
        data = MCS_single(H, lb, new_v);
        v = data[0];
    }
}

template <typename G_t>
int relation_edges_vertices(G_t &G){
    if(boost::num_vertices(G) == 0)
        return -1;
    return (int)(2*boost::num_edges(G)/boost::num_vertices(G));
}

}
}

#endif
