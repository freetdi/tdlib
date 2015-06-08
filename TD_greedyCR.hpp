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
// Offers functionality to compute a tree decomposition of exact width. 
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
//
// These functions are most likely to be interesting for outside use:
//
// void CR_greedy_decomp(G_t &G, T_t &T, int lb)
// void CR_greedy_decomp(G_t &G, T_t &T)
//

#include <boost/graph/adjacency_list.hpp>
#include "simple_graph_algos.hpp"

namespace treedec{
    
#ifndef TD_SUBSETS
#define TD_SUBSETS

void subsets(std::set<unsigned int> &X, int size, int k, unsigned int idx, std::vector<unsigned int> &sub, std::vector<std::set<unsigned int> > &subs){
    if(k==0){
        std::set<unsigned int> subS;
	std::copy(sub.begin(), sub.end(), std::inserter(subS, subS.end()));
        subs.push_back(subS);
        return;
    }
    
    unsigned int i = idx;
    std::set<unsigned int>::iterator sIt = X.begin();
    std::advance(sIt, i);
    for(; i<X.size();i++){
        sub.push_back(*sIt);
        subsets(X,X.size(),k-1,i+1,sub, subs);
        sub.pop_back();
        sIt++;
    }
} 

#endif

#ifndef FUNC_NEW_ROBBER_SPACE
#define FUNC_NEW_ROBBER_SPACE

//computes the new robber space with respect to X and R and saves it in newR
template <typename G_t>
void new_robber_space(G_t G, std::set<unsigned int> &X, std::set<unsigned int> &R, std::set<unsigned int> &newR){
    //G \ X
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;

    for(std::set<unsigned int>::iterator sIt = X.begin(); sIt != X.end(); sIt++){
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            if (G[*vIt].id == *sIt){
                boost::clear_vertex(*vIt, G);
                boost::remove_vertex(*vIt, G);
                break;
            }
        }
    }

    //compute new robber space
    std::vector<std::set<unsigned int> > components;    
    get_components(G, components);   
    
    for(unsigned int i = 0; i < components.size(); i++){
        for(std::set<unsigned int>::iterator sIt1 = components[i].begin(); sIt1 != components[i].end(); sIt1++){
            if(R.find(*sIt1) != R.end()){
                for(std::set<unsigned int>::iterator sIt2 = components[i].begin(); sIt2 != components[i].end(); sIt2++){
                    newR.insert(*sIt2);
                }
                break;
            }
        }
    }
}

#endif

#ifndef FUNC_IS_MONOTONE
#define FUNC_IS_MONOTONE

//checks if robber_space(G\(X ^ X')) == robber_space(G\X)
template <typename G_t>
bool is_monotone(G_t G, std::set<unsigned int> &X, std::set<unsigned int> &X_prime, std::set<unsigned int> &R, std::set<unsigned int> &newR){ 
    if(X == X_prime)
        return false;

    //robber_space(G\(X ^ X'))
    std::set<unsigned int> intersection;  
    std::set_intersection(X.begin(), X.end(), X_prime.begin(), X_prime.end(), std::inserter(intersection, intersection.begin()));   
    std::set<unsigned int> newR1;
    new_robber_space(G, intersection, R, newR1);
    
    //robber_space(G\X)    
    std::set<unsigned int> newR2;
    new_robber_space(G, X, R, newR2);
    
    if(newR1 == newR2){
        newR = newR1;
        for(std::set<unsigned int>::iterator sIt = X_prime.begin(); sIt != X_prime.end(); sIt++)
            newR.erase(*sIt);
              
        return true;
    }
    return false;
}

#endif

#ifndef CR_GLUE_BAGS
#define CR_GLUE_BAGS

template <typename T_t>
void CR_glue_bags(T_t &T, std::set<unsigned int> bag1, std::set<unsigned int> &bag2){
    typename boost::graph_traits<T_t>::vertex_iterator vIt, vEnd;
    
    for(boost::tie(vIt, vEnd) = boost::vertices(T); vIt != vEnd; vIt++){
        if(T[*vIt].bag == bag1){
            typename boost::graph_traits<T_t>::vertex_descriptor b = boost::add_vertex(T);
            T[b].bag = bag2;
            boost::add_edge(*vIt, b, T);
            return;
        }
    }
    
    typename boost::graph_traits<T_t>::vertex_descriptor b1 = boost::add_vertex(T);
    T[b1].bag = bag1;
    typename boost::graph_traits<T_t>::vertex_descriptor b2 = boost::add_vertex(T);
    T[b2].bag = bag2;
    boost::add_edge(b1, b2, T);
}

#endif

#ifndef CR_NEW_ROBBER_COMPONENTS
#define CR_NEW_ROBBER_COMPONENTS

template<typename G_t>
void new_robber_components(G_t G, std::set<unsigned int> &newR, std::vector<std::set<unsigned int> > &newRcomps){
    std::set<unsigned int> cR;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        cR.insert(G[*vIt].id);
    for(std::set<unsigned int>::iterator sIt = newR.begin(); sIt != newR.end(); sIt++)
        cR.erase(*sIt);
        
    for(unsigned int t = 0; t < cR.size(); t++){
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            if(cR.find(G[*vIt].id) !=  cR.end()){
                boost::clear_vertex(*vIt, G);
                boost::remove_vertex(*vIt, G);
                break;
            }
        }
    }

    std::vector<std::set<unsigned int> > components;
    get_components(G, components);  
    
    newRcomps = components;
}

#endif


template <typename G_t>
bool hunt(G_t &G, std::set<unsigned int> &X, std::set<unsigned int> &R, unsigned int k, std::vector<std::vector<std::set<unsigned int> > > &monotone_turns){
    std::set<unsigned int> U = X;
    for(std::set<unsigned int>::iterator sIt = R.begin(); sIt != R.end(); sIt++)
        U.insert(*sIt);
    
    std::vector<unsigned int> sub;
    std::vector<std::set<unsigned int> > subs;
    subsets(U, U.size(), k, 0, sub, subs);
     
    for(unsigned int i = 0; i < subs.size(); i++){  
	std::set<unsigned int> newR;
        if(!is_monotone(G, X, subs[i], R, newR))
            continue;
        
        
        //new robber space components
        std::vector<std::set<unsigned int> > newRcomps;
        new_robber_components(G, newR, newRcomps);
        
        unsigned int idx = monotone_turns.size();

        bool all_monotone = true;
        for(unsigned int j = 0; j < newRcomps.size(); j++){
            if(newRcomps[j].size() == 0)
                continue;

            all_monotone = all_monotone && hunt(G, subs[i], newRcomps[j], k, monotone_turns);
            if(!all_monotone)
                break;
        }
        if(!all_monotone){
            for(unsigned int l = 0; l < monotone_turns.size()-idx; l++)
                monotone_turns.pop_back();
            continue;
        }

	std::vector<std::set<unsigned int> > monotone_turn;
        monotone_turn.push_back(X);
        monotone_turn.push_back(subs[i]);
        monotone_turns.push_back(monotone_turn);

        return true;
    }
    return false;
}

template <typename G_t, typename T_t>
void CR_greedy_decomp(G_t &G, T_t &T, int lb){
    if(boost::num_vertices(G) == 0){
        boost::add_vertex(T);
        return;
    }
    
    std::set<unsigned int> vertices;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        vertices.insert(G[*vIt].id);

    std::vector<std::vector<std::set<unsigned int> > > monotone_turns;
    unsigned int k = (unsigned int)++lb;
    
    while(true){
        std::set<unsigned int> emptySet;
        if(hunt(G, emptySet, vertices, k, monotone_turns))
            break;
        monotone_turns.clear();
        k++;
    }   

    for(unsigned int j = monotone_turns.size(); j > 0; j--)
        CR_glue_bags(T, monotone_turns[j-1][0], monotone_turns[j-1][1]);
}

template <typename G_t, typename T_t>
void CR_greedy_decomp(G_t &G, T_t &T){
    int lb = -1;
    CR_greedy_decomp(G, T, lb);
}

}
