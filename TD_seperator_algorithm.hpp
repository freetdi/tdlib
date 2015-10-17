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
// Offers functionality to preprocess a graph, such that after
// the repeated application of reduction rules, which in case the
// input graph has tree-width at most 3 allow us to determine it's tree-width exactly 
// and in addition compute the corresponding tree decomposition. If the tree-width 
// is larger, the reduction rules return a possibly smaller instance of the same 
// tree-width as the original graph, a partial tree decomposition and a lower bound
// with respect to tree-width, such that
// further algorithms can be applied to the resulting graph.
// Currently, the minDegree-heuristic will be applied to the resulting graph,
// if the input graph can't be fully preprocessed.
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
// These function is most likely to be interesting for outside use:
//
// void seperator_algorithm(G_t &G, T_t &T)
//

#include <set>
#include "TD_NetworkFlow.hpp"
#include "TD_misc.hpp"
#include "simple_graph_algos.hpp"

namespace treedec{

//collects all subsets of X of size k in subs
template <typename T>
void subsets(std::set<T> &X, int size, int k, int idx, std::vector<T> &sub, std::vector<std::set<T> > &subs){
    if(k==0){
        typename std::set<T> subS;
        for(unsigned int i = 0; i < sub.size(); i++)
            subS.insert(sub[i]);
        subs.push_back(subS);
        return;
    }
    
    int i = idx;
    typename std::set<T>::iterator sIt = X.begin();
    std::advance(sIt, i);
    for(; i<size;i++){
        sub.push_back(*sIt);
        subsets(X,size,k-1,i+1,sub, subs);
        sub.pop_back();
        sIt++;
    }
} 

template <typename T>
void disjoint_subsets(std::set<T> &X, unsigned int min_card, unsigned int max_card, std::vector<T> &sub, std::vector<std::set<T> > &subsX, std::vector<std::vector<std::set<T> > > &subsY){
    for(unsigned int i = min_card; i <=max_card; i++){
        subsets(X, X.size(), i, 0, sub, subsX);
    }
    
    for(unsigned int i = 0; i < subsX.size(); i++){  
        std::set<T> difference; 
        std::set_difference(X.begin(), X.end(), subsX[i].begin(), subsX[i].end(), std::inserter(difference, difference.begin()));
            
        unsigned int maximum = (difference.size() > max_card)? max_card : difference.size();
        
        std::vector<std::set<T> > subsXY;
        for(unsigned int t = 1; t <= maximum; t++){
            subsets(difference, difference.size(), t, 0, sub, subsXY);
        }
        subsY.push_back(subsXY);
    }
}

template <typename T>
void superset(T &X, T &V, unsigned int size){
    typename T::iterator sIt = V.begin();
    while(X.size() != size)
        X.insert(*(sIt++));
}

//Finds a nearly balanced seperator S of W by doing an extended deepth-first-search
//which finds a minimal X-Y-seperator. The sets X and Y are all possible disjoint subsets of
//W of size 1 to 2k
template <typename G_t>
bool nearly_balanced_seperator(G_t &G, typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &W, typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &S, std::vector<bool> &disabled, unsigned int k){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> sub;
    std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > subsX;
    std::vector<std::vector<std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > > subsY;

    std::set<unsigned int> W_;
    treedec::descriptor_bag_to_id_bag(G, W_, W);

    disjoint_subsets(W, 1, 2*k, sub, subsX, subsY);

    for(unsigned int i = 0; i < subsX.size(); i++){
        for(unsigned int j = 0; j < subsY[i].size(); j++){
            //there cannot exist a X-Y-seperator if there is an edge between X and Y
            if(is_edge_between_sets(G, subsX[i], subsY[i][j]))
                continue;

            std::vector<bool> disabled_(disabled);
            std::set<unsigned int> X, Y, sX, sY, X_Y, Z;

            treedec::descriptor_bag_to_id_bag(G, X, subsX[i]);
            treedec::descriptor_bag_to_id_bag(G, Y, subsY[i][j]);
            
            std::set_union(X.begin(), X.end(), Y.begin(), Y.end(), std::inserter(X_Y, X_Y.begin()));            
            std::set_difference(W_.begin(), W_.end(), X_Y.begin(), X_Y.end(), std::inserter(Z, Z.begin()));
            
            //do the extended deepth-first-search on the neighbours of vertices in X and Y 
            get_neighbourhood_new(G, disabled_, subsX[i], sX);
            get_neighbourhood_new(G, disabled_, subsY[i][j], sY);
 
            for(std::set<unsigned int>::iterator sIt = X_Y.begin(); sIt != X_Y.end(); sIt++)
                disabled_[*sIt] = true;
                
            //Z must be a subset of S, a S-S-seperator must contain S
            sX.insert(Z.begin(), Z.end());
            sY.insert(Z.begin(), Z.end());

            std::set<unsigned int> S_;

            if(!seperate_vertices(G, disabled_, sX, sY, S_, k+1))
                continue;
            
            //S now is a sX-sY-seperator. check if S holds the remaining property of a nearly balanced seperator
            std::set<unsigned int> isS_W;
            std::set_intersection(S_.begin(), S_.end(), W_.begin(), W_.end(), std::inserter(isS_W, isS_W.begin()));
	
            if(isS_W == Z){
                treedec::id_bag_to_descriptor_bag(G, S, S_);
                return true;
            }
        }
    }
    return false;
}

//glues the "bag" with "glueBag" in the current tree decomposition "T"
template <typename T_t>
void sep_glue_bag(std::set<unsigned int> &bag, std::set<unsigned int> &glueBag, T_t &T){
    if(boost::num_vertices(T) == 0)
        boost::add_vertex(T);
    
    typename boost::graph_traits<T_t>::vertex_iterator vertexIt, vertexEnd;    
    for(boost::tie(vertexIt, vertexEnd) = boost::vertices(T); vertexIt != vertexEnd; vertexIt++){
        if(T[*vertexIt].bag == glueBag){
            typename boost::graph_traits<T_t>::vertex_descriptor t_dec_node = boost::add_vertex(T);
            T[t_dec_node].bag = bag;
            boost::add_edge(t_dec_node, *vertexIt, T);
            return;
        }
    }
}

//the main procedure of the seperator algorithm
template <typename G_t, typename T_t>
bool sep_decomp(G_t &G, T_t &T, typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &W, typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &parent, typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> &vertices, std::vector<bool> &disabled, unsigned int k){    
    //tw > k - one could replace this with a better lower bound (see TD_lower_bounds)
    if(boost::num_edges(G) > k*boost::num_vertices(G))
        return false;
    
    //passing if V(G) is a subset of W 
    if(std::includes(W.begin(), W.end(), vertices.begin(), vertices.end()))
        return true;

    std::set<unsigned int> bag1, bag2;
    treedec::descriptor_bag_to_id_bag(G, bag2, parent);

    //trivial decomposition
    if(vertices.size() < 4*k + 2){
        treedec::descriptor_bag_to_id_bag(G, bag1, vertices);
        sep_glue_bag(bag1, bag2, T);
        return true;
    }

    //choose a superset of W of size 3k + 1
    superset(W, vertices, 3*k + 1);
    
    typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> S;
       
    //if a nearly balanced seperator S of W' exists, proceed with the graphs induced by the resulting
    //components and the seperator recursively, add a bag containing the union of W and S to the decomposition,
    //connected with the bag, created in the "parent-call" of the procedure
    if(nearly_balanced_seperator(G, W, S, disabled, k)){
        std::vector<typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> > C;

        for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = S.begin(); sIt != S.end(); sIt++)
            disabled[G[*sIt].id] = true;

        get_components_provided_map(G, C, disabled);

        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> union_W_S;
        std::set_union(W.begin(), W.end(), S.begin(), S.end(), std::inserter(union_W_S, union_W_S.begin()));
        
        //create a bag (W' v S) and connect it with the bag containing parent
        treedec::descriptor_bag_to_id_bag(G, bag1, union_W_S);
        sep_glue_bag(bag1, bag2, T); 
        
        for(unsigned int i = 0; i < C.size(); i++){
            typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> union_C_i_S, is_C_i_W, newW;
            std::set_union(C[i].begin(), C[i].end(), S.begin(), S.end(), std::inserter(union_C_i_S, union_C_i_S.begin()));
    
            std::set_intersection(C[i].begin(), C[i].end(), W.begin(), W.end(), std::inserter(is_C_i_W, is_C_i_W.begin()));
            std::set_union(is_C_i_W.begin(), is_C_i_W.end(), S.begin(), S.end(), std::inserter(newW, newW.begin()));

            std::set<unsigned int> union_C_i_S_;
            treedec::descriptor_bag_to_id_bag(G, union_C_i_S_, union_C_i_S);

            std::vector<bool> disabled_(boost::num_vertices(G), true);
            for(typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor>::iterator sIt = union_C_i_S.begin(); sIt != union_C_i_S.end(); sIt++)
                disabled_[G[*sIt].id] = false;
            
            //reject if no seperator can be found in one of the ongoing recursive calls
            if(!sep_decomp(G, T, newW, union_W_S, union_C_i_S, disabled_, k))
                return false;
        }
        return true;
    }
    return false;
}

//starts the seperator algorithm, and tries k = 0,1,2,.. 
//until the whole graph could be decomposed
template <typename G_t, typename T_t>
void seperator_algorithm(G_t &G, T_t &T){
    unsigned int k = 0;
    bool finished = false;

    typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> vertices;
    typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
    for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++)
        vertices.insert(*vIt);
    
    while(!finished){
        std::vector<bool> disabled(boost::num_vertices(G), false);
        typename std::set<typename boost::graph_traits<G_t>::vertex_descriptor> emptySet, parent;
        finished = sep_decomp(G, T, emptySet, parent, vertices, disabled, k);
        k++;

        if(!finished)
            T.clear();
    }
}

}
