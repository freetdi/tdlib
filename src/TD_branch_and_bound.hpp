
#ifndef TD_BRANCH_AND_BOUND
#define TD_BRANCH_AND_BOUND

#include "TD_elimination_orderings.hpp"

namespace treedec{

template <typename G_t>
struct bbt_bag{
    G_t graph;
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> elim_ordering;
    std::vector<bool> used_so_far;
    unsigned int eo_width;
    bool active;
};


template <typename G_t, typename BB_t>
unsigned int branching_operator(typename boost::graph_traits<BB_t>::vertex_descriptor v,
          G_t &G, BB_t &BBT, std::stack<typename boost::graph_traits<BB_t>::vertex_descriptor> &leafs; int ub)
{
    if(!BBT[v].active){
        return 0;
    }
    else{
        //Create a new child of v in BBT for each not eliminated vertex according to BBT[v].elim_ordering.

        unsigned int min = UINT_MAX;
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ext_elim_ordering = BBT[v].elim_ordering;
        ext_elim_ordering.push_back(0);
        unsigned int idx = ext_elim_ordering.size()-1;

        for(unsigned int i = 0; i < BBT[v].used_so_far; i++){
            if(!BBT[v].used_so_far[i]){
                ext_elim_ordering[idx] = i;
                G_t G_;
                boost::copy_graph(BBT[v].graph, G_);
                unsigned int w = noboost::eliminate_vertex(i, G_);
                if(w > ub){
                   continue;
                }
                else{

                    //Todo: Compute a lb for this ordering. If lb > ub, then do not add the child.

                    typename boost::graph_traits<BB_t>::vertex_descriptor child = boost::add_vertex(BBT);
                    boost::add_edge(v, child, BBT);
                    BBT[child].graph = G_;
                    BBT[child].elim_ordering = ext_elim_ordering;
                    BBT[child].used_so_far = BBT[v].used_so_far;
                    BBT[child].used_so_far[i] = true;
                    BBT[child].eo_width = = w;
                    BBT[child].active = true;

                    leafs.push(child);

                    if(w < min){
                        min = w;
                    }
                }
            }
        }
        BBT[v].active = false;

        return w;
    }
}

template <typename G_t>
void branch_and_bound(G_t &G){
    //Branch and bound tree.
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, bbt_bag<G_t> > BB_t;

    //Store the currently active leafs o the branch and bound tree in 'leafs'.
    std::stack<typename boost::graph_traits<BB_t>::vertex_descriptor> leafs;

    //Create the root of the branch and bound tree.
    typename boost::graph_traits<BB_t>::vertex_descriptor v = boost::add_vertex(T);
    BBT[v].graph = G;
    BBT[v].used_so_far.assign(boost::num_vertices(G), false);
    BBT[v].eo_width = 0;
    BBT[v].active = true;

    leaf.push(v);

    //Use the upper bound provided by the minFill-heuristic.
    unsigned int ub = treedec::ub::minFill(G);

    //The minimum caused width of a current leaf node in BBT.
    unsigned int minmin = 0;
}

} //namespace treedec

#endif //TD_BRANCH_AND_BOUND
