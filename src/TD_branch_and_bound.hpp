
#ifndef TD_BRANCH_AND_BOUND
#define TD_BRANCH_AND_BOUND

#include <stack>
#include "TD_lower_bounds.hpp"
#include "TD_upper_bounds.hpp"
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


//TODO: pruning... (and reduce data of pruned branches).

unsigned int elim_orderings_count = 0;
unsigned int lb_computed_count = 0;

template <typename G_t, typename BB_t>
unsigned int branching_operator(typename boost::graph_traits<BB_t>::vertex_descriptor v,
          G_t &G, BB_t &BBT, std::stack<typename boost::graph_traits<BB_t>::vertex_descriptor> &leafs, int ub)
{
    if(!BBT[v].active){
        return UINT_MAX;
    }
    else if(BBT[v].elim_ordering.size() == boost::num_vertices(G)){
        return BBT[v].eo_width;
    }
    else{
        //Create a new child of v in BBT for each not eliminated vertex according to BBT[v].elim_ordering.
        unsigned int min = UINT_MAX;
        std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ext_elim_ordering = BBT[v].elim_ordering;
        ext_elim_ordering.push_back(0);
        unsigned int idx = ext_elim_ordering.size()-1;

        for(unsigned int i = 0; i < BBT[v].used_so_far.size(); i++){ elim_orderings_count++;
            if(!BBT[v].used_so_far[i]){
                ext_elim_ordering[idx] = i;

                G_t G_;
                boost::copy_graph(BBT[v].graph, G_);
                unsigned int w = noboost::eliminate_vertex(i, G_); //TODO: add an edge between v and w if |N(v) ^ N(w)| > ub after this command.
                w = (w < BBT[v].eo_width) ? BBT[v].eo_width : w;

                if(w > ub){
                   continue;
                }
                else{ lb_computed_count++;
                    unsigned int lb = (unsigned int)treedec::lb::deltaC_least_c(G_); //deltaC_least_c(H, ub) would fasten this.

                    //Since this partial elimiation ordering causes higher width as the currently best
                    //known upper bound, discard it.
                    if(lb > ub){
                        continue;
                    }

/*
                    std::cout << "elim " << i << " (id=" << G[i].id << "): ";
                    std::cout << "w=" <<  w << " (ub= " << ub << ")" << std::endl;
*/

                    typename boost::graph_traits<BB_t>::vertex_descriptor child = boost::add_vertex(BBT);
                    boost::add_edge(v, child, BBT);
                    BBT[child].graph = G_;
                    BBT[child].elim_ordering = ext_elim_ordering;
                    BBT[child].used_so_far = BBT[v].used_so_far;
                    BBT[child].used_so_far[i] = true;
                    BBT[child].eo_width = w;
                    BBT[child].active = true;

                    leafs.push(child);

                    if(w < min){
                        min = w;
                    }

/*
                    std::cout << "child: " << "";
                    for(unsigned int i = 0; i < ext_elim_ordering.size(); i++){
                        std::cout << ext_elim_ordering[i] << " ";
                    } std::cout << std::endl;
*/
                }
            }
        }
        BBT[v].active = false;

        return min;
    }
}

template <typename G_t>
int branch_and_bound(G_t &G){
    if(boost::num_vertices(G) == 0){
        return -1;
    }
    else if(boost::num_edges(G) == 0){
        return 0;
    }

    elim_orderings_count = 0;
    lb_computed_count = 0;

    //Branch and bound tree.
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, bbt_bag<G_t> > BB_t;
    BB_t BBT;

    //Store the currently active leafs o the branch and bound tree in 'leafs'.
    std::stack<typename boost::graph_traits<BB_t>::vertex_descriptor> leafs;

    //Create the root of the branch and bound tree.
    typename boost::graph_traits<BB_t>::vertex_descriptor root = boost::add_vertex(BBT);
    BBT[root].graph = G; //TODO: add an edge between v and w if |N(v) ^ N(w)| > ub
    BBT[root].used_so_far.assign(boost::num_vertices(G), false);
    BBT[root].eo_width = 0;
    BBT[root].active = true;

    leafs.push(root);

    /* TODO: use best heuristics */

    //Use the lower bound provided by the deltaC_least_c-heuristic.
    unsigned int lb = (unsigned int)treedec::lb::deltaC_least_c(G);

    //Use the upper bound provided by the minFill-heuristic.
    unsigned int ub = treedec::ub::minFill(G);

    if(lb == ub){
        

    //Postorder traversal.
    while(!leafs.empty()){
        //The minimum caused width of a current leaf node in BBT.
        unsigned int minmax = UINT_MAX;

        typename boost::graph_traits<BB_t>::vertex_descriptor cur = leafs.top();

        leafs.pop();
        unsigned int minmax_cur = branching_operator(cur, G, BBT, leafs, ub);
        if(minmax_cur < minmax){
            minmax = minmax_cur;
        }

        //If a partial elimination ordering is a complete elimination ordering
        //and causes lower width than the currently best elimination ordering, update ub.
        if(BBT[cur].elim_ordering.size() == boost::num_vertices(G)){
            ub = (ub < minmax)? ub : minmax;

            //The treewidth has been computed.
            if(ub == lb){
                break;
            }
        }
    }

    std::cout << "computed elim_orderings: " << elim_orderings_count << std::endl;
    std::cout << "computed lbs: " << lb_computed_count << std::endl;

    std::cout << "tw: " << ub << std::endl;
    return ub;
}

} //namespace treedec

#endif //TD_BRANCH_AND_BOUND
