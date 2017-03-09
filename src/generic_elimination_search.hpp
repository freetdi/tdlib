// Lukas Larisch, 2014 - 2017
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

#ifndef TD_GENERIC_ELIM_SEARCH
#define TD_GENERIC_ELIM_SEARCH

#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>

#include "algo.hpp"
#include "generic_elimination_search_configs.hpp"
#include "generic_elimination_search_overlay.hpp"

#include <iostream>

namespace treedec{

namespace gen_search{

template <typename G_t, typename Olay_t, typename CFG_t>
class generic_elimination_search_base : public treedec::algo::draft::algo1{
public:
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;

    //TODO: better use iterators for elim_vertices
    generic_elimination_search_base(overlay<G_t, Olay_t> &Overlay_input,
                                    std::vector<vd> &best_ordering_input, std::vector<vd> &current_ordering_input,
                                    unsigned g_lb, unsigned g_ub, unsigned depth_input, unsigned nodes_generated_input, unsigned orderings_generated_input)
      : algo1(CFG_t::name()), best_ordering(best_ordering_input),
        global_lb(g_lb), global_ub(g_ub), depth(depth_input), nodes_generated(nodes_generated_input), orderings_generated(orderings_generated_input),
        Overlay(Overlay_input), current_ordering(current_ordering_input)
    {}

    virtual void do_it() = 0;

    void elimination_ordering(std::vector<vd> &ordering) { ordering = best_ordering; }

    unsigned global_lower_bound(){ return global_lb; }
    unsigned global_upper_bound(){ return global_ub; }

    unsigned get_nodes_generated(){ return nodes_generated; }
    unsigned get_orderings_generated(){ return orderings_generated; }

protected:
    overlay<G_t, Olay_t> &Overlay;
    std::vector<vd> &best_ordering;
    std::vector<vd> &current_ordering;

    unsigned global_lb; //lb for the original graph
    unsigned global_ub; //ub for the original graph

    unsigned depth;

    unsigned nodes_generated;
    unsigned orderings_generated;
};


template <typename G_t, typename Olay_t, typename CFG_t>
class generic_elimination_search_DFS : public generic_elimination_search_base<G_t, Olay_t, CFG_t>{
    typedef generic_elimination_search_base<G_t, Olay_t, CFG_t> baseclass;

    typedef typename baseclass::vd vd;

public:
    generic_elimination_search_DFS(overlay<G_t, Olay_t> &Overlay_input,
                                   std::vector<vd> &best_ordering_input, std::vector<vd> &current_ordering_input,
                                   unsigned g_lb, unsigned g_ub, unsigned l_lb, unsigned l_ub, unsigned depth_input, unsigned generated_nodes_input, unsigned generated_orderings_input)
      : generic_elimination_search_base<G_t, Olay_t, CFG_t>(Overlay_input, best_ordering_input, current_ordering_input,
                                                            g_lb, g_ub, depth_input, generated_nodes_input, generated_orderings_input),
        local_lb(l_lb), local_ub(l_ub), max_nodes_generated(UINT_MAX), max_orderings_generated(UINT_MAX){}

    void do_it();

    void set_max_nodes_generated(unsigned num){ max_nodes_generated = num; }
    void set_max_orderings_generated(unsigned num){ max_orderings_generated = num; }

    unsigned global_lower_bound_bagsize(){ return baseclass::global_lb; }
    unsigned global_upper_bound_bagsize(){ return baseclass::global_ub; }

protected:
    unsigned local_lb;  //lb for the current graph
    unsigned local_ub;  //ub for the current graph

    unsigned max_nodes_generated;
    unsigned max_orderings_generated;
};


template <typename G_t, typename Olay_t, typename CFG_t>
void generic_elimination_search_DFS<G_t, Olay_t, CFG_t>::do_it()
{
    if(baseclass::nodes_generated % 1000 == 0){
        std::cout << "#: " << baseclass::nodes_generated << std::endl;
    }

    baseclass::timer_on();

    if(baseclass::depth == 0){
        unsigned tmp_global_lb = CFG_t::initial_lb_algo(baseclass::Overlay.underlying());
        baseclass::global_lb = (tmp_global_lb > baseclass::global_lb)? tmp_global_lb : baseclass::global_lb;
        std::cout << "initial lb: " << baseclass::global_lb << std::endl;

        unsigned tmp_global_ub = CFG_t::initial_ub_algo(baseclass::Overlay.underlying(), baseclass::best_ordering);

        baseclass::global_ub = (tmp_global_ub < baseclass::global_ub)? tmp_global_ub : baseclass::global_ub;
        std::cout << "initial ub: " << baseclass::global_ub << std::endl;

        if(baseclass::global_lb == baseclass::global_ub){
            std::cout << "finished: initial lower bound == initial upper bound" << std::endl;
            ++baseclass::orderings_generated; //not necessary..
            //returns now
        }
    }

    if(baseclass::depth == boost::num_vertices(baseclass::Overlay.underlying())){
        if(local_ub < baseclass::global_ub){ //this should be always true?!
            std::cout << "found better ordering of width " << local_ub << std::endl;
            baseclass::global_ub = local_ub;
            baseclass::best_ordering = baseclass::current_ordering;
            std::cout << "updated global_ub to " << baseclass::global_ub << std::endl;

            ++baseclass::orderings_generated; //ifdef stats?!

            std::vector<vd> tmp_ordering(boost::num_vertices(baseclass::Overlay.underlying()));
            unsigned ref_width = CFG_t::refiner(baseclass::Overlay.underlying(), baseclass::best_ordering, tmp_ordering);

           if(ref_width < baseclass::global_ub){
               std::cout << "found better ordering after refinement " << ref_width << std::endl;
                baseclass::global_ub = ref_width;
                baseclass::best_ordering = tmp_ordering;
                std::cout << "updated global_ub to " << baseclass::global_ub << std::endl;
           }

        }
        else{
            unreachable(); //should be the case?
        }
    }
    else{
/*
        local_lb = CFG_t::lb_algo(baseclass::Overlay.underlying());
        if(local_lb > baseclass::global_ub){
            //can be seen as pruning this branch
            //std::cout << "prune branch since local_lb is greater than the current best solution" << std::endl;
            baseclass::timer_off();
            return;
        }
*/
        unsigned idx = 0;

        //search starts here
        while(true){
            vd elim_vertex = CFG_t::next(baseclass::Overlay.underlying(), baseclass::Overlay.active(), idx);
            if(elim_vertex == CFG_t::INVALID_VERTEX()){
                break;
            }


            if(baseclass::nodes_generated > max_nodes_generated || baseclass::orderings_generated > max_orderings_generated){
                break;
            }

            //eliminate elim_vertex:
            //  -do not delete outedges
            //  -store added edges in a container, such that the changes can be undone this way
            //  -if the underlying adj-list graph backend adds edges at the end of the vertex container,
            //   we just have to remember how many edges per node we added and can that pop_back these edges


            unsigned step_width = baseclass::Overlay.eliminate(elim_vertex)+1;

            if(step_width < baseclass::global_ub){
                unsigned next_local_ub = (step_width > local_ub)? step_width : local_ub; //local ub is the current width of the ordering

                baseclass::current_ordering[baseclass::depth] = elim_vertex;

                //can be seen as a recursion
                generic_elimination_search_DFS nextStep(baseclass::Overlay,
                                                        baseclass::best_ordering,
                                                        baseclass::current_ordering,
                                                        baseclass::global_lb,
                                                        baseclass::global_ub,
                                                        local_lb,
                                                        next_local_ub,
                                                        baseclass::depth+1,
                                                        baseclass::nodes_generated+1,
                                                        baseclass::orderings_generated
                );

                nextStep.max_nodes_generated = max_nodes_generated;
                nextStep.max_orderings_generated = max_orderings_generated;

                nextStep.do_it();

                baseclass::nodes_generated = nextStep.nodes_generated; //ifdef stats?!
                baseclass::orderings_generated = nextStep.orderings_generated;

                if(nextStep.global_ub < baseclass::global_ub){
                    baseclass::global_ub = nextStep.global_ub; //may have improved

                    //this branch has already width global_ub, so we cant improve here (or we found the exact solution)
                    if(local_ub >= baseclass::global_ub || baseclass::global_lb == baseclass::global_ub){
                        baseclass::Overlay.undo_eliminate(elim_vertex);
                        break; //returns now
                    }
                }
            }
            else{
                //std::cout << "prune branch, since the current branch has higher width than the current best solution" << std::endl;
            }

            baseclass::Overlay.undo_eliminate(elim_vertex);

        }
    }

    baseclass::timer_off();
}

} //namespace gen_search

} //namespace treedec

#endif //TD_GENERIC_ELIM_SEARCH
