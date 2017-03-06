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

#include <iostream>

namespace treedec{

namespace gen_search{

template <typename G_t, typename CFG_t>
class generic_elimination_search_base : public treedec::algo::draft::algo1{
public:
    //TODO: better use iterators for elim_vertices
    generic_elimination_search_base(G_t &G_input,
                                    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &best_ordering_input,
                                    std::vector<bool> &active_input,
                                    unsigned g_lb, unsigned g_ub, unsigned depth_input, unsigned nodes_generated_input, unsigned orderings_generated_input)
      : algo1(CFG_t::name()), G(G_input), best_ordering(best_ordering_input), active(active_input),
        global_lb(g_lb), global_ub(g_ub), depth(depth_input), nodes_generated(nodes_generated_input), orderings_generated(orderings_generated_input)
    {}

    virtual void do_it() = 0;
    virtual void elimination_ordering(std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &ordering) = 0;

    unsigned global_lower_bound(){ return global_lb; }
    unsigned global_upper_bound(){ return global_ub; }

protected:
    G_t &G;
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &best_ordering;
    std::vector<bool> active; //if next() works in O(1), this is not needed, or?! TODO: develop a better interface for next()

    unsigned global_lb; //lb for the original graph
    unsigned global_ub; //ub for the original graph

    unsigned depth;

public:
    unsigned get_nodes_generated(){ return nodes_generated; }
    unsigned get_orderings_generated(){ return orderings_generated; }

protected:
    unsigned nodes_generated;
    unsigned orderings_generated;
};

template <typename G_t, typename CFG_t>
class generic_elimination_search_DFS : public generic_elimination_search_base<G_t, CFG_t>{
    typedef generic_elimination_search_base<G_t, CFG_t> baseclass;

    typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;

public:
    generic_elimination_search_DFS(G_t &G,
                                   std::vector<vd> &best_ordering_input,
                                   std::vector<bool> &active_input,
                                   unsigned g_lb, unsigned g_ub, unsigned l_lb, unsigned l_ub, unsigned depth_input, unsigned generated_nodes_input, unsigned generated_orderings_input)
      : generic_elimination_search_base<G_t, CFG_t>(G, best_ordering_input, active_input, g_lb, g_ub, depth_input, generated_nodes_input, generated_orderings_input),
        local_lb(l_lb), local_ub(l_ub), max_nodes_generated(UINT_MAX), max_orderings_generated(UINT_MAX){}

    void do_it();

    void elimination_ordering(std::vector<vd> &ordering);

    void set_max_nodes_generated(unsigned num);
    void set_max_orderings_generated(unsigned num);

    unsigned global_lower_bound_bagsize();
    unsigned global_upper_bound_bagsize();

private:
    unsigned local_lb;  //lb for the current graph
    unsigned local_ub;  //ub for the current graph

    unsigned max_nodes_generated;
    unsigned max_orderings_generated;
};




template <typename G_t>
unsigned eliminate(
  typename boost::graph_traits<G_t>::vertex_descriptor elim_vertex,
  G_t &G,
  std::vector<bool> &active,
  std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &changes_container)
{
    active[elim_vertex]= false;

    unsigned actual_degree = 0;

    typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;
    for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(elim_vertex, G); nIt1 != nEnd; ++nIt1){
        if(!active[*nIt1]){
            continue;
        }

        ++actual_degree;

        nIt2 = nIt1;
        ++nIt2;
        for(; nIt2 != nEnd; ++nIt2){
            if(!active[*nIt2]){
                continue;
            }
            if(!boost::edge(*nIt1, *nIt2, G).second){
                boost::add_edge(*nIt1, *nIt2, G);
                boost::add_edge(*nIt2, *nIt1, G);
                changes_container.push_back(*nIt1);
                changes_container.push_back(*nIt2);
            }
        }
    }
    return actual_degree;
}

template <typename G_t>
void undo_eliminate(
  typename boost::graph_traits<G_t>::vertex_descriptor elim_vertex,
  G_t &G,
  std::vector<bool> &active,
  std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &changes_container)
{
    active[elim_vertex]= true;
    unsigned c = changes_container.size() >> 1;
    for(unsigned i = 0; i < c; ++i){
        typename boost::graph_traits<G_t>::vertex_descriptor v1 = changes_container.back();
        changes_container.pop_back();
        typename boost::graph_traits<G_t>::vertex_descriptor v2 = changes_container.back();
        changes_container.pop_back();

        boost::remove_edge(v1, v2, G);
        boost::remove_edge(v2, v1, G);
    }
}


template <typename G_t, typename CFG_t>
void generic_elimination_search_DFS<G_t, CFG_t>::do_it()
{
    baseclass::timer_on();

    //global_ub may have changes in the meantime, so we may cut of the search for this branch here
    if(local_ub >= baseclass::global_ub || baseclass::global_lb == baseclass::global_ub){
        baseclass::timer_off();
        return;
    }

    //std::cout << "depth: " << baseclass::depth << std::endl;
    //std::cout << "local ub: " << local_ub << std::endl;
    //std::cout << "global ub: " << baseclass::global_ub << std::endl;

    if(baseclass::depth == 0){
        unsigned tmp_global_lb = CFG_t::initial_lb_algo(baseclass::G);
        baseclass::global_lb = (tmp_global_lb > baseclass::global_lb)? tmp_global_lb : baseclass::global_lb;
        std::cout << "initial lb: " << baseclass::global_lb << std::endl;

        unsigned tmp_global_ub = CFG_t::initial_ub_algo(baseclass::G, baseclass::best_ordering);

        baseclass::global_ub = (tmp_global_ub < baseclass::global_ub)? tmp_global_ub : baseclass::global_ub;
        std::cout << "initial ub: " << baseclass::global_ub << std::endl;

        if(baseclass::global_lb == baseclass::global_ub){
            std::cout << "finished: initial lower bound == initial upper bound" << std::endl;
            ++baseclass::orderings_generated; //not necessary..

            baseclass::timer_off();
            return;
        }
    }

    if(baseclass::depth == boost::num_vertices(baseclass::G)){
        if(local_ub < baseclass::global_ub){
            //TODO: baseclass::best_ordering = local_ordering; TODO:propagate ordering
            std::cout << "found better ordering of width " << local_ub << std::endl;
            baseclass::global_ub = local_ub;
            std::cout << "updated global_ub to " << local_ub << std::endl;

            ++baseclass::orderings_generated; //ifdef stats?!
        }
        else{
            unreachable(); //should be the case?
        }
    }
    else{
        local_lb = CFG_t::lb_algo(baseclass::G);
        if(local_lb > baseclass::global_ub){
            //can be seen as pruning this branch
            //std::cout << "prune branch since local_lb is greater than the current best solution" << std::endl;
            baseclass::timer_off();
            return;
        }

        unsigned idx = 0;

        //search starts here
        while(true){
            typename boost::graph_traits<G_t>::vertex_descriptor elim_vertex = CFG_t::next(baseclass::G, baseclass::active, idx);
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


            std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> changes_container;
            unsigned step_width = eliminate(elim_vertex, baseclass::G, baseclass::active, changes_container)+1;

            //std::cout << "depth: " << baseclass::depth << std::endl;
            //std::cout << "elim_vertex: " << elim_vertex << std::endl;

            if(step_width <= baseclass::global_ub){
                unsigned next_local_ub = (step_width > local_ub)? step_width : local_ub; //local ub is the current width of the ordering

                //can be seen as a recursion
                generic_elimination_search_DFS nextStep(baseclass::G,
                                                        baseclass::best_ordering,
                                                        baseclass::active,
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
                    baseclass::best_ordering[baseclass::depth] = elim_vertex;

/* not necessary if best_ordering is a &
                    for(unsigned i = baseclass::depth+1; i < boost::num_vertices(baseclass::G); ++i){
                        baseclass::best_ordering[i] = nextStep.best_ordering[i];
                    }
*/
                }
            }
            else{
                //std::cout << "prune branch, since the current branch has higher width than the current best solution" << std::endl;
            }
            undo_eliminate(elim_vertex, baseclass::G, baseclass::active, changes_container);
        }
    }

    baseclass::timer_off();
}


template <typename G_t, typename CFG_t>
void generic_elimination_search_DFS<G_t, CFG_t>::elimination_ordering(std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> &elimination_ordering)
{
    elimination_ordering = baseclass::best_ordering;
}


template <typename G_t, typename CFG_t>
void generic_elimination_search_DFS<G_t, CFG_t>::set_max_nodes_generated(unsigned num)
{
    max_nodes_generated = num;
}

template <typename G_t, typename CFG_t>
void generic_elimination_search_DFS<G_t, CFG_t>::set_max_orderings_generated(unsigned num)
{
    max_orderings_generated = num;
}


template <typename G_t, typename CFG_t>
unsigned generic_elimination_search_DFS<G_t, CFG_t>::global_lower_bound_bagsize()
{
    return baseclass::global_lb;
}

template <typename G_t, typename CFG_t>
unsigned generic_elimination_search_DFS<G_t, CFG_t>::global_upper_bound_bagsize()
{
    return baseclass::global_ub;
}


template <typename G_t>
void generic_elimination_search_test1(G_t &G){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering(boost::num_vertices(G), 9);
    std::vector<bool> active(boost::num_vertices(G), true);

    generic_elimination_search_DFS<G_t, configs::CFG_DFS_1<G_t> > //TODO: constructor...
       generic_elim_DFS_test
              (G,
               ordering,
               active,
               0,                         //global_lb
               boost::num_vertices(G),    //global_ub
               0,
               0,
               0,
               1,			//nodes generated
               0			//orderings generated
       );

    generic_elim_DFS_test.do_it();

    std::cout << "lower bound: " << generic_elim_DFS_test.global_lower_bound_bagsize() << std::endl;
    std::cout << "upper bound: " << generic_elim_DFS_test.global_upper_bound_bagsize() << std::endl;
    std::cout << "nodes generated: " << generic_elim_DFS_test.get_nodes_generated() << std::endl;
    std::cout << "orderings generated: " << generic_elim_DFS_test.get_orderings_generated() << std::endl;

    std::cout << "ordering: " << std::endl;
    for(unsigned i = 0; i < ordering.size(); ++i){
        std::cout << ordering[i] << " ";
    } std::cout << std::endl;

    //std::cout << "time: " << generic_elim_DFS_test.get_runtime() << std::endl << std::endl;

    unsigned w = treedec::get_width_of_elimination_ordering(G, ordering);
    std::cout << "width of elimination ordering (check): " << w+1 << std::endl;
    assert(w == generic_elim_DFS_test.global_upper_bound_bagsize());
}

template <typename G_t>
void generic_elimination_search_test2(G_t &G){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering(boost::num_vertices(G));
    std::vector<bool> active(boost::num_vertices(G), true);

    generic_elimination_search_DFS<G_t, configs::CFG_DFS_1<G_t> > //TODO: constructor...
       generic_elim_DFS_test
              (G,
               ordering,
               active,
               0,                         //global_lb
               boost::num_vertices(G),    //global_ub
               0,
               0,
               0,
               1,			//nodes generated
               0			//orderings generated
       );

    generic_elim_DFS_test.set_max_nodes_generated(100000);
    generic_elim_DFS_test.set_max_orderings_generated(100);

    generic_elim_DFS_test.do_it();

    std::cout << "lower bound: " << generic_elim_DFS_test.global_lower_bound_bagsize() << std::endl;
    std::cout << "upper bound: " << generic_elim_DFS_test.global_upper_bound_bagsize() << std::endl;
    std::cout << "nodes generated: " << generic_elim_DFS_test.get_nodes_generated() << std::endl;
    std::cout << "orderings generated: " << generic_elim_DFS_test.get_orderings_generated() << std::endl;

    std::cout << "ordering: " << std::endl;
    for(unsigned i = 0; i < ordering.size(); ++i){
        std::cout << ordering[i] << " ";
    } std::cout << std::endl;

    //std::cout << "time: " << generic_elim_DFS_test.get_runtime() << std::endl << std::endl;

    unsigned w = treedec::get_width_of_elimination_ordering(G, ordering);
    std::cout << "width of elimination ordering (check): " << w+1 << std::endl;
    assert(w == generic_elim_DFS_test.global_upper_bound_bagsize());
}

template <typename G_t>
void generic_elimination_search_test3(G_t &G){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering(boost::num_vertices(G));
    std::vector<bool> active(boost::num_vertices(G), true);

    generic_elimination_search_DFS<G_t, configs::CFG_DFS_2<G_t> > //TODO: constructor...
       generic_elim_DFS_test
              (G,
               ordering,
               active,
               0,                         //global_lb
               boost::num_vertices(G),    //global_ub
               0,
               0,
               0,
               1,			//nodes generated
               0			//orderings generated
       );

    generic_elim_DFS_test.do_it();

    std::cout << "lower bound: " << generic_elim_DFS_test.global_lower_bound_bagsize() << std::endl;
    std::cout << "upper bound: " << generic_elim_DFS_test.global_upper_bound_bagsize() << std::endl;
    std::cout << "nodes generated: " << generic_elim_DFS_test.get_nodes_generated() << std::endl;
    std::cout << "orderings generated: " << generic_elim_DFS_test.get_orderings_generated() << std::endl;

    std::cout << "ordering: " << std::endl;
    for(unsigned i = 0; i < ordering.size(); ++i){
        std::cout << ordering[i] << " ";
    } std::cout << std::endl;

    //std::cout << "time: " << generic_elim_DFS_test.get_runtime() << std::endl << std::endl;

    unsigned w = treedec::get_width_of_elimination_ordering(G, ordering);
    std::cout << "width of elimination ordering (check): " << w+1 << std::endl;
    assert(w == generic_elim_DFS_test.global_upper_bound_bagsize());
}

template <typename G_t>
void generic_elimination_search_test4(G_t &G){
    std::vector<typename boost::graph_traits<G_t>::vertex_descriptor> ordering(boost::num_vertices(G));
    std::vector<bool> active(boost::num_vertices(G), true);

    generic_elimination_search_DFS<G_t, configs::CFG_DFS_2<G_t> > //TODO: constructor...
       generic_elim_DFS_test
              (G,
               ordering,
               active,
               0,                         //global_lb
               boost::num_vertices(G),    //global_ub
               0,
               0,
               0,
               1,			//nodes generated
               0			//orderings generated
       );

    generic_elim_DFS_test.set_max_nodes_generated(10000);
    generic_elim_DFS_test.set_max_orderings_generated(100);

    generic_elim_DFS_test.do_it();

    std::cout << "lower bound: " << generic_elim_DFS_test.global_lower_bound_bagsize() << std::endl;
    std::cout << "upper bound: " << generic_elim_DFS_test.global_upper_bound_bagsize() << std::endl;
    std::cout << "nodes generated: " << generic_elim_DFS_test.get_nodes_generated() << std::endl;
    std::cout << "orderings generated: " << generic_elim_DFS_test.get_orderings_generated() << std::endl;

    std::cout << "ordering: " << std::endl;
    for(unsigned i = 0; i < ordering.size(); ++i){
        std::cout << ordering[i] << " ";
    } std::cout << std::endl;

    //std::cout << "time: " << generic_elim_DFS_test.get_runtime() << std::endl << std::endl;

    unsigned w = treedec::get_width_of_elimination_ordering(G, ordering);
    std::cout << "width of elimination ordering (check): " << w+1 << std::endl;
    assert(w == generic_elim_DFS_test.global_upper_bound_bagsize());
}


} //namespace gen_search

} //namespace treedec

#endif //TD_GENERIC_ELIM_SEARCH
