// Lukas Larisch, 2014 - 2017
// Felix Salfelder 2017
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

#ifndef TD_GENERIC_ELIM_SEARCH_H
#define TD_GENERIC_ELIM_SEARCH_H

#include <vector>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>

#include "generic_base.hpp"
#include "generic_elimination_search_overlay.hpp"
#include "graph.hpp"

#include "generic_elimination_search_configs.hpp"
#include "marker_util.hpp"

#include <iostream>


namespace treedec {

namespace gen_search {

template <typename G_t, template<class G, class ...> class CFGT_t>
generic_elimination_search_base<G_t, CFGT_t>::
    generic_elimination_search_base(internal_graph_type &Overlay_input, // BUG: exposes graph type
            std::vector<BOOL>& active, // BUG need normal constructor
                                    std::vector<vd> &best_ordering_input,
                                    std::vector<vd> &current_ordering_input,
                                    unsigned g_lb, unsigned g_ub,
                                    unsigned depth_input, unsigned nodes_generated_input,
                                    unsigned orderings_generated_input)
      : algo1(CFG_t::name()),
        Overlay(Overlay_input),
        _active(active),
        best_ordering(best_ordering_input),
        current_ordering(current_ordering_input),
        global_lb(g_lb), global_ub(g_ub), depth(depth_input),
        nodes_generated(nodes_generated_input), orderings_generated(orderings_generated_input),
        _marker(boost::num_vertices(Overlay_input))
{
#ifdef DEBUG_NOTYET // perhaps not here.
    auto p=boost::edges(Overlay); // not implemented
    for(; p.first!=p.second; ++p.first){
        /// test self loops? no. only G.
        // assert(source!=target)...
    }
#endif
}

template <typename G_t, template<class G, class...> class CFGT_t>
generic_elimination_search_base<G_t, CFGT_t>::generic_elimination_search_base(
    generic_elimination_search_base<G_t, CFGT_t>& o)
    : baseclass(o),
    Overlay(o.Overlay),
    _active(o._active),
    best_ordering(o.best_ordering),
    current_ordering(o.current_ordering),
    global_lb(o.global_lb),
    global_ub(o.global_ub),
    depth(o.depth), // ??
    nodes_generated(o.nodes_generated),
     orderings_generated(o.orderings_generated),
    _marker(boost::num_vertices(o.Overlay))
{ untested();
}

template <typename G_t, template<class G, class ...> class CFGT_t>
unsigned generic_elimination_search_base<G_t, CFGT_t>::eliminate(
		typename generic_elimination_search_base<G_t, CFGT_t>::vertex_descriptor elim_vertex)
{
	using draft::concat_iterator;
	active()[elim_vertex] = false;

	unsigned actual_degree = 0;

	// auto p=boost::adjacent_vertices(elim_vertex, *this);
	auto p=adjacent_vertices(elim_vertex);
	for(; p.first!=p.second; ++p.first){
		assert(active()[*p.first]);

		++actual_degree;

		_marker.clear();
                // mark_smaller_neighbours(_marker, *p.first, *this); doesntwork
                mark_smaller_neighbours(_marker, *p.first, Overlay, active());

                auto q=adjacent_vertices(elim_vertex);
                for(; q.first!=q.second; ++q.first){
                    assert(active()[*q.first]);
                    if(*q.first>=*p.first){
                        // skip. TODO: more efficient skip
                    }else if(_marker.is_marked(*q.first)){
                        // done.
                    }else{
                        // std::cerr << "ae " << *p.first << *q.first << "\n";
                        assert(!edge(*p.first, *q.first, Overlay).second);
                        Overlay.add_edge(*p.first, *q.first);

                        // treedec graph iface enforces this. (.. should)
                        assert(edge(*p.first, *q.first, Overlay).second);
                        assert(edge(*q.first, *p.first, Overlay).second);
                    }
		}
	}

	Overlay.commit();
	return actual_degree;
}

template <typename G_t, template<class G, class ...> class CFGT_t>
void generic_elimination_search_base<G_t, CFGT_t>::undo_eliminate(
		typename generic_elimination_search_base<G_t, CFGT_t>::vertex_descriptor elim_vertex)
{
    Overlay.reset(1);
    active()[elim_vertex] = true;
}


template <typename G_t, template<class G, class ...> class CFGT_t>
class generic_elimination_search_DFS : public generic_elimination_search_base<G_t, CFGT_t>{
	 typedef CFGT_t<G_t> CFG_t;
    typedef generic_elimination_search_base<G_t, CFGT_t> baseclass;
    typedef typename baseclass::vertex_descriptor vd;

public:
    // BUG, too messy
    generic_elimination_search_DFS(overlay<G_t, G_t> &Overlay_input,
            std::vector<BOOL>& active,
            std::vector<vd> &best_ordering_input, std::vector<vd> &current_ordering_input,
            unsigned g_lb, unsigned g_ub, unsigned l_lb,
            unsigned l_ub, unsigned depth_input,
            unsigned generated_nodes_input, unsigned generated_orderings_input)
      : baseclass(Overlay_input,
              active,
                  best_ordering_input,
                  current_ordering_input,
		  g_lb, g_ub, depth_input,
                  generated_nodes_input,
                  generated_orderings_input),
		  local_lb(l_lb), local_ub(l_ub), max_nodes_generated(UINT_MAX),
                  max_orderings_generated(UINT_MAX)
    { untested();
    }
private: // recursion
    generic_elimination_search_DFS(baseclass& base, unsigned a, unsigned b)
        : baseclass(base), local_lb(a), local_ub(b),
          max_nodes_generated(UINT_MAX),
          max_orderings_generated(UINT_MAX)
    { untested();

        ++baseclass::depth;
        ++baseclass::nodes_generated;
    }
public:

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


template <typename G_t, template<class G, class...> class CFGT_t>
void generic_elimination_search_DFS<G_t, CFGT_t>::do_it()
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
            vd elim_vertex = CFG_t::next(baseclass::Overlay.underlying(), baseclass::active(), idx);
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


            unsigned step_width = baseclass::eliminate(elim_vertex)+1;

            if(step_width < baseclass::global_ub){
                unsigned next_local_ub = (step_width > local_ub)? step_width : local_ub; //local ub is the current width of the ordering

                baseclass::current_ordering[baseclass::depth] = elim_vertex;

                //can be seen as a recursion
                generic_elimination_search_DFS nextStep(*this,
                                                        local_lb,
                                                        next_local_ub);

                nextStep.max_nodes_generated = max_nodes_generated;
                nextStep.max_orderings_generated = max_orderings_generated;

                nextStep.do_it();

                baseclass::nodes_generated = nextStep.nodes_generated; //ifdef stats?!
                baseclass::orderings_generated = nextStep.orderings_generated;

                if(nextStep.global_ub < baseclass::global_ub){
                    baseclass::global_ub = nextStep.global_ub; //may have improved

                    //this branch has already width global_ub, so we cant improve here (or we found the exact solution)
                    if(local_ub >= baseclass::global_ub || baseclass::global_lb == baseclass::global_ub){
                        baseclass::undo_eliminate(elim_vertex);
                        break; //returns now
                    }
                }
            }
            else{
                //std::cout << "prune branch, since the current branch has higher width than the current best solution" << std::endl;
            }

            baseclass::undo_eliminate(elim_vertex);

        }
    }

    baseclass::timer_off();
}


template <typename G_t, template<class G, class...> class CFGT_t>
typename generic_elimination_search_base<G_t, CFGT_t>::adj_range
inline generic_elimination_search_base<G_t, CFGT_t>::adjacent_vertices(
        typename generic_elimination_search_base<G_t, CFGT_t>::vertex_descriptor v) const
{
        active_filter p(active());
        std::pair<overlay_adjacency_iterator, overlay_adjacency_iterator>
            q=boost::adjacent_vertices(v, Overlay);
        overlay_adjacency_iterator f=q.first;
        overlay_adjacency_iterator s=q.second;
        typedef boost::filter_iterator<active_filter, overlay_adjacency_iterator> FilterIter;
        FilterIter b(p, f, s);
        FilterIter e(p, s, s);

        return std::make_pair(b, e);
}


} //namespace gen_search

} //namespace treedec


#endif //TD_GENERIC_ELIM_SEARCH_H
// vim:ts=8:sw=4:et
