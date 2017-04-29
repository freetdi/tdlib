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
#include <boost/property_map/property_map.hpp>

#include "generic_base.hpp"
#include "generic_elimination_search_overlay.hpp"
#include "graph.hpp"

//#include "generic_elimination_search_configs.hpp"
#include "marker_util.hpp"

#include <iostream>


namespace treedec {

namespace gen_search {

template <typename G_t, class CFG_t, template<class G, class ...> class CFGT_t>
generic_elimination_search_base<G_t, CFG_t, CFGT_t>::
    generic_elimination_search_base(G_t const &g,
                                    unsigned g_lb, unsigned g_ub,
                                    unsigned depth, unsigned nodes_generated,
                                    unsigned orderings_generated)
      : algo1(CFG_t::name()),
        _active(*new std::vector<BOOL>(boost::num_vertices(g), true)),
        _g(*new internal_graph_type(g,
                    boost::make_iterator_property_map(&_active[0],
                               boost::typed_identity_property_map<vertex_descriptor>() ))),
        _best_ordering    (*new std::vector<vd>  (boost::num_vertices(g))),
        _current_ordering (*new std::vector<vd>  (boost::num_vertices(g))),
        _global_lb(g_lb),
        _global_ub(g_ub),
        _depth(depth),
        _nodes_generated(nodes_generated),
        _orderings_generated(orderings_generated),
        _need_cleanup(3)
{
#ifdef DEBUG
    auto p=boost::edges(g);
    for(; p.first!=p.second; ++p.first){
        /// test self loops? no. only G.
        // assert(source!=target)...
    }
#endif
}

// strange: still takes overlay as arg..
template <typename G_t, class CFG_t, template<class G, class ...> class CFGT_t>
generic_elimination_search_base<G_t, CFG_t, CFGT_t>::
    generic_elimination_search_base(internal_graph_type &Overlay, // BUG: exposes graph type
                                    unsigned g_lb, unsigned g_ub,
                                    unsigned depth, unsigned nodes_generated,
                                    unsigned orderings_generated)
      : algo1(CFG_t::name()),
        _g(Overlay),
        _active(*(new std::vector<BOOL>(boost::num_vertices(Overlay), true))),
        _best_ordering    (*(new std::vector<vd>  (boost::num_vertices(Overlay)))),
        _current_ordering (*(new std::vector<vd>  (boost::num_vertices(Overlay)))),
        _global_lb(g_lb),
        _global_ub(g_ub),
        _depth(depth),
        _nodes_generated(nodes_generated),
        _orderings_generated(orderings_generated),
        _need_cleanup(1)
{
#ifdef DEBUG_NOTYET // perhaps not here.
    auto p=boost::edges(_g); // not implemented
    for(; p.first!=p.second; ++p.first){
        /// test self loops? no. only G.
        // assert(source!=target)...
    }
#endif
}

// BUG: optional active etc.
template <typename G_t, class CFG_t, template<class G, class ...> class CFGT_t>
generic_elimination_search_base<G_t, CFG_t, CFGT_t>::
    generic_elimination_search_base(internal_graph_type &Overlay, // BUG: exposes graph type
            std::vector<BOOL>& active, // BUG need normal constructor
                                    std::vector<vd> &best_ordering,
                                    std::vector<vd> &current_ordering,
                                    unsigned g_lb, unsigned g_ub,
                                    unsigned depth, unsigned nodes_generated,
                                    unsigned orderings_generated)
      : algo1(CFG_t::name()),
        _g(Overlay),
        _active(active),
        _best_ordering(best_ordering),
        _current_ordering(current_ordering),
        _global_lb(g_lb),
        _global_ub(g_ub),
        _depth(depth),
        _nodes_generated(nodes_generated),
        _orderings_generated(orderings_generated),
        _need_cleanup(0)
{
#ifdef DEBUG_NOTYET // perhaps not here.
    auto p=boost::edges(_g); // not implemented
    for(; p.first!=p.second; ++p.first){
        /// test self loops? no. only G.
        // assert(source!=target)...
    }
#endif
}

// recursion.
template <typename G_t, class CFG_t, template<class G, class...> class CFGT_t>
generic_elimination_search_base<G_t, CFG_t, CFGT_t>::generic_elimination_search_base(
    generic_elimination_search_base<G_t, CFG_t, CFGT_t>& o)
    : baseclass(o), // good idea?!
      _g(o._g),
      _active(o._active),
      _best_ordering(o._best_ordering),
      _current_ordering(o._current_ordering),
      _global_lb(o._global_lb),
      _global_ub(o._global_ub),
      _depth(o._depth), // ??
      _nodes_generated(o._nodes_generated),
      _orderings_generated(o._orderings_generated),
      _need_cleanup(0)
{}


// BUG: CFG_t is a generic_search config...
// FIXME: do not inherit for IS-IMPLEMENTED-IN-TERMS-OF, see
// http://www.gotw.ca/publications/mill07.htm
template <typename G_t, class OC, template<class G, class ...> class CFGT_t>
class generic_elimination_search_DFS
    : public generic_elimination_search_base<G_t, OC, CFGT_t> {
private:
    // typedef CFG_t UC?
    typedef CFGT_t<G_t> UC;
    typedef generic_elimination_search_base<G_t, OC, CFGT_t> baseclass;
    typedef typename baseclass::vertex_descriptor vd;

public:
    // too many args
    generic_elimination_search_DFS(overlay<G_t, G_t> &Overlay,
            std::vector<BOOL>& active,
            std::vector<vd> &best_ordering, std::vector<vd> &current_ordering)
      : baseclass(Overlay,
              active,
                  best_ordering,
                  current_ordering,
		  0,
                  boost::num_vertices(Overlay),
                  0, 0, 0),
		  local_lb(0),
                  local_ub(0),
                  max_nodes_generated(1),
                  max_orderings_generated(UINT_MAX)
    {}

    generic_elimination_search_DFS(overlay<G_t, G_t> &Overlay)
      : baseclass(Overlay,
		  0,
                  boost::num_vertices(Overlay),
                  0, 0, 0),
		  local_lb(0),
                  local_ub(0),
                  max_nodes_generated(1),
                  max_orderings_generated(UINT_MAX)
    {}

    generic_elimination_search_DFS(G_t const &g)
      : baseclass(g, 0, boost::num_vertices(g),
                  0, 0, 0),
		  local_lb(0),
                  local_ub(0),
                  max_nodes_generated(1),
                  max_orderings_generated(UINT_MAX)
    {}

private: // recursion
    generic_elimination_search_DFS(baseclass& base, unsigned lb, unsigned ub)
        : baseclass(base), local_lb(lb), local_ub(ub),
          max_nodes_generated(UINT_MAX),
          max_orderings_generated(UINT_MAX)
    {
        ++baseclass::_depth;
        ++baseclass::_nodes_generated;
    }

public:
    void do_it();

    void set_max_nodes_generated(unsigned num){ max_nodes_generated = num; }
    void set_max_orderings_generated(unsigned num){ max_orderings_generated = num; }

    unsigned global_lower_bound_bagsize(){ return baseclass::_global_lb; }
    unsigned global_upper_bound_bagsize(){ return baseclass::_global_ub; }

protected:
    unsigned local_lb;  //lb for the current graph
    unsigned local_ub;  //ub for the current graph

    unsigned max_nodes_generated;
    unsigned max_orderings_generated;
};


template <typename G_t, class CFG_t, template<class G, class...> class CFGT_t>
void generic_elimination_search_DFS<G_t, CFG_t, CFGT_t>::do_it()
{
    if(baseclass::_nodes_generated % 1000 == 0){
        std::cout << "#: " << baseclass::_nodes_generated << std::endl;
    }

    baseclass::timer_on();

    if(baseclass::_depth == 0){
        unsigned tmp_global_lb = CFG_t::initial_lb_algo(baseclass::_g.underlying());
        if(tmp_global_lb > baseclass::_global_lb){
            baseclass::_global_lb = tmp_global_lb;
        }

        std::cout << "initial lb: " << baseclass::_global_lb << std::endl;

        unsigned tmp_global_ub = CFG_t::initial_ub_algo(baseclass::_g.underlying(), baseclass::_best_ordering);

        // TODO: proper range container...
        if(tmp_global_ub < baseclass::_global_ub){
            baseclass::_global_ub = tmp_global_ub;
        }
        std::cout << "initial ub: " << baseclass::_global_ub << std::endl;

        // baseclass::bagsize_range().size()==1...?
        if(baseclass::_global_lb == baseclass::_global_ub){
            std::cout << "finished: initial lower bound == initial upper bound" << std::endl;
            ++baseclass::_orderings_generated; //not necessary..
            //returns now
        }
    }

    if(baseclass::_depth == boost::num_vertices(baseclass::_g.underlying())){
        if(local_ub < baseclass::_global_ub){ //this should be always true?!
            std::cout << "found better ordering of width " << local_ub << std::endl;
            baseclass::_global_ub = local_ub;
            baseclass::_best_ordering = baseclass::_current_ordering;
            std::cout << "updated global_ub to " << baseclass::_global_ub << std::endl;

            ++baseclass::_orderings_generated; //ifdef stats?!

            std::vector<vd> tmp_ordering(boost::num_vertices(baseclass::_g.underlying()));
            unsigned ref_width = CFG_t::refiner(baseclass::_g.underlying(), baseclass::_best_ordering, tmp_ordering);

            if(ref_width < baseclass::_global_ub){
                std::cout << "found better ordering after refinement " << ref_width << std::endl;
                baseclass::_global_ub = ref_width;
                baseclass::_best_ordering = tmp_ordering;
                std::cout << "updated global_ub to " << baseclass::_global_ub << std::endl;
           }

        }
        else{
            unreachable(); //should be the case?
        }
    }
    else{
/*
        local_lb = CFG_t::lb_algo(baseclass::_g.underlying());
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
            vd elim_vertex = CFG_t::next(baseclass::_g.underlying(), baseclass::active(), idx);
            if(elim_vertex == CFG_t::INVALID_VERTEX()){
                break;
            }

            if(baseclass::_nodes_generated > max_nodes_generated || baseclass::_orderings_generated > max_orderings_generated){
                break;
            }

            //eliminate elim_vertex:
            //  -do not delete outedges
            //  -store added edges in a container, such that the changes can be undone this way
            //  -if the underlying adj-list graph backend adds edges at the end of the vertex container,
            //   we just have to remember how many edges per node we added and can that pop_back these edges


//            baseclass::eliminate(elim_vertex);
//            unsigned step_width = baseclass::degree(elim_vertex)+1;
            unsigned step_width = baseclass::_g.degree(elim_vertex)+1;

            if(step_width < baseclass::_global_ub){
                baseclass::eliminate(elim_vertex);
                unsigned next_local_ub = (step_width > local_ub)? step_width : local_ub; //local ub is the current width of the ordering

                baseclass::_current_ordering[baseclass::_depth] = elim_vertex;

                //can be seen as a recursion
                generic_elimination_search_DFS nextStep(*this,
                                                        local_lb,
                                                        next_local_ub);

                nextStep.max_nodes_generated = max_nodes_generated;
                nextStep.max_orderings_generated = max_orderings_generated;

                nextStep.do_it();

                baseclass::_nodes_generated = nextStep._nodes_generated; //ifdef stats?!
                baseclass::_orderings_generated = nextStep._orderings_generated;

 //               baseclass::undo_eliminate();

                if(nextStep._global_ub < baseclass::_global_ub){
                    baseclass::_global_ub = nextStep._global_ub; //may have improved

                    //this branch has already width global_ub, so we cant improve here (or we found the exact solution)
                    if(local_ub >= baseclass::_global_ub || baseclass::_global_lb == baseclass::_global_ub){
                        baseclass::undo_eliminate();
                        break; //returns now
                    }
                }
            }
            else{
                //std::cout << "prune branch, since the current branch has higher width than the current best solution" << std::endl;
            }
        }
    }

    baseclass::timer_off();
}


#if 0
template <typename G_t, class CFG_t, template<class G, class...> class CFGT_t>
typename generic_elimination_search_base<G_t, CFG_t, CFGT_t>::adj_range
inline generic_elimination_search_base<G_t, CFG_t, CFGT_t>::adjacent_vertices(
        typename generic_elimination_search_base<G_t, CFG_t, CFGT_t>::vertex_descriptor v) const
{
        active_filter p(active());
        std::pair<overlay_adjacency_iterator, overlay_adjacency_iterator>
            q=boost::adjacent_vertices(v, _g);
        overlay_adjacency_iterator f=q.first;
        overlay_adjacency_iterator s=q.second;
        typedef boost::filter_iterator<active_filter, overlay_adjacency_iterator> FilterIter;
        FilterIter b(p, f, s);
        FilterIter e(p, s, s);

        return std::make_pair(b, e);
}
#endif


} //namespace gen_search

} //namespace treedec


#endif //TD_GENERIC_ELIM_SEARCH_H

// vim:ts=8:sw=4:et
