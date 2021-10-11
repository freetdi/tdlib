// Lukas Larisch, 2014 - 2016
// Felix Salfelder, 2016
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
// Foundation, 51 Franklin Street - Suite 500, Boston, MA 02110-1335, USA.
//
//
//
// greedy heuristics

#ifndef TREEDEC_FILL_IN_OBSOLETE_HPP
#define TREEDEC_FILL_IN_OBSOLETE_HPP

#ifndef TREEDEC_ELIMINATION_ORDERINGS_HPP
#error "not intended to be used like that."
#endif

#include "../algo.hpp"
#include "greedy_base.hpp"
#include "obsolete_greedy_base.hpp"

namespace treedec{
namespace obsolete {
// the fillIn heuristic.
template <typename G_t, template<class G, class...> class CFGT_t=algo::default_config>
class fillIn : public treedec::impl::greedy_heuristic_base<G_t, CFGT_t>{
public: //types
    typedef treedec::impl::greedy_heuristic_base<G_t, CFGT_t> baseclass;
    typedef typename treedec::obsolete::FILL<G_t> fill_type;

    struct fill_update_cb : public graph_callback<G_t>{
        typedef typename baseclass::vertex_descriptor vertex_descriptor;

        fill_update_cb(fill_type* d, G_t const& g) :
            _fill(d), G(g){}

        void operator()(vertex_descriptor v){ untested();
            _fill->q_eval(v);
        }
        void operator()(vertex_descriptor s, vertex_descriptor t) { untested();
            assert(s < t); // likely not. is this necessary below?
            // e has just been inserted.
            BOOST_AUTO(cni, common_out_edges(s, t, G));
            BOOST_AUTO(i, cni.first);
            BOOST_AUTO(e, cni.second);
            for(; i!=e; ++i){ untested();
                assert(*i != s);
                assert(*i != t);
    //            no. maybe theres only half an edge.
    //            assert(boost::edge(boost::source(edg, G), *i, G).second);
    //            assert(boost::edge(boost::target(edg, G), *i, G).second);

                // BUG: *i might be within 1-neighborhood.
                _fill->q_decrement(*i);
            }
        }
    private:
        fill_type* _fill;
        G_t const& G;
    }; // update_cb

public: // construct
    fillIn(G_t &g, unsigned ub=UINT_MAX, bool ignore_isolated_vertices=false)
        : baseclass(g, ub, ignore_isolated_vertices),
         _fill(baseclass::_g), _cb(fill_update_cb(&_fill, baseclass::_g))
    { untested();
    }

    fillIn(G_t &G, bool ignore_isolated_vertices, unsigned ub=-1u)
        : baseclass(G, ub, ignore_isolated_vertices),
          _fill(baseclass::_g), _cb(fill_update_cb(&_fill, baseclass::_g))
    { untested();
    }

public: // implementation
    void initialize(){ untested();
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(baseclass::_g); vIt != vEnd; ++vIt){
            if(boost::out_degree(*vIt, baseclass::_g) == 0){
                if(!baseclass::_iiv){
                    (*baseclass::_o)[baseclass::_i++] = *vIt;
                }
                else{ untested();
                    --baseclass::_num_vert;
                }
            }
        }
    }

    // obs::fillIn::
    void next(typename baseclass::vertex_descriptor &c){
        _fill.check();
        boost::tie(c, baseclass::_min) = _fill.pick_min(0, -1, true);
        _fill.check();
    }

    // obs::fillIn::
    void eliminate(typename baseclass::vertex_descriptor v){
        _fill.mark_neighbors(v, baseclass::_min);

        baseclass::_current_N->resize(boost::out_degree(v, baseclass::_g));

        make_clique_and_detach(v, baseclass::_g, *baseclass::_current_N, &_cb);

        _fill.unmark_neighbours(*baseclass::_current_N);
    }

    void postprocessing(){ untested();
        for(; baseclass::_i < baseclass::_num_vert; ++baseclass::_i){ untested();
            auto v = _fill.pick_min(0, 0, true).first;
            (*baseclass::_o)[baseclass::_i] = v;
        }
    }

private:
    fill_type _fill;
    fill_update_cb _cb;

}; // fillIn

} // obsolete
} // treedec

#endif // guard

// vim:ts=8:sw=4:et
