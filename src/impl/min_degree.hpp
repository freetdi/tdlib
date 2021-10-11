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
// mindegree heuristic

#ifndef TREEDEC_MIN_DEGREE_HPP
#define TREEDEC_MIN_DEGREE_HPP

#ifndef TREEDEC_ELIMINATION_ORDERINGS_HPP
#error "not intended to be used like that."
#endif

#include "../algo.hpp"
#include "greedy_base.hpp"
#include "obsolete_greedy_base.hpp"

namespace treedec{

namespace impl{

template <typename G_t, template<class G, class...> class CFG=algo::default_config>
class minDegree : public greedy_heuristic_base<G_t, CFG>{
public:
    typedef greedy_heuristic_base<G_t, CFG> baseclass;

    typedef typename deg_chooser<G_t>::type degs_type;

    minDegree(G_t &g,
                    unsigned ub=UINT_MAX, bool ignore_isolated_vertices=false)
        : baseclass(g, ub, ignore_isolated_vertices),
         _degs(baseclass::_g)
    { untested();
    }

    minDegree(G_t &G, bool ignore_isolated_vertices)
        : baseclass(G, -1u, ignore_isolated_vertices),
          _degs(baseclass::_g)
    { untested();
    }

#if 0 // base
    void get_elimination_ordering(){ untested();
        // incomplete()
    }
#endif

    void initialize(){ untested();
        auto zerodegbag1=MOVE(_degs.detach_bag(0));
        BOOST_AUTO(it, zerodegbag1.begin());

        if(!baseclass::_iiv){
            for(; it!=zerodegbag1.end(); ++it){ untested();
                (*baseclass::_o)[baseclass::_i++] = *it;
            }
        }else{ untested();
            baseclass::_num_vert -= zerodegbag1.size();
        }

        baseclass::_min = 1;
    }

    void next(typename baseclass::vertex_descriptor &c){
        if(baseclass::_min>1){
            --baseclass::_min;
        }

        boost::tie(c, baseclass::_min) = _degs.pick_min(baseclass::_min, baseclass::_num_vert);
    }

    // md::
    void eliminate(typename baseclass::vertex_descriptor v){
        typename baseclass::adjacency_iterator I, E;
        for(boost::tie(I, E) = boost::adjacent_vertices(v, baseclass::_g); I!=E; ++I){
            assert(*I!=v); // no self loops...
            typename baseclass::vertex_descriptor w=*I;
            _degs.unlink(w);
        }

        baseclass::_current_N->resize(boost::out_degree(v, baseclass::_g));

        make_clique_and_detach(v, baseclass::_g, *baseclass::_current_N);

        redegree(NULL, baseclass::_g, *baseclass::_current_N, _degs);
        _degs.unlink(v, baseclass::_min);
        _degs.flush();
    }

    void postprocessing(){ untested();
        auto zerodegbag=MOVE(_degs.detach_bag(0));
        BOOST_AUTO(it, zerodegbag.begin());

        for(; it!=zerodegbag.end(); ++it){ untested();
            (*baseclass::_o)[baseclass::_i++] = *it;
        }
    }

private:
    degs_type _degs;

}; // minDegree

} // impl
} // treedec

#endif // guard

// vim:ts=8:sw=4:et
