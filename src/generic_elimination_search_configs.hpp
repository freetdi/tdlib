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

#ifndef TD_GENERIC_ELIM_SEARCH_CONFIGS
#define TD_GENERIC_ELIM_SEARCH_CONFIGS

#include <vector>
#include <limits.h>

#include "lower_bounds.hpp"
#include "elimination_orderings.hpp"
#include "postprocessing.hpp"
#include "generic_elimination_search.hpp"

#include <iostream>

// "virtual overloads" for algos derived from gen_search_base.

namespace treedec{

namespace gen_search{

namespace configs{
using treedec::gen_search::generic_elimination_search_DFS;
template <typename G_t, template<class G, class ...> class cfg>
struct CFG_DFS_1;
template <typename G_t, template<class G, class ...> class cfg>
struct CFG_DFS_2;
template <typename G_t, template<class G, class ...> class cfg>
struct CFG_DFS_3;

template <typename G_t, template<class G, class ...> class cfg>
struct CFG_DFS_p17;




/*
    -initial_lb_algo = deltaC_least_c
    -initial_ub_algo = NONE
    -lb_algo = NONE
    -next = all nodes "from left to right"
    -refiner = NONE
*/

template <typename G_t, template<class G, class ...> class cfg>
struct CFG_DFS_1 : generic_elimination_search_DFS<G_t, CFG_DFS_1<G_t, cfg>, cfg> {
    typedef generic_elimination_search_DFS<G_t, CFG_DFS_1<G_t, cfg>, cfg> baseclass;
    CFG_DFS_1(G_t const& G) : baseclass(G)
    {}

    typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;

    static unsigned INVALID_VERTEX()
    {
        return UINT_MAX;
    }

    static const std::string name()
    {
        return "CFG_DFS_1";
    }

    static unsigned initial_lb_algo(const G_t &G)
    {
        G_t H(G);
        return treedec::lb::deltaC_least_c(H)+1;
    }

    static unsigned initial_ub_algo(const G_t &G, std::vector<vd> &O)
    {
        for(unsigned i = 0; i < boost::num_vertices(G); ++i){
            O[i] = i;
        }
        return boost::num_vertices(G);
    }


    static unsigned lb_algo(const G_t &G){ //aka no lb algo
        return 0;
    }

    static vd next(const G_t & /*G*/, const std::vector<BOOL> &active, unsigned &idx)
    {
        for(; idx < active.size(); ++idx){
            if(active[idx]){
                return idx++;
            }
        }

        return INVALID_VERTEX();
    }

    static unsigned refiner(const G_t &G, std::vector<vd> &orig_elim, std::vector<vd> &new_elim) //aka no refiner
    {
        G_t H(G);
        treedec::minimalChordal(H, orig_elim, new_elim);
        G_t H2(G);
        return treedec::get_bagsize_of_elimination_ordering(H2, new_elim); //not necessary
    }

};

/*
    -initial_lb_algo = deltaC_least_c
    -initial_ub_algo = minDegree
    -lb_algo = NONE
    -next = all nodes "from left to right"
*/
template <typename G_t, template<class G, class ...> class CFGT>
struct CFG_DFS_2 : generic_elimination_search_DFS<G_t, CFG_DFS_2<G_t, CFGT>, CFGT> {
    typedef generic_elimination_search_DFS<G_t, CFG_DFS_2<G_t, CFGT>, CFGT> baseclass;
    CFG_DFS_2(G_t const& G) : baseclass(G)
    {}

    CFG_DFS_2(G_t const& G, unsigned m, unsigned n) : baseclass(G, m, n)
    {}

    typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;

    static unsigned INVALID_VERTEX()
    {
        return UINT_MAX;
    }

    static const std::string name()
    {
        return "CFG_DFS_2";
    }

    static unsigned initial_lb_algo(const G_t &G)
    {
        G_t H(G);
        return treedec::lb::deltaC_least_c(H)+1;
    }

    static unsigned initial_ub_algo(const G_t &G, std::vector<vd> &O)
    {
        G_t H(G);
        return treedec::minDegree_ordering(H, O)+1;
    }


    static unsigned lb_algo(const G_t &G){ //aka no lb algo
        return 0;
    }

    static vd next(const G_t & /*G*/, const std::vector<BOOL> &active, unsigned &idx)
    {
        for(; idx < active.size(); ++idx){
            if(active[idx]){
                return idx++;
            }
        }

        return INVALID_VERTEX();
    }

    static unsigned refiner(const G_t &G, std::vector<vd> &orig_elim, std::vector<vd> &new_elim)
    {
        G_t H(G);
        treedec::minimalChordal(H, orig_elim, new_elim);
        G_t H2(G);
        return treedec::get_bagsize_of_elimination_ordering(H2, new_elim); //not necessary
    }
};

/*
    -initial_lb_algo = deltaC_least_c
    -initial_ub_algo = fillIn
    -lb_algo = NONE
    -next = all nodes "from left to right"
*/
template <typename G_t, template<class G, class ...> class cfg>
struct CFG_DFS_3 : generic_elimination_search_DFS<G_t, CFG_DFS_3<G_t, cfg>, cfg> {
    typedef generic_elimination_search_DFS<G_t, CFG_DFS_3<G_t, cfg>, cfg> baseclass;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
    CFG_DFS_3(G_t const& G) : baseclass(G)
    {}

    static unsigned INVALID_VERTEX()
    {
        return UINT_MAX;
    }

    static const std::string name()
    {
        return "CFG_DFS_3";
    }

    static unsigned initial_lb_algo(const G_t &G)
    {
        G_t H(G);
        return treedec::lb::deltaC_least_c(H)+1;
    }

    static unsigned initial_ub_algo(const G_t &G, std::vector<vd> &O)
    {
        G_t H(G);
        return treedec::fillIn_ordering(H, O)+1;
    }


    static unsigned lb_algo(const G_t &G){ //aka no lb algo
        return 0;
    }

    static vd next(const G_t & /*G*/, const std::vector<BOOL> &active, unsigned &idx)
    {
        for(; idx < active.size(); ++idx){
            if(active[idx]){
                return idx++;
            }
        }

        return INVALID_VERTEX();
    }

    static unsigned refiner(const G_t &G, std::vector<vd> &orig_elim, std::vector<vd> &new_elim) //aka no refiner
    {
        G_t H(G);
        treedec::minimalChordal(H, orig_elim, new_elim);
        G_t H2(G);
        return treedec::get_bagsize_of_elimination_ordering(H2, new_elim); //TODO: not necessary
    }
};


/* PACE 2017 config
    -initial_lb_algo = deltaC_least_c
    -initial_ub_algo = minDegree
    -lb_algo = NONE
    -next = all nodes "from left to right"
*/
template <typename G_t, template<class G, class ...> class CFGT>
struct CFG_DFS_p17 : generic_elimination_search_DFS<G_t, CFG_DFS_p17<G_t, CFGT>, CFGT> {
    typedef generic_elimination_search_DFS<G_t, CFG_DFS_p17<G_t, CFGT>, CFGT> baseclass;
    CFG_DFS_p17(G_t const& G) : baseclass(G)
    {}

    CFG_DFS_p17(G_t const& G, unsigned m, unsigned n) : baseclass(G, m, n)
    {}

    typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;

    static unsigned INVALID_VERTEX()
    {
        return UINT_MAX;
    }

    static const std::string name()
    {
        return "CFG_DFS_p17";
    }

    static unsigned initial_lb_algo(const G_t &G)
    {
        G_t H(G);
        return treedec::lb::deltaC_least_c(H)+1;
    }

    static unsigned initial_ub_algo(const G_t &G, std::vector<vd> &O)
    {
        G_t H(G);
        return treedec::minDegree_ordering(H, O)+1;
    }


    static unsigned lb_algo(const G_t &G){ //aka no lb algo
        return 0;
    }

    static vd next(const G_t & /*G*/, const std::vector<BOOL> &active, unsigned &idx)
    {
        for(; idx < active.size(); ++idx){
            if(active[idx]){
                return idx++;
            }
        }

        return INVALID_VERTEX();
    }

    static unsigned refiner(const G_t &G, std::vector<vd> &orig_elim, std::vector<vd> &new_elim)
    {
        G_t H(G);
        treedec::minimalChordal(H, orig_elim, new_elim);
        G_t H2(G);
        return treedec::get_bagsize_of_elimination_ordering(H2, new_elim); //not necessary
    }
};

} //namespace configs

} //namespace gen_search

} //namespace treedec

#endif //TD_GENERIC_ELIM_SEARCH
// vim:ts=8:sw=4:et
