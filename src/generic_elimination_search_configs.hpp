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

#include <iostream>

namespace treedec{

namespace gen_search{

namespace configs{

/*
    -initial_lb_algo = deltaC_least_c
    -initial_ub_algo = NONE
    -lb_algo = NONE
    -next = all nodes "from left to right"
    -refiner = NONE
*/
template <typename G_t>
struct CFG_DFS_1{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;

    static const unsigned INVALID_VERTEX()
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

    static vd next(const G_t &G, const std::vector<bool> &active, unsigned &idx)
    {
        for(; idx < active.size(); ++idx){
            if(active[idx]){
                return idx++;
            }
        }

        return INVALID_VERTEX();
    }

/*
    static unsigned refiner(const G_t &G, std::vector<vd> &orig_elim, std::vector<vd> &new_elim) //aka no refiner
    {
        return boost::num_vertices(G);
    }
*/

    static unsigned refiner(const G_t &G, std::vector<vd> &orig_elim, std::vector<vd> &new_elim) //aka no refiner
    {
        G_t H(G);
        return treedec::minimalChordal(H, orig_elim, new_elim)+1;
    }

};

/*
    -initial_lb_algo = deltaC_least_c
    -initial_ub_algo = minDegree
    -lb_algo = NONE
    -next = all nodes "from left to right"
*/
template <typename G_t>
struct CFG_DFS_2{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;

    static const unsigned INVALID_VERTEX()
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

    static vd next(const G_t &G, const std::vector<bool> &active, unsigned &idx)
    {
        for(; idx < active.size(); ++idx){
            if(active[idx]){
                return idx++;
            }
        }

        return INVALID_VERTEX();
    }

/*
    {
        return boost::num_vertices(G);
    }
*/

    static unsigned refiner(const G_t &G, std::vector<vd> &orig_elim, std::vector<vd> &new_elim)
    {
        G_t H(G);
        return treedec::minimalChordal(H, orig_elim, new_elim)+1;
    }
};

/*
    -initial_lb_algo = deltaC_least_c
    -initial_ub_algo = fillIn
    -lb_algo = NONE
    -next = all nodes "from left to right"
*/
template <typename G_t>
struct CFG_DFS_3{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;

    static const unsigned INVALID_VERTEX()
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

    static vd next(const G_t &G, const std::vector<bool> &active, unsigned &idx)
    {
        for(; idx < active.size(); ++idx){
            if(active[idx]){
                return idx++;
            }
        }

        return INVALID_VERTEX();
    }

/*
    static unsigned refiner(const G_t &G, std::vector<vd> &orig_elim, std::vector<vd> &new_elim) //aka no refiner
    {
        return boost::num_vertices(G);
    }
*/

    static unsigned refiner(const G_t &G, std::vector<vd> &orig_elim, std::vector<vd> &new_elim) //aka no refiner
    {
        G_t H(G);
        return treedec::minimalChordal(H, orig_elim, new_elim)+1;
    }
};



/* AKA minDegree
    -initial_lb_algo = NONE
    -initial_ub_algo = NONE
    -lb_algo = NONE
    -ub_algo = NONE
    -next = "minDegree"
*/
/* this is just an example, to not use this - very inefficient
template <typename G_t>
struct CFG_DFS_4{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;

    static const unsigned INVALID_VERTEX()
    {
        return UINT_MAX;
    }

    static const std::string name()
    {
        return "CFG_DFS_4";
    }

    static unsigned initial_lb_algo(const G_t &G)
    {
        return 0;
    }

    static unsigned initial_ub_algo(const G_t &G, std::vector<vd> &O)
    {
        return boost::num_vertices(G);
    }


    static unsigned lb_algo(const G_t &G){ //aka no lb algo
        return 0;
    }

    static vd next(const G_t &G, const std::vector<bool> &active, unsigned &idx)
    {
        unsigned min = UINT_MAX;
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            if(!active[*vIt]){
                continue;
            }

            unsigned deg = 0;

            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, G); nIt != nEnd; nIt++){
                if(!active[*nIt]){
                    continue;
                }
                ++deg;
            }

            min = (min < deg)? min : deg;
        }

        for(; idx < active.size(); ++idx){
            if(!active[idx]){
                continue;
            }

            unsigned deg = 0;

            typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
            for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(idx, G); nIt != nEnd; nIt++){
                if(!active[*nIt]){
                    continue;
                }
                ++deg;
            }

            if(deg == min){
                return idx++;
            }
        }

        return INVALID_VERTEX();
    }

    static unsigned refiner(const G_t &G, std::vector<vd> &orig_elim, std::vector<vd> &new_elim) //aka no refiner
    {
        return boost::num_vertices(G);
    }
};
*/


} //namespace configs

} //namespace gen_search

} //namespace treedec

#endif //TD_GENERIC_ELIM_SEARCH
