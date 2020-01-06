// (c) 2016, 2017, 2020 Felix Salfelder
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option) any
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
// some example graphs
//
#ifndef EXAMPLES_GALA_GRAPHS_H
#define EXAMPLES_GALA_GRAPHS_H

#ifdef HAVE_STX_BTREE_SET_H
// not currently used anywhere.
#include <gala/examples/svbs.h>
#include <gala/examples/svbs_random.h>
typedef simplegraph_vector_bs svbs;
#endif

template<class T, class...>
using tdDEGS = misc::DEGS<T>;


// possibly oriented directed simple loopless graph
template<class G>
struct odsvv_config : gala::graph_cfg_default<G> {
    static constexpr bool is_directed=true;
    static constexpr bool force_simple=true;
    static constexpr bool force_loopless=true; // not used yet?
    static constexpr bool force_symmetric=false;
    static constexpr bool force_oriented=true; // not used, how?!
    typedef tdDEGS<G> degs_type; // obsolete.
};

// parsed raw data.
template<class G>
struct dpvv_config : uvv_config<G> {
    static constexpr bool force_simple=false;
    static constexpr bool force_symmetric=false;
    static constexpr bool is_directed=true;
};

typedef gala::graph<std::vector, std::vector, uint16_t, dpvv_config> sg_dpvv16;
typedef gala::graph<std::vector, std::vector, uint32_t, dpvv_config> sg_dpvv32;

typedef gala::graph<std::vector, std::vector, uint16_t, odsvv_config> sg_odsvv16;


//typedef gala::graph<std::set, std::vector, uint16_t> ssg16;
#endif
