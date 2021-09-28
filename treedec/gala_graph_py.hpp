// Felix Salfelder, 2021
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
// Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
//

#ifndef GRAPH_GALA_PYTHON_H
#define GRAPH_GALA_PYTHON_H

#include <treedec/treedec_traits.hpp>
#include <boost/graph/graph_traits.hpp>
#include <gala/graph.h>
#include <gala/boost.h>
#include <gala/td.h>

// undirected simple loopless graph
namespace detail {
template<class G>
struct uvv_config : gala::graph_cfg_default<G> {
    static constexpr bool is_directed=false;
    static constexpr bool force_simple=true;
    static constexpr bool force_loopless=true; // not used yet?
    // static constexpr bool force_symmetric=true; // meaninngless (undirected)
    // typedef tdDEGS<G> degs_type; // obsolete.
};
}

typedef gala::graph<std::vector, std::vector, uint16_t, detail::uvv_config> _gsgvvu16;
typedef gala::graph<std::vector, std::vector, uint32_t, detail::uvv_config> _gsgvvu32;
typedef gala::graph<std::vector, std::vector, uint64_t, detail::uvv_config> _gsgvvu64;

typedef boost::property<treedec::bag_t, std::vector<uint16_t> > uint16_bag_p;
typedef boost::property<treedec::bag_t, std::vector<uint32_t> > uint32_bag_p;
typedef boost::property<treedec::bag_t, std::vector<uint64_t> > uint64_bag_p;

// #define tree_directedness_ boost::bidirectionalS //  broken?
#define tree_directedness_ boost::undirectedS

typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              tree_directedness_, uint16_bag_p> _gsgvvu16_treedec;

typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              tree_directedness_, uint32_bag_p> _gsgvvu32_treedec;

typedef boost::adjacency_list<boost::vecS, boost::vecS,
                              tree_directedness_, uint64_bag_p> _gsgvvu64_treedec;

#undef tree_directedness_

// NB: this is not needed when using ordinary properties.
// TREEDEC_TREEDEC_BAG_TRAITS(_gsgvvu32_treedec, bag); [..]
#endif
