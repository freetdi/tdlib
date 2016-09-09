// Felix Salfelder 2016
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
//
// Configure degree tracker.

#ifndef TD_DEGREE_CONFIG_H
#define TD_DEGREE_CONFIG_H
#include <boost/graph/graph_traits.hpp>
#include "trace.hpp"
#include <set>

namespace misc {
namespace detail {

//	TODO:: amend namespace

template<class G_t>
struct deg_config{
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vd_type;

	 // use vector?!
    typedef std::set<vd_type> bag_type;
    static void alloc_init(size_t){
    }
    static unsigned num_threads(){return 1;}

    template <typename C_t>
    static vd_type pick(C_t const &C){
        return *C.begin();
    }
    template <typename C_t>
    static vd_type pick_and_erase(C_t &){
		 incomplete();
    }
};

}
}

#endif
