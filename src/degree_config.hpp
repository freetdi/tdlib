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

#ifndef TREEDEC_DEGREE_CONFIG_HPP
#define TREEDEC_DEGREE_CONFIG_HPP
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
        assert(C.begin() != C.end());
        return *C.begin();
    }
    template <typename C_t>
    static vd_type pick_and_erase(C_t &)
	 {
		 incomplete();
    }
};

}

}

// use this. cleanup later.
namespace treedec{

namespace degs{

template<class G_t>
struct default_config : misc::detail::deg_config<G_t> {};


template<class G_t>
struct mapped_config : treedec::degs::default_config<G_t> {
    typedef typename boost::graph_traits<G_t>::vertices_size_type vertices_size_type;
    typedef typename boost::graph_traits<G_t>::vertex_descriptor vertex_descriptor;
    typedef typename boost::property_map<G_t, boost::vertex_index_t>::type idmap_type;
    typedef boost::iterator_property_map<vertex_descriptor*,
        idmap_type, vertex_descriptor, vertex_descriptor&> degree_type;
};

} // degs

} // treedec

#endif
