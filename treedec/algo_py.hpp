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
//

#ifndef TREEDEC_ALGO_PY_H
#define TREEDEC_ALGO_PY_H
#include "util_py.hpp"

template<class what, class G>
static boost::shared_ptr<what> algo_factory( G& g )
{
	return boost::shared_ptr<what>(new what(g));
}

template<class T, StringLiteral n>
static std::string make_string(T const&)
{
	return std::string(n.value);
}

#define defgt(A, G, H) \
		    .def("get_treedec", static_cast<void ( A ## G:: * )(H ## _treedec&)> \
					     ( &A ## G::get_tree_decomposition< H ## _treedec >) )

#define COMMON_ALGO_IFACE(A, G) \
	py::class_< A ## G, boost::noncopyable>(# A # G, py::no_init) \
			 .def("__init__", py::make_constructor(algo_factory<A ## G, G>)) \
			  defgt(A, G, G) \
			 .def("__repr__", &make_string<A ## G, # A # G > ) \
			 .def("bagsize", &A ## G::bagsize) \
			 .def("lower_bound_bagsize", &A ## G::lower_bound_bagsize) \
			 .def("do_it", static_cast<void ( A ## G:: *)()> \
				        (&A ## G ::do_it) )

#endif
