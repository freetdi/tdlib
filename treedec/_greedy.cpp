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

// various greedy algorithm wrappers

#include "boost_compat.h"
#include <boost/graph/graph_traits.hpp>
#include "graph_py.hpp"
#include "directed_view.hpp"
#include "induced_subgraph.hpp"
#include "elimination_orderings.hpp"
// treedec::impl::bmdo<sg_dvv16> A(g16, _elimord);

#include "marker.hpp"
#include "impl/fill_in.hpp"
#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include <vector>
#include <set>
#include "graph_traits.hpp"
#include "boost_graph.hpp"
#include "algo_py.hpp"

using treedec::impl::fillIn;
using treedec::algo::default_config;
using treedec::impl::bmdo;

typedef fillIn<_balsvu, default_config> _fi_balsvu;
typedef fillIn<_balsvd, default_config> _fi_balsvd;
typedef fillIn<_balvvu, default_config> _fi_balvvu;
typedef fillIn<_balvvd, default_config> _fi_balvvd;

typedef treedec::impl::bmdo<_balvvu> _bmd_balvvu;
typedef treedec::impl::bmdo<_balvvd> _bmd_balvvd;
typedef treedec::impl::bmdo<_balsvu> _bmd_balsvu;
typedef treedec::impl::bmdo<_balsvd> _bmd_balsvd;

#ifdef HAVE_GALA_GRAPH_H
typedef fillIn<_gsgvvu16, default_config> _fi_gsgvvu16;
typedef fillIn<_gsgvvu32, default_config> _fi_gsgvvu32;
typedef fillIn<_gsgvvu64, default_config> _fi_gsgvvu64;

typedef bmdo<_gsgvvu16> _bmd_gsgvvu16;
typedef bmdo<_gsgvvu32> _bmd_gsgvvu32;
typedef bmdo<_gsgvvu64> _bmd_gsgvvu64;
#endif

BOOST_PYTHON_MODULE(_greedy)
{ itested();
	namespace py = boost::python;

	COMMON_ALGO_IFACE(_fi, _balvvu);
	COMMON_ALGO_IFACE(_fi, _balvvd);
	COMMON_ALGO_IFACE(_fi, _balsvu);
	COMMON_ALGO_IFACE(_fi, _balsvd);

	COMMON_ALGO_IFACE(_bmd, _balvvu)
#ifdef HAVE_GALA_GRAPH_H
		    defgt(_bmd, _balvvu, _gsgvvu16)
		    defgt(_bmd, _balvvu, _gsgvvu32)
		    defgt(_bmd, _balvvu, _gsgvvu64)
#endif
	;
	COMMON_ALGO_IFACE(_bmd, _balvvd);
	COMMON_ALGO_IFACE(_bmd, _balsvu);
	COMMON_ALGO_IFACE(_bmd, _balsvd);

#ifdef HAVE_GALA_GRAPH_H
	COMMON_ALGO_IFACE(_fi, _gsgvvu16)
		    .def("get_treedec", &_fi_gsgvvu16::get_tree_decomposition<_balvvu_treedec>);
	COMMON_ALGO_IFACE(_fi, _gsgvvu32)
		    .def("get_treedec", &_fi_gsgvvu32::get_tree_decomposition<_balvvu_treedec>);
	COMMON_ALGO_IFACE(_fi, _gsgvvu64)
		    .def("get_treedec", &_fi_gsgvvu64::get_tree_decomposition<_balvvu_treedec>);

	COMMON_ALGO_IFACE(_bmd, _gsgvvu16)
		defgt(_bmd, _gsgvvu16, _balvvu)
		defgt(_bmd, _gsgvvu16, _gsgvvu32)
		defgt(_bmd, _gsgvvu16, _gsgvvu64);
	COMMON_ALGO_IFACE(_bmd, _gsgvvu32)
		defgt(_bmd, _gsgvvu32, _balvvu)
		defgt(_bmd, _gsgvvu32, _gsgvvu16)
		defgt(_bmd, _gsgvvu32, _gsgvvu64);
	COMMON_ALGO_IFACE(_bmd, _gsgvvu64)
		defgt(_bmd, _gsgvvu64, _balvvu)
		defgt(_bmd, _gsgvvu64, _gsgvvu16)
		defgt(_bmd, _gsgvvu64, _gsgvvu32);

#endif
}
