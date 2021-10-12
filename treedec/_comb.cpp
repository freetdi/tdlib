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

// various algorithm combination wrappers

#include "boost_compat.h"
#include <boost/graph/graph_traits.hpp>
#include "graph_py.hpp"
#include "directed_view.hpp"
#include "induced_subgraph.hpp"
#include "elimination_orderings.hpp"
#include "marker.hpp"
#include "elimination_orderings.hpp"
#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include <vector>
#include <set>
#include "graph_traits.hpp"
#include "boost_graph.hpp"
#include "algo_py.hpp"
#include "combinations.hpp"

using treedec::pending::PP_FI_TM;
using treedec::pending::PP_FI;
using treedec::algo::default_config;

typedef PP_FI<_balsvu, default_config> _ppfi_balsvu;
typedef PP_FI<_balsvd, default_config> _ppfi_balsvd;
typedef PP_FI<_balvvu, default_config> _ppfi_balvvu;
typedef PP_FI<_balvvd, default_config> _ppfi_balvvd;

typedef PP_FI_TM<_balsvu, default_config> _ppfitm_balsvu;
typedef PP_FI_TM<_balsvd, default_config> _ppfitm_balsvd;
typedef PP_FI_TM<_balvvu, default_config> _ppfitm_balvvu;
typedef PP_FI_TM<_balvvd, default_config> _ppfitm_balvvd;

BOOST_PYTHON_MODULE(_comb)
{ itested();
	namespace py = boost::python;

//	to_python_converter<std::pair<int, bool>, pair2tuple<int, bool> >(); // hack
//	to_python_converter<std::pair<int, int>, pair2tuple<int, int> >();
//	to_python_converter<std::pair<long, long>, pair2tuple<long, long> >();
//	to_python_converter<std::pair<unsigned, unsigned>, pair2tuple<unsigned, unsigned> >();
//	to_python_converter<std::pair<unsigned long, unsigned long>, pair2tuple<unsigned long, unsigned long> >();

	COMMON_ALGO_IFACE(_ppfitm, _balvvu);
	COMMON_ALGO_IFACE(_ppfitm, _balvvd);
	COMMON_ALGO_IFACE(_ppfitm, _balsvu);
	COMMON_ALGO_IFACE(_ppfitm, _balsvd);

	COMMON_ALGO_IFACE(_ppfi, _balvvu);
	COMMON_ALGO_IFACE(_ppfi, _balvvd);
	COMMON_ALGO_IFACE(_ppfi, _balsvu);
	COMMON_ALGO_IFACE(_ppfi, _balsvd);
}
