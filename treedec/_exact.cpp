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

#include "boost_compat.h"

#include <boost/python.hpp>
#include <boost/python/list.hpp>
#include <vector>
#include <set>

#include "graph_py.hpp"
#include "exact_cutset.hpp"
#include "exact_ta.hpp"
#include "util_py.hpp"
#include "algo_py.hpp"
#include "tr_iodec.h"

namespace py = boost::python;
using treedec::exact_ta;
using treedec::draft::exact_decomposition;
using treedec::algo::default_config;

typedef treedec::draft::exact_cutset<_balsvu, default_config> _cutset_balsvu;
typedef treedec::draft::exact_cutset<_balsvd, default_config> _cutset_balsvd;
typedef treedec::draft::exact_cutset<_balvvu, default_config> _cutset_balvvu;
typedef treedec::draft::exact_cutset<_balvvd, default_config> _cutset_balvvd;

template<class G>
struct cfg_8 : public default_config<G> {
	static constexpr unsigned max_vertex_index=255u;
};

template<class G, class... rest>
  struct cfg32 : public default_config<G>{
      static constexpr unsigned max_vertex_index=31;
};

template<class G, class... rest>
  struct cfg64 : public default_config<G>{
      static constexpr unsigned max_vertex_index=63;
};

template<class G, class... rest>
  struct cfg128 : public default_config<G>{
      static constexpr unsigned max_vertex_index=127;
};

template<class G, class... rest>
  struct cfg192 : public default_config<G>{
      static constexpr unsigned max_vertex_index=191;
};

template<class G, class... rest>
  struct cfg256 : public default_config<G>{
      static constexpr unsigned max_vertex_index=255;
};

template<class G, class... rest>
  struct cfg512 : public default_config<G>{
      static constexpr unsigned max_vertex_index=511;
};

template<class G, class... rest>
  struct cfg1024 : public default_config<G>{
      static constexpr unsigned max_vertex_index=1023;
};




typedef exact_ta<_balvvu, cfg_8> _ta_balvvu;
typedef exact_ta<_balsvu, cfg_8> _ta_balsvu;

typedef exact_decomposition<_balvvu, default_config, exact_ta> _ex17_balvvu;
typedef exact_decomposition<_balsvu, default_config, exact_ta> _ex17_balsvu;


typedef TR<_balvvu, default_config> _tr_balvvu;
typedef TR<_balsvu, default_config> _tr_balsvu;

typedef TR<_balvvu, cfg32> _tr_32_balvvu;
typedef TR<_balsvu, cfg32> _tr_32_balsvu;

typedef TR<_balvvu, cfg64> _tr_64_balvvu;
typedef TR<_balsvu, cfg64> _tr_64_balsvu;

typedef TR<_balvvu, cfg128> _tr_128_balvvu;
typedef TR<_balsvu, cfg128> _tr_128_balsvu;

typedef TR<_balvvu, cfg192> _tr_192_balvvu;
typedef TR<_balsvu, cfg192> _tr_192_balsvu;

typedef TR<_balvvu, cfg256> _tr_256_balvvu;
typedef TR<_balsvu, cfg256> _tr_256_balsvu;

typedef TR<_balvvu, cfg512> _tr_512_balvvu;
typedef TR<_balsvu, cfg512> _tr_512_balsvu;

typedef TR<_balvvu, cfg1024> _tr_1024_balvvu;
typedef TR<_balsvu, cfg1024> _tr_1024_balsvu;


BOOST_PYTHON_MODULE(_exact)
{ itested();

//	to_python_converter<std::pair<int, bool>, pair2tuple<int, bool> >(); // hack
//	to_python_converter<std::pair<int, int>, pair2tuple<int, int> >();
//	to_python_converter<std::pair<long, long>, pair2tuple<long, long> >();
//	to_python_converter<std::pair<unsigned, unsigned>, pair2tuple<unsigned, unsigned> >();
//	to_python_converter<std::pair<unsigned long, unsigned long>, pair2tuple<unsigned long, unsigned long> >();

	COMMON_ALGO_IFACE(_cutset, _balvvu);
	// COMMON_ALGO_IFACE(_cutset, _balvvd); does not compile.
	COMMON_ALGO_IFACE(_cutset, _balsvu);
	// COMMON_ALGO_IFACE(_cutset, _balsvd); does not compile.

	COMMON_ALGO_IFACE(_ta, _balvvu);
	COMMON_ALGO_IFACE(_ta, _balsvu);

	COMMON_ALGO_IFACE(_ex17, _balvvu);
	COMMON_ALGO_IFACE(_ex17, _balsvu);

        COMMON_ALGO_IFACE(_tr_32, _balsvu);
        COMMON_ALGO_IFACE(_tr_32, _balvvu);

	COMMON_ALGO_IFACE(_tr_64, _balsvu);
        COMMON_ALGO_IFACE(_tr_64, _balvvu);

	COMMON_ALGO_IFACE(_tr_128, _balsvu);
        COMMON_ALGO_IFACE(_tr_128, _balvvu);

	COMMON_ALGO_IFACE(_tr_192, _balsvu);
        COMMON_ALGO_IFACE(_tr_192, _balvvu);

	COMMON_ALGO_IFACE(_tr_256, _balsvu);
        COMMON_ALGO_IFACE(_tr_256, _balvvu);

	COMMON_ALGO_IFACE(_tr_512, _balsvu);
        COMMON_ALGO_IFACE(_tr_512, _balvvu);

	COMMON_ALGO_IFACE(_tr_1024, _balsvu);
        COMMON_ALGO_IFACE(_tr_1024, _balvvu);
}
