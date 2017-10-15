
// graph iface gala specialisations
#pragma once

#include "graph_traits.hpp"
#include <gala/cbset.h>
#include <gala/graph.h>

namespace treedec {

		// incomplete: assuming cbs. fix later.
		// incomplete: assuming cbs. fix later.
		// incomplete: assuming cbs. fix later.
		// incomplete: assuming cbs. fix later.
		// incomplete: assuming cbs. fix later.
template< template<class T, typename... > class ECT, \
          template<class T, typename... > class VCT, \
          class VDP, \
          template<class GG> class CFG >
struct graph_helper<gala::graph<ECT, VCT, VDP, CFG> > {
	using G=gala::graph<ECT, VCT, VDP, CFG>;

	template<class S>
	static void close_neighbourhood(S& c, G const& g)
	{ itested();
		S verts(c);
		for(auto const&v : verts){
			c.merge(g.out_edges(v));
		}
	}

	template<class S>
	static void open_neighbourhood(S& c, G const& g)
	{
		S verts(c);
		close_neighbourhood(c, g);
		c.carve(verts);
		trim(c); // hack.
	}

	template<class S>
	static void saturate(S& c, G const& g) {
		// T cnb = closedNeighbor(c);
		S cnb(c);
		close_neighbourhood(cnb, g);

		assert(cnb.howmany() == cnb.howmany());
		S onb = cbset::diff(cnb, c);
		for (auto const&v : onb) {
			if (cbset::isSubset(g.out_edges(v), cnb)) { itested();
				cbset::insert(c, v);
//				trace1("sat", c);
			}else{
			}
		}
	}
};

} // treedec
